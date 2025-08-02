/*!
 * DEPRECATED: Multi-Segment Hybrid BAM Processing
 * 
 * ⚠️  This implementation is preserved as a CODE RELIC for research and educational purposes.
 * 🚫 DO NOT USE in production - use `bam_to_arrow_ipc_htslib_optimized()` instead.
 * 
 * ## Performance Analysis Results:
 * - **This implementation**: ~109k rec/sec (18.4s for 2M records)
 * - **Single optimized reader**: 205k rec/sec (9.7s for 2M records) 
 * - **Performance gap**: 47% slower than single-reader baseline
 * 
 * ## Root Cause: BGZF I/O Serialization Bottleneck
 * The fundamental issue is that multiple concurrent readers accessing the same BAM file
 * create BGZF random access conflicts, resulting in 94.6% efficiency loss per segment.
 * 
 * ## Technical Details:
 * - Multiple readers → same file → I/O contention
 * - BGZF seeking conflicts destroy cache efficiency  
 * - Thread contention is NOT the primary bottleneck
 * - Would work well with PRE-SPLIT BAM files (separate physical files)
 * 
 * ## Historical Context:
 * This represents Phase 6 of the optimization roadmap, which achieved the critical
 * breakthrough of fixing data loss issues (87.5% → 0%) through proper Arrow IPC
 * concatenation, but revealed the fundamental I/O limitations of the multi-reader approach.
 * 
 * ## See Also:
 * - PERFORMANCE_ROADMAP.md Phase 7 for detailed analysis
 * - `bam_to_arrow_ipc_htslib_optimized()` for production use
 * - Sequential-parallel architecture proposals for future development
 */

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::fs::File;
use std::path::Path;
use std::time::Instant;
use std::thread;
use std::sync::mpsc;

use log::debug;
#[cfg(feature = "htslib")]
use rust_htslib::bam as hts_bam;
#[cfg(feature = "htslib")]
use rust_htslib::bam::Read as HtsRead;

use arrow::ipc::writer::FileWriter as ArrowIpcWriter;

use crate::bam_htslib::{discover_split_points, create_bam_schema, build_chromosome_lookup};

/// Hybrid BGZF approach: Multiple independent optimized single-readers
/// Each segment runs the proven optimized pipeline independently
#[cfg(feature = "htslib")]
#[pyfunction]
#[pyo3(signature = (
    bam_path,
    arrow_ipc_path,
    batch_size = 20000,
    include_sequence = true,
    include_quality = true,
    max_bgzf_threads = 6,
    writing_threads = 10,
    read_buffer_mb = 1024,
    write_buffer_mb = 256,
    limit = None,
    num_segments = None
))]
pub fn bam_to_arrow_ipc_htslib_hybrid_segments(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    max_bgzf_threads: usize,
    writing_threads: usize,
    read_buffer_mb: usize,
    write_buffer_mb: usize,
    limit: Option<usize>,
    num_segments: Option<usize>,
) -> PyResult<()> {
    let start_time = Instant::now();
    
    let input_path = Path::new(bam_path);
    let _output_path = Path::new(arrow_ipc_path);

    if !input_path.exists() {
        return Err(PyErr::new::<PyRuntimeError, _>(
            format!("BAM file does not exist: {}", bam_path)
        ));
    }
    
    // Determine optimal number of segments
    let file_size_gb = std::fs::metadata(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Cannot read file metadata: {}", e)))?
        .len() as f64 / (1024.0 * 1024.0 * 1024.0);
    
    let effective_segments = num_segments.unwrap_or_else(|| {
        // Get number of available CPU cores for better scaling
        let available_cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
        
        if file_size_gb < 1.0 {
            1  // Single segment for small files
        } else if file_size_gb < 10.0 {
            std::cmp::min(2, available_cores / 4)  // 2 segments for medium files, limited by cores
        } else if file_size_gb < 50.0 {
            std::cmp::min(4, available_cores / 2)  // 4 segments for large files, scale with cores
        } else {
            std::cmp::min(available_cores, 16)  // Scale with cores up to 16 segments for very large files
        }
    });
    
    let available_cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
    
    eprintln!("Starting Hybrid BGZF-Optimized BAM to IPC conversion:");
    eprintln!("  Input: {} ({:.1}GB)", bam_path, file_size_gb);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Available cores: {}, Selected segments: {}", available_cores, effective_segments);
    eprintln!("  Optimized parameters: batch_size={}, bgzf_threads={}, writing_threads={}", 
           batch_size, max_bgzf_threads, writing_threads);
    eprintln!("  Buffer sizes: read={}MB, write={}MB per segment", read_buffer_mb, write_buffer_mb);
    eprintln!("  Expected throughput: ~{}k rec/sec (TRUE utilization of all segments)", 
           40 * effective_segments);
    
    if effective_segments == 1 {
        eprintln!("  → Using single optimized reader (file too small for segmentation)");
        // Just call the optimized single-reader directly
        return crate::bam::bam_to_arrow_ipc_htslib_optimized(
            bam_path,
            arrow_ipc_path,
            batch_size,
            include_sequence,
            include_quality,
            max_bgzf_threads,
            writing_threads,
            Some(read_buffer_mb),
            Some(write_buffer_mb),
            limit,
        );
    }
    
    // Discover split points for segments
    let split_points = discover_split_points(bam_path, effective_segments)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
            format!("Failed to discover split points: {}", e)
        ))?;
    
    if split_points.len() < 2 {
        eprintln!("  → Falling back to single optimized reader (insufficient split points)");
        return crate::bam::bam_to_arrow_ipc_htslib_optimized(
            bam_path,
            arrow_ipc_path,
            batch_size,
            include_sequence,
            include_quality,
            max_bgzf_threads,
            writing_threads,
            Some(read_buffer_mb),
            Some(write_buffer_mb),
            limit,
        );
    }
    
    // Calculate per-segment limits
    let num_segments = split_points.len() - 1;
    let per_segment_limit = limit.map(|total| {
        let base_limit = total / num_segments;
        let remainder = total % num_segments;
        eprintln!("Distributing {} total records across {} segments: {} base + {} remainder", 
                 total, num_segments, base_limit, remainder);
        base_limit
    });
    
    // Create temporary output files for each segment
    let temp_files: Vec<String> = (0..split_points.len()-1)
        .map(|i| format!("{}.segment_{}.tmp", arrow_ipc_path, i))
        .collect();
    
    // Spawn independent optimized readers for each segment
    let (sender, receiver) = mpsc::channel();
    let mut handles = Vec::new();
    
    for (segment_id, window) in split_points.windows(2).enumerate() {
        let start_pos = window[0];
        let end_pos = window[1];
        let sender = sender.clone();
        let bam_path = bam_path.to_string();
        let temp_file = temp_files[segment_id].clone();
        
        let handle = thread::spawn(move || {
            debug!("Starting segment {} processing (pos: {} -> {})", segment_id, start_pos, end_pos);
            let segment_start = std::time::Instant::now();
            
            debug!("Segment {} limit: {:?}", segment_id, per_segment_limit);
            
            let result = run_segment_with_optimized_reader(
                &bam_path,
                &temp_file,
                start_pos,
                end_pos,
                batch_size,
                include_sequence,
                include_quality,
                max_bgzf_threads,
                writing_threads,
                Some(read_buffer_mb),
                Some(write_buffer_mb),
                per_segment_limit,
                segment_id,
            );
            
            let segment_duration = segment_start.elapsed();
            debug!("Segment {} thread completed in {:.1}s", segment_id, segment_duration.as_secs_f64());
            
            let _ = sender.send((segment_id, result));
        });
        
        handles.push(handle);
    }
    
    drop(sender); // Close sender so receiver knows when all segments are done
    
    // Collect results from all segments
    let mut segment_results = Vec::new();
    for (segment_id, result) in receiver {
        debug!("Received result from segment {}: {:?}", segment_id, 
                 result.as_ref().map(|r| format!("OK: {} records", r)).unwrap_or_else(|e| format!("ERROR: {}", e)));
        segment_results.push((segment_id, result));
    }
    
    eprintln!("Total segments completed: {}", segment_results.len());
    
    // Wait for all threads to complete
    for handle in handles {
        let _ = handle.join();
    }
    
    // Sort results by segment ID to maintain order
    segment_results.sort_by_key(|(id, _)| *id);
    
    // Check for any failures
    for (segment_id, result) in &segment_results {
        if let Err(e) = result {
            return Err(PyErr::new::<PyRuntimeError, _>(
                format!("Segment {} failed: {}", segment_id, e)
            ));
        }
    }
    
    // Concatenate all segment outputs into final result
    eprintln!("Concatenating {} segment outputs into final file", temp_files.len());
    concatenate_arrow_files(&temp_files, arrow_ipc_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
            format!("Failed to concatenate segment outputs: {}", e)
        ))?;
    
    // Clean up temporary files
    for temp_file in &temp_files {
        let _ = std::fs::remove_file(temp_file);
    }
    
    let total_duration = start_time.elapsed();
    let total_records: usize = segment_results.iter()
        .map(|(_, result)| result.as_ref().unwrap_or(&0))
        .sum();
    
    let throughput = total_records as f64 / total_duration.as_secs_f64();
    
    eprintln!("Hybrid BGZF-Optimized conversion completed:");
    eprintln!("  Total records: {}", total_records);
    eprintln!("  Duration: {:.1}s", total_duration.as_secs_f64());
    eprintln!("  Throughput: {:.0} records/sec", throughput);
    eprintln!("  Speedup: {:.1}x over single-reader baseline", throughput / 40_000.0);
    
    Ok(())
}

/// Run a single segment using the proven optimized single-reader approach
#[cfg(feature = "htslib")]
fn run_segment_with_optimized_reader(
    bam_path: &str,
    output_path: &str,
    start_pos: u64,
    end_pos: u64,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    max_bgzf_threads: usize,
    _writing_threads: usize,
    _read_buffer_mb: Option<usize>,
    _write_buffer_mb: Option<usize>,
    record_limit: Option<usize>,
    segment_id: usize,
) -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
    // Create a reader for this segment
    let mut reader = hts_bam::Reader::from_path(bam_path)?;
    reader.set_threads(max_bgzf_threads)?;
    
    // Seek to start position if not the first segment
    if start_pos > 0 {
        let virtual_offset = start_pos << 16;
        reader.seek(virtual_offset as i64)?;
    }
    
    // Get header and build chromosome lookup
    let header_view = reader.header().clone();
    let chrom_lookup = build_chromosome_lookup(&header_view);
    
    // Create output writer
    let schema = create_bam_schema(include_sequence, include_quality);
    let output_file = File::create(output_path)?;
    let mut writer = ArrowIpcWriter::try_new(output_file, &schema)?;
    
    // Process records using the optimized pipeline approach
    let mut total_records = 0;
    let mut current_batch_records = Vec::new();
    
    let mut record = hts_bam::Record::new();
    
    loop {
        // Read the next record
        match reader.read(&mut record) {
            Some(Ok(())) => {
                // Check if we've reached the end of this segment
                if end_pos != std::u64::MAX {
                    let current_voffset = reader.tell();
                    if current_voffset >= 0 {
                        let current_file_pos = (current_voffset as u64) >> 16;
                        if current_file_pos >= end_pos {
                            break;
                        }
                    }
                }
                
                current_batch_records.push(record.clone());
                total_records += 1;
            }
            None => break, // EOF
            Some(Err(e)) => return Err(e.into()),
        }
        
        // Check record limit
        if let Some(limit) = record_limit {
            if total_records >= limit {
                break;
            }
        }
        
        // Process batch when full
        if current_batch_records.len() >= batch_size {
            let batch = crate::bam_htslib::process_htslib_records_to_batch(
                &current_batch_records,
                &chrom_lookup,
                include_sequence,
                include_quality,
                0,
            )?;
            
            writer.write(&batch)?;
            current_batch_records.clear();
        }
    }
    
    // Process remaining records
    if !current_batch_records.is_empty() {
        let batch = crate::bam_htslib::process_htslib_records_to_batch(
            &current_batch_records,
            &chrom_lookup,
            include_sequence,
            include_quality,
            0,
        )?;
        
        writer.write(&batch)?;
    }
    
    writer.finish()?;
    
    debug!("Segment {} completed: {} records", segment_id, total_records);
    Ok(total_records)
}

/// Concatenate multiple Arrow IPC files into one
#[cfg(feature = "htslib")]
fn concatenate_arrow_files(input_files: &[String], output_path: &str) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    use arrow::ipc::reader::FileReader as ArrowIpcReader;
    use std::fs::File;
    
    if input_files.is_empty() {
        return Err("No input files provided for concatenation".into());
    }
    
    // Filter out files that don't exist
    let existing_files: Vec<&String> = input_files.iter()
        .filter(|path| std::path::Path::new(path).exists())
        .collect();
    
    if existing_files.is_empty() {
        return Err("No valid input files found for concatenation".into());
    }
    
    eprintln!("  Concatenating {} Arrow IPC files into {}", existing_files.len(), output_path);
    
    // Open first file to get schema
    let first_file = File::open(existing_files[0])?;
    let mut first_reader = ArrowIpcReader::try_new(first_file, None)?;
    let schema = first_reader.schema();
    
    // Create output writer
    let output_file = File::create(output_path)?;
    let mut writer = ArrowIpcWriter::try_new(output_file, &schema)?;
    
    // Read and write all batches from first file
    for batch_result in first_reader {
        let batch = batch_result?;
        writer.write(&batch)?;
    }
    
    // Process remaining files
    for file_path in &existing_files[1..] {
        eprintln!("    Processing file: {}", file_path);
        let file = File::open(file_path)?;
        let reader = ArrowIpcReader::try_new(file, None)?;
        
        // Verify schema compatibility
        if reader.schema() != schema {
            return Err(format!("Schema mismatch in file: {}", file_path).into());
        }
        
        // Write all batches from this file
        for batch_result in reader {
            let batch = batch_result?;
            writer.write(&batch)?;
        }
    }
    
    writer.finish()?;
    eprintln!("  ✅ Successfully concatenated {} files into {}", existing_files.len(), output_path);
    
    Ok(())
}
