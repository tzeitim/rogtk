/*!
 * DEPRECATED: Hybrid Minimal Fix BAM Processing  
 * 
 * ⚠️  CODE RELIC - Preserved for research purposes only
 * 🚫 Use `bam_to_arrow_ipc_htslib_optimized()` for production (205k+ rec/sec)
 * 
 * This was a minimal fix attempt for the hybrid approach, but the fundamental
 * BGZF I/O bottleneck cannot be resolved through parameter tuning alone.
 * See PERFORMANCE_ROADMAP.md Phase 7 for detailed analysis.
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
use arrow::ipc::reader::FileReader as ArrowIpcReader;

use crate::bam_htslib::{discover_split_points, create_bam_schema, build_chromosome_lookup};

/// Minimal optimization: Fix segment concatenation ONLY
/// This is the original hybrid implementation with ONE critical fix:
/// Use ALL segment outputs instead of just the first one
#[cfg(feature = "htslib")]
#[pyfunction]
#[pyo3(signature = (
    bam_path,
    arrow_ipc_path,
    batch_size = 15000,
    include_sequence = true,
    include_quality = true,
    max_bgzf_threads = 4,
    writing_threads = 6,
    read_buffer_mb = Some(2048),
    write_buffer_mb = Some(128),
    limit = None,
    num_segments = None
))]
pub fn bam_to_arrow_ipc_htslib_hybrid_minimal_fix(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    max_bgzf_threads: usize,
    writing_threads: usize,
    read_buffer_mb: Option<usize>,
    write_buffer_mb: Option<usize>,
    limit: Option<usize>,
    num_segments: Option<usize>,
) -> PyResult<()> {
    let start_time = Instant::now();
    
    let input_path = Path::new(bam_path);

    if !input_path.exists() {
        return Err(PyErr::new::<PyRuntimeError, _>(
            format!("BAM file does not exist: {}", bam_path)
        ));
    }
    
    // Determine optimal number of segments (same as original)
    let file_size_gb = std::fs::metadata(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Cannot read file metadata: {}", e)))?
        .len() as f64 / (1024.0 * 1024.0 * 1024.0);
    
    let effective_segments = num_segments.unwrap_or_else(|| {
        if file_size_gb < 1.0 {
            1
        } else if file_size_gb < 10.0 {
            2
        } else if file_size_gb < 50.0 {
            4
        } else {
            8
        }
    });
    
    debug!("Starting Minimal-Fix Hybrid BAM to IPC conversion:");
    debug!("  Input: {} ({:.1}GB)", bam_path, file_size_gb);
    debug!("  Output: {}", arrow_ipc_path);
    debug!("  Segments: {} (with PROPER concatenation)", effective_segments);
    
    if effective_segments == 1 {
        return crate::bam::bam_to_arrow_ipc_htslib_optimized(
            bam_path,
            arrow_ipc_path,
            batch_size,
            include_sequence,
            include_quality,
            max_bgzf_threads,
            writing_threads,
            read_buffer_mb,
            write_buffer_mb,
            limit,
        );
    }
    
    // Discover split points (same as original)
    let split_points = discover_split_points(bam_path, effective_segments)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
            format!("Failed to discover split points: {}", e)
        ))?;
    
    if split_points.len() < 2 {
        return crate::bam::bam_to_arrow_ipc_htslib_optimized(
            bam_path,
            arrow_ipc_path,
            batch_size,
            include_sequence,
            include_quality,
            max_bgzf_threads,
            writing_threads,
            read_buffer_mb,
            write_buffer_mb,
            limit,
        );
    }
    
    // Calculate per-segment limits (same as original)
    let num_segments = split_points.len() - 1;
    let per_segment_limit = limit.map(|total| total / num_segments);
    
    // Create temporary output files for each segment
    let temp_files: Vec<String> = (0..num_segments)
        .map(|i| format!("{}.segment_{}.tmp", arrow_ipc_path, i))
        .collect();
    
    // Spawn segments (identical to original)
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
            
            let result = run_original_segment(
                &bam_path,
                &temp_file,
                start_pos,
                end_pos,
                batch_size,
                include_sequence,
                include_quality,
                max_bgzf_threads,
                writing_threads,
                read_buffer_mb,
                write_buffer_mb,
                per_segment_limit,
                segment_id,
            );
            
            let segment_duration = segment_start.elapsed();
            debug!("Segment {} completed in {:.1}s", segment_id, segment_duration.as_secs_f64());
            
            let _ = sender.send((segment_id, result));
        });
        
        handles.push(handle);
    }
    
    drop(sender);
    
    // Collect results (same as original)
    let mut segment_results = Vec::new();
    let mut total_records = 0;
    
    for (segment_id, result) in receiver {
        match &result {
            Ok(count) => {
                total_records += count;
                debug!("Segment {} completed: {} records", segment_id, count);
            }
            Err(e) => {
                debug!("Segment {} failed: {}", segment_id, e);
            }
        }
        segment_results.push((segment_id, result));
    }
    
    // Wait for all threads
    for handle in handles {
        let _ = handle.join();
    }
    
    // Sort results by segment ID
    segment_results.sort_by_key(|(id, _)| *id);
    
    // Check for failures
    for (segment_id, result) in &segment_results {
        if let Err(e) = result {
            return Err(PyErr::new::<PyRuntimeError, _>(
                format!("Segment {} failed: {}", segment_id, e)
            ));
        }
    }
    
    // THE CRITICAL FIX: Properly concatenate ALL segment files instead of using just the first
    concatenate_arrow_files_simple(&temp_files, arrow_ipc_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
            format!("Failed to concatenate segments: {}", e)
        ))?;
    
    // Clean up temporary files
    for temp_file in &temp_files {
        let _ = std::fs::remove_file(temp_file);
    }
    
    let total_duration = start_time.elapsed();
    let throughput = total_records as f64 / total_duration.as_secs_f64();
    let efficiency = (throughput / (153_000.0 * effective_segments as f64)) * 100.0;
    
    debug!("Minimal-Fix Hybrid conversion completed:");
    debug!("  Total records: {}", total_records);
    debug!("  Duration: {:.1}s", total_duration.as_secs_f64());
    debug!("  Throughput: {:.0} records/sec", throughput);
    debug!("  Efficiency: {:.1}% per segment", efficiency);
    debug!("  FIXED: Used ALL {} segments instead of just the first", num_segments);
    
    Ok(())
}

/// Run original segment processing (identical to original implementation)
#[cfg(feature = "htslib")]
fn run_original_segment(
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
    // Identical to original segment processing
    let mut reader = hts_bam::Reader::from_path(bam_path)?;
    reader.set_threads(max_bgzf_threads)?;
    
    if start_pos > 0 {
        let virtual_offset = start_pos << 16;
        reader.seek(virtual_offset as i64)?;
    }
    
    let header_view = reader.header().clone();
    let chrom_lookup = build_chromosome_lookup(&header_view);
    
    let schema = create_bam_schema(include_sequence, include_quality);
    let output_file = File::create(output_path)?;
    let mut writer = ArrowIpcWriter::try_new(output_file, &schema)?;
    
    let mut total_records = 0;
    let mut current_batch_records = Vec::new();
    let mut record = hts_bam::Record::new();
    
    loop {
        match reader.read(&mut record) {
            Some(Ok(())) => {
                // Check segment boundary
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
            None => break,
            Some(Err(e)) => return Err(e.into()),
        }
        
        if let Some(limit) = record_limit {
            if total_records >= limit {
                break;
            }
        }
        
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
    
    debug!("Original segment {} completed: {} records", segment_id, total_records);
    Ok(total_records)
}

/// Simple and fast Arrow file concatenation
#[cfg(feature = "htslib")]
fn concatenate_arrow_files_simple(input_files: &[String], output_path: &str) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    if input_files.is_empty() {
        return Err("No input files to concatenate".into());
    }
    
    // Filter out non-existent files
    let existing_files: Vec<&String> = input_files.iter().filter(|f| Path::new(f).exists()).collect();
    
    if existing_files.is_empty() {
        return Err("No valid input files found".into());
    }
    
    // If only one file, just copy it
    if existing_files.len() == 1 {
        std::fs::copy(existing_files[0], output_path)?;
        return Ok(());
    }
    
    // Open first file to get schema
    let first_file = File::open(existing_files[0])?;
    let first_reader = ArrowIpcReader::try_new(first_file, None)?;
    let schema = first_reader.schema();
    
    // Create output writer
    let output_file = File::create(output_path)?;
    let mut writer = ArrowIpcWriter::try_new(output_file, &schema)?;
    
    // Process first file
    for batch_result in first_reader {
        let batch = batch_result?;
        writer.write(&batch)?;
    }
    
    // Process remaining files
    for input_file in &existing_files[1..] {
        let file = File::open(input_file)?;
        let reader = ArrowIpcReader::try_new(file, None)?;
        
        for batch_result in reader {
            let batch = batch_result?;
            writer.write(&batch)?;
        }
    }
    
    writer.finish()?;
    Ok(())
}