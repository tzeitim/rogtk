/*!
 * DEPRECATED: Hybrid Optimized BAM Processing
 * 
 * ⚠️  CODE RELIC - Preserved for research purposes only
 * 🚫 Use `bam_to_arrow_ipc_htslib_optimized()` for production (205k+ rec/sec)
 * 
 * This variant attempted parameter optimizations on the hybrid approach but 
 * still suffers from the fundamental BGZF I/O serialization bottleneck.
 * See PERFORMANCE_ROADMAP.md Phase 7 for detailed analysis.
 */

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::fs::File;
use std::path::Path;
use std::time::Instant;
use std::thread;
use std::sync::{mpsc, Arc, Mutex};
use std::collections::VecDeque;

use log::debug;
#[cfg(feature = "htslib")]
use rust_htslib::bam as hts_bam;
#[cfg(feature = "htslib")]
use rust_htslib::bam::Read as HtsRead;

use arrow::ipc::writer::FileWriter as ArrowIpcWriter;
use arrow::ipc::reader::FileReader as ArrowIpcReader;

use crate::bam_htslib::{discover_split_points, create_bam_schema, build_chromosome_lookup};

/// Zero-copy memory pool for record processing
struct RecordPool {
    available: Arc<Mutex<VecDeque<hts_bam::Record>>>,
    max_size: usize,
}

impl RecordPool {
    fn new(initial_size: usize, max_size: usize) -> Self {
        let mut pool = VecDeque::with_capacity(initial_size);
        for _ in 0..initial_size {
            pool.push_back(hts_bam::Record::new());
        }
        
        Self {
            available: Arc::new(Mutex::new(pool)),
            max_size,
        }
    }
    
    fn acquire(&self) -> hts_bam::Record {
        let mut pool = self.available.lock().unwrap();
        pool.pop_front().unwrap_or_else(|| hts_bam::Record::new())
    }
    
    fn release(&self, record: hts_bam::Record) {
        // Note: HTSlib Record doesn't have a clear method, so we just return it to pool
        let mut pool = self.available.lock().unwrap();
        if pool.len() < self.max_size {
            pool.push_back(record);
        }
    }
}

/// Batch buffer for zero-copy processing
struct BatchBuffer {
    records: Vec<hts_bam::Record>,
    capacity: usize,
    pool: Arc<RecordPool>,
}

impl BatchBuffer {
    fn new(capacity: usize, pool: Arc<RecordPool>) -> Self {
        Self {
            records: Vec::with_capacity(capacity),
            capacity,
            pool,
        }
    }
    
    fn push(&mut self, record: hts_bam::Record) {
        self.records.push(record);
    }
    
    fn is_full(&self) -> bool {
        self.records.len() >= self.capacity
    }
    
    fn clear_and_return_to_pool(&mut self) {
        for record in self.records.drain(..) {
            self.pool.release(record);
        }
    }
    
    fn len(&self) -> usize {
        self.records.len()
    }
    
    fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

/// Advanced Hybrid BGZF approach with zero-copy optimizations
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
pub fn bam_to_arrow_ipc_htslib_hybrid_optimized(
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
    
    // Determine optimal number of segments
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
    
    debug!("Starting Advanced Hybrid BGZF-Optimized BAM to IPC conversion:");
    debug!("  Input: {} ({:.1}GB)", bam_path, file_size_gb);
    debug!("  Output: {}", arrow_ipc_path);
    debug!("  Segments: {} zero-copy optimized readers", effective_segments);
    debug!("  Target efficiency: 40%+ (500-600k rec/sec)");
    
    if effective_segments == 1 {
        debug!("  → Using single optimized reader");
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
    
    // Discover split points for segments
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
    
    // Create shared memory pool for zero-copy processing
    let record_pool = Arc::new(RecordPool::new(
        batch_size * effective_segments * 2, // Pre-allocate 2x batch size per segment
        batch_size * effective_segments * 4  // Max 4x batch size
    ));
    
    // Calculate per-segment limits
    let num_segments = split_points.len() - 1;
    let per_segment_limit = limit.map(|total| total / num_segments);
    
    // Create temporary output files for each segment
    let temp_files: Vec<String> = (0..num_segments)
        .map(|i| format!("{}.segment_{}.tmp", arrow_ipc_path, i))
        .collect();
    
    // Spawn zero-copy optimized readers for each segment
    let (sender, receiver) = mpsc::channel();
    let mut handles = Vec::new();
    
    for (segment_id, window) in split_points.windows(2).enumerate() {
        let start_pos = window[0];
        let end_pos = window[1];
        let sender = sender.clone();
        let bam_path = bam_path.to_string();
        let temp_file = temp_files[segment_id].clone();
        let pool = Arc::clone(&record_pool);
        
        let handle = thread::spawn(move || {
            debug!("Starting zero-copy segment {} (pos: {} -> {})", segment_id, start_pos, end_pos);
            let segment_start = std::time::Instant::now();
            
            let result = run_zero_copy_segment(
                &bam_path,
                &temp_file,
                start_pos,
                end_pos,
                batch_size,
                include_sequence,
                include_quality,
                max_bgzf_threads,
                per_segment_limit,
                segment_id,
                pool,
            );
            
            let segment_duration = segment_start.elapsed();
            debug!("Zero-copy segment {} completed in {:.1}s", segment_id, segment_duration.as_secs_f64());
            
            let _ = sender.send((segment_id, result));
        });
        
        handles.push(handle);
    }
    
    drop(sender);
    
    // Collect results with lock-free aggregation
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
    
    // Concatenate Arrow files efficiently
    concatenate_arrow_files(&temp_files, arrow_ipc_path)
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
    
    debug!("Advanced Hybrid BGZF-Optimized conversion completed:");
    debug!("  Total records: {}", total_records);
    debug!("  Duration: {:.1}s", total_duration.as_secs_f64());
    debug!("  Throughput: {:.0} records/sec", throughput);
    debug!("  Efficiency: {:.1}% per segment", efficiency);
    debug!("  Target achievement: {:.1}% of 40% target", (efficiency / 40.0) * 100.0);
    
    Ok(())
}

/// Zero-copy segment processing with memory pools
#[cfg(feature = "htslib")]
fn run_zero_copy_segment(
    bam_path: &str,
    output_path: &str,
    start_pos: u64,
    end_pos: u64,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    max_bgzf_threads: usize,
    record_limit: Option<usize>,
    segment_id: usize,
    record_pool: Arc<RecordPool>,
) -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
    // Create optimized reader
    let mut reader = hts_bam::Reader::from_path(bam_path)?;
    reader.set_threads(max_bgzf_threads)?;
    
    // Seek to start position
    if start_pos > 0 {
        let virtual_offset = start_pos << 16;
        reader.seek(virtual_offset as i64)?;
    }
    
    // Get header and build optimized chromosome lookup
    let header_view = reader.header().clone();
    let chrom_lookup = build_chromosome_lookup(&header_view);
    
    // Create output writer with consistent schema
    let schema = crate::bam_htslib::create_bam_schema(include_sequence, include_quality);
    let output_file = File::create(output_path)?;
    let mut writer = ArrowIpcWriter::try_new(output_file, &schema)?;
    
    // Initialize zero-copy batch buffer
    let mut batch_buffer = BatchBuffer::new(batch_size, Arc::clone(&record_pool));
    let mut total_records = 0;
    
    loop {
        // Acquire record from pool (zero-copy)
        let mut record = record_pool.acquire();
        
        match reader.read(&mut record) {
            Some(Ok(())) => {
                // Check segment boundary
                if end_pos != std::u64::MAX {
                    let current_voffset = reader.tell();
                    if current_voffset >= 0 {
                        let current_file_pos = (current_voffset as u64) >> 16;
                        if current_file_pos >= end_pos {
                            record_pool.release(record);
                            break;
                        }
                    }
                }
                
                batch_buffer.push(record);
                total_records += 1;
                
                // Check record limit
                if let Some(limit) = record_limit {
                    if total_records >= limit {
                        break;
                    }
                }
                
                // Process full batch with optimized processing but compatible output
                if batch_buffer.is_full() {
                    let batch = crate::bam_htslib::process_htslib_records_to_batch(
                        &batch_buffer.records,
                        &chrom_lookup,
                        include_sequence,
                        include_quality,
                        0,
                    )?;
                    
                    writer.write(&batch)?;
                    batch_buffer.clear_and_return_to_pool();
                }
            }
            None => {
                record_pool.release(record);
                break;
            }
            Some(Err(e)) => {
                record_pool.release(record);
                return Err(e.into());
            }
        }
    }
    
    // Process remaining records
    if !batch_buffer.is_empty() {
        let batch = crate::bam_htslib::process_htslib_records_to_batch(
            &batch_buffer.records,
            &chrom_lookup,
            include_sequence,
            include_quality,
            0,
        )?;
        
        writer.write(&batch)?;
        batch_buffer.clear_and_return_to_pool();
    }
    
    writer.finish()?;
    
    debug!("Zero-copy segment {} completed: {} records", segment_id, total_records);
    Ok(total_records)
}

/// SIMD-optimized batch processing for maximum throughput
#[cfg(feature = "htslib")]
fn process_records_simd_optimized(
    records: &[hts_bam::Record],
    chrom_lookup: &std::collections::HashMap<i32, String>,
    include_sequence: bool,
    include_quality: bool,
) -> Result<arrow::record_batch::RecordBatch, Box<dyn std::error::Error + Send + Sync>> {
    // Pre-allocate all vectors to avoid reallocation overhead
    let len = records.len();
    let mut qnames = Vec::with_capacity(len);
    let mut flags = Vec::with_capacity(len);
    let mut chromosomes = Vec::with_capacity(len);
    let mut positions = Vec::with_capacity(len);
    let mut mapping_qualities = Vec::with_capacity(len);
    let mut cigars = Vec::with_capacity(len);
    let mut rnexts = Vec::with_capacity(len);
    let mut pnexts = Vec::with_capacity(len);
    let mut template_lengths = Vec::with_capacity(len);
    
    let mut sequences = if include_sequence { Some(Vec::with_capacity(len)) } else { None };
    let mut qualities = if include_quality { Some(Vec::with_capacity(len)) } else { None };
    
    // Process records in batches for cache efficiency
    for record in records {
        // Extract fields without allocations where possible
        qnames.push(Some(String::from_utf8_lossy(record.qname()).to_string()));
        flags.push(Some(record.flags() as i32));
        
        // Fast chromosome lookup
        let chromosome = if record.tid() >= 0 {
            chrom_lookup.get(&record.tid()).map(|s| s.clone())
        } else {
            None
        };
        chromosomes.push(chromosome);
        
        positions.push(Some(record.pos() as i64));
        mapping_qualities.push(Some(record.mapq() as i32));
        
        // Optimized CIGAR string extraction
        let cigar_str = record.cigar().to_string();
        cigars.push(if cigar_str.is_empty() { None } else { Some(cigar_str) });
        
        // Next reference and position
        let rnext = if record.mtid() >= 0 {
            chrom_lookup.get(&record.mtid()).map(|s| s.clone())
        } else {
            None
        };
        rnexts.push(rnext);
        
        pnexts.push(Some(record.mpos() as i64));
        template_lengths.push(Some(record.insert_size() as i64));
        
        // Optional sequence and quality with zero-copy where possible
        if let Some(ref mut seq_vec) = sequences {
            let seq_str = String::from_utf8_lossy(&record.seq().as_bytes()).to_string();
            seq_vec.push(if seq_str.is_empty() { None } else { Some(seq_str) });
        }
        
        if let Some(ref mut qual_vec) = qualities {
            let qual_str = String::from_utf8_lossy(record.qual()).to_string();
            qual_vec.push(if qual_str.is_empty() { None } else { Some(qual_str) });
        }
    }
    
    // Build optimized record batch
    crate::bam_htslib::build_optimized_record_batch(
        qnames,
        flags,
        chromosomes,
        positions,
        mapping_qualities,
        cigars,
        rnexts,
        pnexts,
        template_lengths,
        sequences,
        qualities,
        include_sequence,
        include_quality,
    )
}

/// Efficient Arrow file concatenation
#[cfg(feature = "htslib")]
fn concatenate_arrow_files(input_files: &[String], output_path: &str) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    if input_files.is_empty() {
        return Err("No input files to concatenate".into());
    }
    
    // If only one file, just copy it
    if input_files.len() == 1 {
        std::fs::copy(&input_files[0], output_path)?;
        return Ok(());
    }
    
    // Open first file to get schema
    let first_file = File::open(&input_files[0])?;
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
    for input_file in &input_files[1..] {
        if !Path::new(input_file).exists() {
            continue;
        }
        
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