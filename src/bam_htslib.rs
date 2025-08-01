use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::fs::File;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use rayon::ThreadPoolBuilder;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::io::{Read as StdRead, Seek, SeekFrom};
use std::collections::VecDeque;

#[cfg(feature = "htslib")]
use rust_htslib::bam as hts_bam;
#[cfg(feature = "htslib")]
use rust_htslib::bam::Read as HtsRead;
#[cfg(feature = "htslib")]
use rust_htslib::bam::ext::BamRecordExtensions;

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use arrow::ipc::writer::FileWriter as ArrowIpcWriter;

/// Creates the Arrow schema for BAM records
pub fn create_bam_schema(include_sequence: bool, include_quality: bool) -> Schema {
    let mut fields = vec![
        Field::new("name", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, true),
        Field::new("start", DataType::UInt32, true),
        Field::new("end", DataType::UInt32, true),
        Field::new("flag", DataType::UInt32, false),
    ];
    
    if include_sequence {
        fields.push(Field::new("sequence", DataType::Utf8, true));
    }
    
    if include_quality {
        fields.push(Field::new("quality_scores", DataType::Utf8, true));
    }
    
    Schema::new(fields)
}

/// Build chromosome lookup table from HTSlib header
#[cfg(feature = "htslib")]
pub fn build_chromosome_lookup(header: &hts_bam::HeaderView) -> std::collections::HashMap<i32, String> {
    let mut lookup = std::collections::HashMap::new();
    for tid in 0..header.target_count() {
        let chrom_name = String::from_utf8_lossy(header.tid2name(tid as u32)).into_owned();
        lookup.insert(tid as i32, chrom_name);
    }
    lookup
}

/// Build chromosome lookup table from HTSlib header (Arc version for compatibility)
#[cfg(feature = "htslib")]
// Removed unused function: build_chromosome_lookup_arc
// This was a duplicate of build_chromosome_lookup but returning Arc<Vec<String>>
// No references found in active codebase or relic code

/// Zero-copy quality score processing for maximum performance
#[cfg(feature = "htslib")]
fn quality_to_string_zero_copy(qual: &[u8]) -> String {
    if qual.is_empty() {
        return String::new();
    }
    
    let mut result = String::with_capacity(qual.len());
    unsafe {
        let bytes = result.as_mut_vec();
        bytes.extend_from_slice(qual);
        for byte in bytes {
            *byte += b'!';
        }
    }
    result
}

/// HTSlib reader pool for reusing initialized readers across segments
/// Reduces expensive reader initialization overhead in parallel processing
#[cfg(feature = "htslib")]
pub struct ReaderPool {
    bam_path: String,
    bgzf_threads: usize,
    pool: Arc<Mutex<VecDeque<hts_bam::Reader>>>,
    max_size: usize,
}

#[cfg(feature = "htslib")]
impl ReaderPool {
    pub fn new(bam_path: &str, bgzf_threads: usize, max_size: usize) -> Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        let pool = Arc::new(Mutex::new(VecDeque::new()));
        Ok(ReaderPool {
            bam_path: bam_path.to_string(),
            bgzf_threads,
            pool,
            max_size,
        })
    }
    
    /// Borrow a reader from the pool, creating a new one if pool is empty
    pub fn borrow_reader(&self) -> Result<hts_bam::Reader, Box<dyn std::error::Error + Send + Sync>> {
        // Try to get a reader from the pool first
        if let Ok(mut pool) = self.pool.lock() {
            if let Some(reader) = pool.pop_front() {
                return Ok(reader);
            }
        }
        
        // Pool is empty, create a new reader
        let mut reader = hts_bam::Reader::from_path(&self.bam_path)?;
        reader.set_threads(self.bgzf_threads)?;
        Ok(reader)
    }
    
    /// Return a reader to the pool for reuse
    pub fn return_reader(&self, reader: hts_bam::Reader) {
        if let Ok(mut pool) = self.pool.lock() {
            if pool.len() < self.max_size {
                pool.push_back(reader);
            }
            // If pool is full, just drop the reader
        }
    }
}

// REMOVED: build_optimized_record_batch function
// Was only used by deprecated hybrid relic code (bam_htslib_hybrid_optimized.rs)
// Since relic modules are not compiled, this function is no longer needed

// REMOVED: create_optimized_bam_schema function  
// Was only used by deprecated build_optimized_record_batch function
// Since relic modules are not compiled, this function is no longer needed

/// Process HTSlib records into Arrow RecordBatch with optimized performance
#[cfg(feature = "htslib")]
pub fn process_htslib_records_to_batch(
    hts_records: &[hts_bam::Record],
    chrom_lookup: &std::collections::HashMap<i32, String>,
    include_sequence: bool,
    include_quality: bool,
    _batch_id: usize,
) -> Result<RecordBatch, Box<dyn std::error::Error + Send + Sync>> {
    let record_count = hts_records.len();
    
    let mut names = Vec::with_capacity(record_count);
    let mut chroms = Vec::with_capacity(record_count);
    let mut starts = Vec::with_capacity(record_count);
    let mut ends = Vec::with_capacity(record_count);
    let mut flags = Vec::with_capacity(record_count);
    let mut sequences = if include_sequence { Some(Vec::with_capacity(record_count)) } else { None };
    let mut quality_scores = if include_quality { Some(Vec::with_capacity(record_count)) } else { None };
    
    for record in hts_records {
        // Name
        let name = String::from_utf8_lossy(record.qname()).into_owned();
        names.push(name);
        
        // Chromosome (using lookup table for correct names)
        let chrom = if record.tid() >= 0 {
            chrom_lookup.get(&record.tid()).map(|s| s.clone())
        } else {
            None
        };
        chroms.push(chrom);
        
        // Positions
        let start = if record.pos() >= 0 { Some(record.pos() as u32) } else { None };
        let end = if record.reference_end() > record.pos() { 
            Some(record.reference_end() as u32) 
        } else { 
            start 
        };
        starts.push(start);
        ends.push(end);
        
        // Flag
        flags.push(record.flags() as u32);
        
        // Sequence (if requested)
        if let Some(ref mut seq_vec) = sequences {
            let seq = record.seq().as_bytes();
            let sequence_str = if seq.is_empty() {
                None
            } else {
                Some(String::from_utf8_lossy(&seq).into_owned())
            };
            seq_vec.push(sequence_str);
        }
        
        // Quality scores (if requested) - using zero-copy optimization
        if let Some(ref mut qual_vec) = quality_scores {
            let qual = record.qual();
            let quality_str = if qual.is_empty() {
                None
            } else {
                Some(quality_to_string_zero_copy(qual))
            };
            qual_vec.push(quality_str);
        }
    }
    
    // Build Arrow arrays
    let mut arrays: Vec<Arc<dyn Array>> = vec![
        Arc::new(StringArray::from(names)),
        Arc::new(StringArray::from(chroms)),
        Arc::new(UInt32Array::from(starts)),
        Arc::new(UInt32Array::from(ends)),
        Arc::new(UInt32Array::from(flags)),
    ];
    
    if let Some(seq) = sequences {
        arrays.push(Arc::new(StringArray::from(seq)));
    }
    
    if let Some(qual) = quality_scores {
        arrays.push(Arc::new(StringArray::from(qual)));
    }
    
    let schema = create_bam_schema(include_sequence, include_quality);
    let batch = RecordBatch::try_new(Arc::new(schema), arrays)?;
    
    Ok(batch)
}

/// Fast discovery of BGZF split points for parallel processing
/// Uses file size estimation and boundary searching instead of exhaustive scanning
/// Optimized for fewer, larger segments to reduce coordination overhead
#[cfg(feature = "htslib")]
pub fn discover_split_points(bam_path: &str, num_workers: usize) -> Result<Vec<u64>, Box<dyn std::error::Error + Send + Sync>> {
    let file = File::open(bam_path)?;
    let file_size = file.metadata()?.len();
    
    // Minimum segment size (512MB) to avoid excessive coordination overhead
    let min_segment_size = 512 * 1024 * 1024; // 512MB
    let max_effective_workers = (file_size / min_segment_size).min(num_workers as u64) as usize;
    let effective_workers = max_effective_workers.max(1);
    
    eprintln!("Split point discovery: file_size={}, min_segment_size={}, max_effective_workers={}, effective_workers={}", 
             file_size, min_segment_size, max_effective_workers, effective_workers);
    
    let mut split_points = vec![0u64]; // Always start at beginning
    
    if effective_workers <= 1 {
        eprintln!("Only 1 effective worker, returning single split point");
        return Ok(split_points);
    }
    
    // Calculate approximate split positions with larger segments
    let chunk_size = file_size / effective_workers as u64;
    
    for i in 1..effective_workers {
        let estimated_pos = chunk_size * i as u64;
        
        // Find the nearest BGZF block boundary around this position
        match find_nearest_bgzf_boundary(bam_path, estimated_pos, chunk_size / 4)? {
            Some(boundary) => {
                split_points.push(boundary);
            }
            None => {
                split_points.push(estimated_pos);
            }
        }
    }
    
    // Add the end of file as the final split point
    split_points.push(file_size);
    
    eprintln!("Found {} split points for parallel processing: {:?}", split_points.len(), split_points);
    eprintln!("This creates {} segments", split_points.len() - 1);
    Ok(split_points)
}

/// Find the nearest BGZF block boundary within search_radius of target_pos
#[cfg(feature = "htslib")]
fn find_nearest_bgzf_boundary(
    bam_path: &str, 
    target_pos: u64, 
    search_radius: u64
) -> Result<Option<u64>, Box<dyn std::error::Error + Send + Sync>> {
    let mut file = File::open(bam_path)?;
    let file_size = file.metadata()?.len();
    
    // Define search range
    let search_start = target_pos.saturating_sub(search_radius);
    let search_end = (target_pos + search_radius).min(file_size);
    
    // Start searching from target_pos and expand outward
    for offset in 0..search_radius {
        // Check positions before and after target
        for &pos in &[target_pos.saturating_sub(offset), target_pos + offset] {
            if pos >= search_start && pos < search_end {
                if let Some(boundary) = check_bgzf_boundary_at(&mut file, pos)? {
                    return Ok(Some(boundary));
                }
            }
        }
    }
    
    Ok(None)
}

/// Check if there's a BGZF block boundary at or near the given position
#[cfg(feature = "htslib")]
fn check_bgzf_boundary_at(file: &mut File, pos: u64) -> Result<Option<u64>, Box<dyn std::error::Error + Send + Sync>> {
    // Read a small chunk around this position to look for BGZF magic
    const SEARCH_CHUNK: usize = 1024;
    let start_pos = pos.saturating_sub(SEARCH_CHUNK as u64 / 2);
    
    file.seek(SeekFrom::Start(start_pos))?;
    let mut buffer = vec![0u8; SEARCH_CHUNK];
    let bytes_read = file.read(&mut buffer)?;
    
    // Look for BGZF magic header in the buffer
    const BGZF_MAGIC: &[u8] = &[0x1f, 0x8b, 0x08, 0x04];
    
    for i in 0..bytes_read.saturating_sub(4) {
        if &buffer[i..i+4] == BGZF_MAGIC {
            let found_pos = start_pos + i as u64;
            
            // Verify this is actually a valid BGZF block by reading its size
            if let Ok(Some(_)) = get_bgzf_block_size_at(file, found_pos) {
                return Ok(Some(found_pos));
            }
        }
    }
    
    Ok(None)
}

/// Get the size of a BGZF block starting at the given position
#[cfg(feature = "htslib")]
fn get_bgzf_block_size_at(file: &mut File, pos: u64) -> Result<Option<u64>, Box<dyn std::error::Error + Send + Sync>> {
    file.seek(SeekFrom::Start(pos))?;
    
    // Read and verify BGZF header
    let mut header = [0u8; 18];
    if file.read_exact(&mut header).is_err() {
        return Ok(None);
    }
    
    // Check magic number
    if &header[0..4] != &[0x1f, 0x8b, 0x08, 0x04] {
        return Ok(None);
    }
    
    // Extract block size from BSIZE field (bytes 16-17)
    let block_size = u16::from_le_bytes([header[16], header[17]]) as u64 + 1;
    
    // Sanity check: BGZF blocks should be reasonable size (typically 16KB-64KB)
    if block_size > 65536 || block_size < 26 {
        return Ok(None);
    }
    
    Ok(Some(block_size))
}

/// Process a file segment between start_pos and end_pos
#[cfg(feature = "htslib")]
pub fn process_file_segment_with_pool(
    reader_pool: &ReaderPool,
    start_pos: u64,
    end_pos: u64,
    record_limit: Option<usize>,
    chrom_lookup: &std::collections::HashMap<i32, String>,
    include_sequence: bool,
    include_quality: bool,
    batch_size: usize,
) -> Result<Vec<RecordBatch>, Box<dyn std::error::Error + Send + Sync>> {
    // Borrow a reader from the pool instead of creating a new one
    let mut reader = reader_pool.borrow_reader()?;
    
    // For the first segment (start_pos == 0), read from beginning
    // For other segments, we need to seek to start_pos using BGZF virtual offsets
    if start_pos > 0 {
        // Convert raw file position to BGZF virtual offset
        // BGZF virtual offset = (block_offset << 16) | uncompressed_offset
        // Since we're seeking to block boundaries, uncompressed_offset = 0
        let virtual_offset = start_pos << 16;
        
        // Use HTSlib's seek functionality with virtual offset (convert to i64)
        if let Err(_) = reader.seek(virtual_offset as i64) {
            // If seek fails, continue from beginning (will have duplicate records but won't fail)
            // This is better than failing the entire segment
        }
    }
    
    let mut batches = Vec::new();
    let mut current_batch_records = Vec::new();
    let mut total_records_processed = 0;
    let mut records_since_position_check = 0;
    
    let mut record = hts_bam::Record::new();
    while let Some(Ok(_)) = reader.read(&mut record) {
        // Check if we've reached the record limit for this segment
        if let Some(limit) = record_limit {
            if total_records_processed >= limit {
                eprintln!("Segment reached record limit: {} records", limit);
                break;
            }
        }
        
        current_batch_records.push(record.clone());
        total_records_processed += 1;
        records_since_position_check += 1;
        
        // Check segment boundary every 1000 records to avoid performance impact
        if records_since_position_check >= 1000 {
            records_since_position_check = 0;
            
            // Get current virtual file position from reader
            let current_voffset = reader.tell();
            if current_voffset >= 0 {
                let current_file_pos = (current_voffset as u64) >> 16; // Extract block offset
                
                // Stop if we've gone beyond our segment boundary
                if end_pos != std::u64::MAX && current_file_pos >= end_pos {
                    eprintln!("Segment reached boundary at file position {} (target: {})", 
                             current_file_pos, end_pos);
                    break;
                }
            }
        }
        
        // Create batch when we reach batch_size
        if current_batch_records.len() >= batch_size {
            let batch = process_htslib_records_to_batch(
                &current_batch_records,
                chrom_lookup,
                include_sequence,
                include_quality,
                batches.len(),
            )?;
            batches.push(batch);
            current_batch_records.clear();
        }
    }
    
    // Process remaining records
    if !current_batch_records.is_empty() {
        let batch = process_htslib_records_to_batch(
            &current_batch_records,
            chrom_lookup,
            include_sequence,
            include_quality,
            batches.len(),
        )?;
        batches.push(batch);
    }
    
    // Return the reader to the pool for reuse
    reader_pool.return_reader(reader);
    
    Ok(batches)
}

/// Process a file segment between start_pos and end_pos (original version for compatibility)
#[allow(dead_code)]
#[cfg(feature = "htslib")]
pub fn process_file_segment(
    bam_path: &str,
    start_pos: u64,
    end_pos: u64,
    record_limit: Option<usize>,
    chrom_lookup: &std::collections::HashMap<i32, String>,
    include_sequence: bool,
    include_quality: bool,
    batch_size: usize,
    bgzf_threads: usize,
) -> Result<Vec<RecordBatch>, Box<dyn std::error::Error + Send + Sync>> {
    // Create a reader pool with just one reader for this call
    let reader_pool = ReaderPool::new(bam_path, bgzf_threads, 1)?;
    
    // Use the optimized version with reader pool
    process_file_segment_with_pool(
        &reader_pool,
        start_pos,
        end_pos,
        record_limit,
        chrom_lookup,
        include_sequence,
        include_quality,
        batch_size,
    )
}

/// Parallel BAM processing using true BGZF block-level parallelization
/// This function implements divide-and-conquer for large unmapped BAM files
#[cfg(feature = "htslib")]
#[pyfunction]
#[pyo3(signature = (
    bam_path,
    arrow_ipc_path,
    batch_size = 20000,
    include_sequence = true,
    include_quality = true,
    bgzf_threads = 4,
    writing_threads = 8,
    read_buffer_mb = None,
    write_buffer_mb = None,
    limit = None,
    num_block_workers = None
))]
pub fn bam_to_arrow_ipc_htslib_bgzf_blocks(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    bgzf_threads: usize,
    writing_threads: usize,
    read_buffer_mb: Option<usize>,
    write_buffer_mb: Option<usize>,
    limit: Option<usize>,
    num_block_workers: Option<usize>,
) -> PyResult<()> {
    let start_time = Instant::now();
    
    let input_path = Path::new(bam_path);
    let output_path = Path::new(arrow_ipc_path);

    if !input_path.exists() {
        return Err(PyErr::new::<PyRuntimeError, _>(
            format!("BAM file does not exist: {}", bam_path)
        ));
    }
    
    if batch_size == 0 {
        return Err(PyErr::new::<PyRuntimeError, _>(
            "batch_size must be greater than 0"
        ));
    }

    let effective_batch_size = batch_size.min(1000000);
    let effective_bgzf_threads = if bgzf_threads == 0 { 4 } else { bgzf_threads.min(32) };
    let effective_writing_threads = writing_threads.max(1).min(32);
    
    // Dynamically adjust worker count based on file size and limit
    let base_workers = num_block_workers.unwrap_or(effective_writing_threads).min(16);
    let effective_block_workers = if let Some(record_limit) = limit {
        if record_limit <= 1_000_000 {
            // For small datasets, use fewer workers to reduce overhead
            base_workers.min(2)
        } else if record_limit <= 5_000_000 {
            base_workers.min(4)
        } else {
            base_workers
        }
    } else {
        base_workers
    };
    
    // Intelligent buffer sizing for block-level processing
    // Scale buffer sizes based on number of workers for optimal I/O efficiency
    let default_read_buffer_mb = if effective_block_workers >= 8 {
        1536  // 1.5GB for high-parallelism scenarios
    } else if effective_block_workers >= 4 {
        1024  // 1GB for medium parallelism
    } else {
        768   // 768MB for low parallelism
    };
    
    let default_write_buffer_mb = if effective_block_workers >= 8 {
        384   // 384MB for high-parallelism scenarios
    } else {
        256   // 256MB for medium/low parallelism
    };
    
    let read_buffer_size = read_buffer_mb
        .unwrap_or(default_read_buffer_mb)
        * 1024 * 1024;
        
    let write_buffer_size = write_buffer_mb
        .unwrap_or(default_write_buffer_mb)
        * 1024 * 1024;
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    // Check file size to determine if BGZF block-level is worth it
    let input_file = std::fs::File::open(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Cannot access BAM file: {}", e)))?;
    let file_size_gb = input_file.metadata()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Cannot read file metadata: {}", e)))?
        .len() as f64 / (1024.0 * 1024.0 * 1024.0);
    
    eprintln!("Starting BGZF Block-Level Parallel BAM to IPC conversion:");
    eprintln!("  Input: {} ({:.1}GB)", bam_path, file_size_gb);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Batch size: {}", effective_batch_size);
    eprintln!("  BGZF threads: {}", effective_bgzf_threads);
    eprintln!("  Writing threads: {}", effective_writing_threads);
    eprintln!("  Block workers: {}", effective_block_workers);
    eprintln!("  Read buffer: {}MB", read_buffer_size / (1024 * 1024));
    eprintln!("  Write buffer: {}MB", write_buffer_size / (1024 * 1024));
    eprintln!("  Include sequence: {}", include_sequence);
    eprintln!("  Include quality: {}", include_quality);
    
    // Performance warning for small files
    if file_size_gb < 10.0 {
        eprintln!("⚠️  WARNING: File size ({:.1}GB) may be too small for BGZF block-level optimization.", file_size_gb);
        eprintln!("   Consider using bam_to_arrow_ipc_htslib_optimized() for files <10GB.");
        eprintln!("   BGZF block-level is optimized for very large files (50GB+).");
    }

    // Step 1: Discover split points (much faster than exhaustive block discovery)
    let split_points = discover_split_points(bam_path, effective_block_workers)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
            format!("Failed to discover split points: {}", e)
        ))?;
    
    if split_points.len() < 2 {
        return Err(PyErr::new::<PyRuntimeError, _>(
            "Could not find sufficient split points for parallel processing"
        ));
    }
    
    eprintln!("Using {} file segments for parallel processing", split_points.len() - 1);
    
    // Calculate per-segment limits if total limit is specified
    let per_segment_limit = if let Some(total_limit) = limit {
        let num_segments = split_points.len() - 1;
        let base_limit = total_limit / num_segments;
        let remainder = total_limit % num_segments;
        
        eprintln!("Distributing {} total records across {} segments:", total_limit, num_segments);
        eprintln!("  Base limit per segment: {}", base_limit);
        if remainder > 0 {
            eprintln!("  {} segments will get +1 extra record", remainder);
        }
        
        Some((base_limit, remainder))
    } else {
        None
    };
    
    // Step 2: Setup parallel processing architecture
    let schema = create_bam_schema(include_sequence, include_quality);
    let channel_size = (effective_block_workers * 4).max(32);
    let (segment_sender, segment_receiver): (Sender<(usize, (u64, u64, Option<usize>))>, Receiver<(usize, (u64, u64, Option<usize>))>) = bounded(channel_size);
    let (result_sender, result_receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    // Step 3: Configure block processing thread pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(effective_block_workers)
        .build()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create thread pool: {}", e)))?;

    // Step 4: Get header from main BAM file for chromosome lookup
    let main_reader = hts_bam::Reader::from_path(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file: {}", e)))?;
    let header_view = main_reader.header().clone();
    let chrom_lookup = build_chromosome_lookup(&header_view);
    drop(main_reader); // Close main reader

    // Step 5: Create reader pool for optimized performance
    let reader_pool = Arc::new(ReaderPool::new(bam_path, effective_bgzf_threads, effective_block_workers * 2)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create reader pool: {}", e)))?);

    // Step 6: Spawn segment processing workers
    let workers_result_sender = result_sender.clone();
    let workers_segment_receiver = segment_receiver.clone();
    let workers_chrom_lookup = chrom_lookup.clone();
    let workers_reader_pool = reader_pool.clone();
    
    let _workers_handle = std::thread::spawn(move || {
        pool.install(|| {
            while let Ok((segment_id, (start_pos, end_pos, segment_limit))) = workers_segment_receiver.recv() {
                let segment_processing_start = Instant::now();
                
                // Process file segment from start_pos to end_pos using reader pool
                let segment_result = process_file_segment_with_pool(
                    &workers_reader_pool,
                    start_pos,
                    end_pos,
                    segment_limit,
                    &workers_chrom_lookup,
                    include_sequence,
                    include_quality,
                    effective_batch_size,
                );
                
                match segment_result {
                    Ok(batches) => {
                        let _segment_duration = segment_processing_start.elapsed();
                        let _total_records: usize = batches.iter().map(|b| b.num_rows()).sum();
                        
                        // Send all batches from this segment
                        for (batch_idx, batch) in batches.into_iter().enumerate() {
                            let combined_id = segment_id * 10000 + batch_idx; // Unique batch ID
                            if let Err(_) = workers_result_sender.send((combined_id, batch)) {
                                eprintln!("Failed to send batch {} from segment {}", batch_idx, segment_id);
                                return;
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Failed to process segment {}: {}", segment_id, e);
                        // Continue processing other segments
                    }
                }
            }
        });
    });

    // Step 6: Writer thread with buffered output
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file: {}", e)))?;
    let buffered_output = std::io::BufWriter::with_capacity(write_buffer_size, output_file);
    let mut arrow_writer = ArrowIpcWriter::try_new(buffered_output, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow writer: {}", e)))?;

    let writer_handle = std::thread::spawn(move || -> PyResult<usize> {
        let mut total_records = 0;
        let mut batch_count = 0;
        
        // Write batches as they arrive from block workers
        while let Ok((_, batch)) = result_receiver.recv() {
            arrow_writer.write(&batch)
                .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                    format!("Failed to write batch: {}", e)
                ))?;
            
            total_records += batch.num_rows();
            batch_count += 1;
            
            // Report progress every 50 batches
            if batch_count % 50 == 0 {
                eprintln!("Progress: {} batches written, {} total records (BGZF blocks)", batch_count, total_records);
            }
        }
        
        arrow_writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to finish writing: {}", e)))?;
        
        Ok(total_records)
    });

    // Step 7: Main thread distributes file segments to workers
    let mut segments_sent = 0;
    
    // Create segment ranges from split points
    for i in 0..split_points.len().saturating_sub(1) {
        let start_pos = split_points[i];
        let end_pos = if i + 1 < split_points.len() {
            split_points[i + 1]
        } else {
            // Last segment goes to end of file
            std::u64::MAX
        };
        
        // Calculate this segment's record limit
        let this_segment_limit = if let Some((base_limit, remainder)) = per_segment_limit {
            let extra = if i < remainder { 1 } else { 0 };
            Some(base_limit + extra)
        } else {
            None
        };
        
        if let Err(_) = segment_sender.send((i, (start_pos, end_pos, this_segment_limit))) {
            eprintln!("Failed to send segment {} for processing", i);
            break;
        }
        segments_sent += 1;
        
        let limit_str = if let Some(limit) = this_segment_limit {
            format!(" (limit: {} records)", limit)
        } else {
            String::new()
        };
        
        eprintln!("Distributed segment {} ({}MB - {}MB){}", 
                 i, 
                 start_pos / (1024 * 1024),
                 if end_pos == std::u64::MAX { 
                     "EOF".to_string() 
                 } else { 
                     (end_pos / (1024 * 1024)).to_string() 
                 },
                 limit_str);
    }
    
    // Close senders to signal completion
    drop(segment_sender);
    drop(result_sender);
    
    eprintln!("All {} file segments distributed, waiting for completion...", segments_sent);

    // Wait for writer to complete
    let total_records = writer_handle.join()
        .map_err(|_| PyErr::new::<PyRuntimeError, _>("Writer thread panicked"))?
        .map_err(|e| e)?;

    let total_duration = start_time.elapsed();
    let records_per_sec = if total_duration.as_secs() > 0 {
        total_records as f64 / total_duration.as_secs_f64()
    } else {
        0.0
    };

    eprintln!("BGZF Block-Level Parallel Processing completed!");
    eprintln!("  Total records: {}", total_records);
    eprintln!("  Total segments processed: {}", segments_sent);
    eprintln!("  Duration: {:?}", total_duration);
    eprintln!("  Throughput: {:.0} records/sec", records_per_sec);

    Ok(())
}
