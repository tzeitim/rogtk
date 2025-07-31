use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use pyo3::Python;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;
use rayon::ThreadPoolBuilder;
use crossbeam_channel::{bounded, Receiver, Sender};
use log::debug;
// GZP imports temporarily disabled - using enhanced I/O buffering approach instead
// use gzp::par::decompress::ParDecompressBuilder;
// use gzp::deflate::Mgzip;
// TEMPORARILY DISABLED: HTSlib imports - requires libclang system dependency
#[cfg(feature = "htslib")]
use rust_htslib::bam as hts_bam;
#[cfg(feature = "htslib")]
use rust_htslib::bam::Read as HtsRead;

use noodles::{bam, sam, bgzf};
use sam::alignment::record::cigar::op::Kind as CigarOpKind;
use sam::alignment::record::QualityScores; 
use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use arrow::ipc::writer::FileWriter as ArrowIpcWriter;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, Encoding};
use parquet::file::properties::WriterProperties;



// Struct to hold reusable buffers
struct ReusableBuffers {
    names: Vec<String>,
    chroms: Vec<Option<String>>,
    starts: Vec<Option<u32>>,
    ends: Vec<Option<u32>>,
    flags: Vec<u32>,
    sequences: Option<Vec<Option<String>>>,
    quality_scores: Option<Vec<Option<String>>>,
    // Intermediate string buffers to avoid allocations
    temp_sequence: String,
    temp_quality: String,
}

impl ReusableBuffers {
    fn new(capacity: usize) -> Self {
        Self {
            names: Vec::with_capacity(capacity),
            chroms: Vec::with_capacity(capacity),
            starts: Vec::with_capacity(capacity),
            ends: Vec::with_capacity(capacity),
            flags: Vec::with_capacity(capacity),
            sequences: Some(Vec::with_capacity(capacity)),
            quality_scores: Some(Vec::with_capacity(capacity)),
            temp_sequence: String::with_capacity(1024),
            temp_quality: String::with_capacity(1024),
        }
    }
    
    fn clear(&mut self) {
        self.names.clear();
        self.chroms.clear();
        self.starts.clear();
        self.ends.clear();
        self.flags.clear();
        if let Some(ref mut seq) = self.sequences {
            seq.clear();
        }
        if let Some(ref mut qual) = self.quality_scores {
            qual.clear();
        }
        // Clear but retain capacity for temp strings
        self.temp_sequence.clear();
        self.temp_quality.clear();
    }
    
    fn shrink_to_fit_if_needed(&mut self) {
        // Periodically shrink buffers if they've grown too large
        if self.names.capacity() > 500000 {  // Higher threshold for shrinking
            self.names.shrink_to_fit();
            self.chroms.shrink_to_fit();
            self.starts.shrink_to_fit();
            self.ends.shrink_to_fit();
            self.flags.shrink_to_fit();
            if let Some(ref mut seq) = self.sequences {
                seq.shrink_to_fit();
            }
            if let Some(ref mut qual) = self.quality_scores {
                qual.shrink_to_fit();
            }
        }
    }
}

fn read_bam_batch_enhanced(
    reader: &mut bam::io::Reader<bgzf::Reader<&mut File>>,
    header: &sam::Header,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    buffers: &mut ReusableBuffers,
) -> PyResult<RecordBatch> {
    // Clear previous data but retain allocated capacity
    buffers.clear();
    
    let mut count = 0;
    let mut record = bam::Record::default();

    while count < batch_size {
        match reader.read_record(&mut record) {
            Ok(0) => break, // EOF
            Ok(_) => {
                extract_record_data_enhanced(
                    &record, 
                    header, 
                    buffers,
                    include_sequence,
                    include_quality,
                )?;
                count += 1;
            }
            Err(e) => {
                return Err(PyErr::new::<PyRuntimeError, _>(
                    format!("Failed to read BAM record at position {}: {}", count, e)
                ));
            }
        }
    }
    
    // Shrink buffers periodically to prevent runaway memory growth
    if count % 10000 == 0 {  // Less frequent shrinking with larger batches
        buffers.shrink_to_fit_if_needed();
    }

    // Create Arrow arrays from buffers
    let name_array = Arc::new(StringArray::from(buffers.names.clone()));
    let chrom_array = Arc::new(StringArray::from(buffers.chroms.clone()));
    let start_array = Arc::new(UInt32Array::from(buffers.starts.clone()));
    let end_array = Arc::new(UInt32Array::from(buffers.ends.clone()));
    let flags_array = Arc::new(UInt32Array::from(buffers.flags.clone()));
    
    let mut arrays: Vec<ArrayRef> = vec![
        name_array,
        chrom_array, 
        start_array,
        end_array,
        flags_array,
    ];
    
    if include_sequence {
        if let Some(ref sequences) = buffers.sequences {
            arrays.push(Arc::new(StringArray::from(sequences.clone())));
        }
    }
    
    if include_quality {
        if let Some(ref quality_scores) = buffers.quality_scores {
            arrays.push(Arc::new(StringArray::from(quality_scores.clone())));
        }
    }

    let schema = create_bam_schema(include_sequence, include_quality);
    RecordBatch::try_new(schema, arrays)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create record batch: {}", e)))
}

fn extract_record_data_enhanced(
    record: &bam::Record,
    header: &sam::Header,
    buffers: &mut ReusableBuffers,
    include_sequence: bool,
    include_quality: bool,
) -> PyResult<()> {
    // 1. Extract read name - exact same pattern as existing extract_record_data
    let name = record.name()
        .map(|bstr| bstr.to_string())
        .unwrap_or_else(|| "unknown".to_string());
    buffers.names.push(name);

    // 2. Extract chromosome/reference name - exact same logic as existing code
    let chrom = if let Some(reference_sequence_id) = record.reference_sequence_id() {
        // reference_sequence_id is a Result, so we need to handle it
        match reference_sequence_id {
            Ok(id) => {
                header.reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| name.to_string())
            }
            Err(_) => None // Invalid reference ID
        }
    } else {
        None // Unmapped read
    };
    buffers.chroms.push(chrom);

    // 3. Extract coordinates (convert to 1-based to match polars-bio) - exact same logic
    let (start_pos, end_pos) = if let Some(alignment_start) = record.alignment_start() {
        match alignment_start {
            Ok(pos) => {
                let start_1based = pos.get() as u32; // noodles Position.get() returns usize
                let end_1based = start_1based + calculate_bam_alignment_length(&record.cigar()) - 1;
                (Some(start_1based), Some(end_1based))
            }
            Err(_) => (None, None) // Invalid position
        }
    } else {
        (None, None) // Unmapped read
    };
    buffers.starts.push(start_pos);
    buffers.ends.push(end_pos);

    // 4. Extract flags - exact same as existing code
    buffers.flags.push(record.flags().bits() as u32); // Convert u16 to u32

    // 5. Extract sequence if requested - exact same complex logic as existing code
    if include_sequence {
        if let Some(ref mut seq_vec) = buffers.sequences {
            // Get the sequence object and keep it alive
            let sequence_obj = record.sequence();
            let seq_len = sequence_obj.len();
            let sequence_bytes = sequence_obj.as_ref();
            buffers.temp_sequence.clear();
            buffers.temp_sequence.reserve(seq_len);
            
            // Each byte contains 2 bases (4 bits each)
            for i in 0..seq_len {
                let byte_idx = i / 2;
                let base_encoded = if i % 2 == 0 {
                    // Even index: higher 4 bits
                    (sequence_bytes[byte_idx] >> 4) & 0x0F
                } else {
                    // Odd index: lower 4 bits
                    sequence_bytes[byte_idx] & 0x0F
                };
                buffers.temp_sequence.push(decode_base(base_encoded));
            }
            
            seq_vec.push(if buffers.temp_sequence.is_empty() { None } else { Some(buffers.temp_sequence.clone()) });
        }
    }

    // 6. Extract quality scores if requested - FIXED: use proper lifetime management
    if include_quality {
        if let Some(ref mut qual_vec) = buffers.quality_scores {
            buffers.temp_quality.clear();
            // Fix: Store quality_scores directly and iterate over it
            let quality_scores = record.quality_scores();
            for q in quality_scores.iter() {
                buffers.temp_quality.push(char::from(q + b'!')); // Convert to ASCII PHRED+33
            }
            qual_vec.push(if buffers.temp_quality.is_empty() { None } else { Some(buffers.temp_quality.clone()) });
        }
    }

    Ok(())
}
// ####################
// ####################
// ####################
// ####################
#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    parquet_path, 
    batch_size = 50000,  // Balanced default batch size for optimal memory/performance trade-off
    include_sequence = true,
    include_quality = true,
    compression = "snappy",
    limit = None
))]
pub fn bam_to_parquet(
    bam_path: &str,
    parquet_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    compression: &str,
    limit: Option<usize>,
) -> PyResult<()> {
    let input_path = Path::new(bam_path);
    let output_path = Path::new(parquet_path);

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

    // Cap batch size to prevent excessive memory usage only for extremely large values
    let effective_batch_size = batch_size.min(1000000);  // Much higher cap
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    let mut file = File::open(input_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file '{}': {}", bam_path, e)))?;
    
    let mut bam_reader = bam::io::Reader::new(&mut file);
    
    let header = bam_reader.read_header()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to read BAM header: {}", e)))?;

    let schema = create_bam_schema(include_sequence, include_quality);
    
    let writer_props = WriterProperties::builder()
        .set_compression(parse_compression(compression))
        .set_encoding(Encoding::PLAIN)
        .build();

    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file '{}': {}", parquet_path, e)))?;

    let mut parquet_writer = ArrowWriter::try_new(output_file, schema.clone(), Some(writer_props))
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Parquet writer: {}", e)))?;

    // Pre-allocate reusable buffers outside the loop
    let mut reusable_buffers = ReusableBuffers::new(effective_batch_size);
    
    let mut batch_count = 0;
    let mut total_records = 0;
    let target_records = limit.unwrap_or(usize::MAX);

    loop {
        // Check for Python interrupts every 10 batches
        if batch_count % 10 == 0 {
            Python::with_gil(|py| {
                py.check_signals().map_err(|e| {
                    eprintln!("Conversion interrupted by user");
                    e
                })
            })?;
        }

        let remaining_records = target_records.saturating_sub(total_records);
        if remaining_records == 0 {
            eprintln!("Reached limit of {} records", target_records);
            break;
        }

        // Adjust batch size if we're approaching the limit
        let current_batch_size = effective_batch_size.min(remaining_records);
        
        let batch = read_bam_batch_enhanced(
            &mut bam_reader, 
            &header, 
            current_batch_size, 
            include_sequence, 
            include_quality,
            &mut reusable_buffers  // Pass reusable buffers
        )?;
        
        if batch.num_rows() == 0 {
            break;
        }

        parquet_writer.write(&batch)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to write batch {}: {}", batch_count, e)))?;

        total_records += batch.num_rows();
        batch_count += 1;
        
        // Reduce frequency: report every 200 batches instead of 50
        if batch_count % 200 == 0 {
            let progress_msg = if let Some(limit_val) = limit {
                format!("Processed {} batches, {} / {} records ({:.1}%) - Batch size: {}", 
                       batch_count, total_records, limit_val, 
                       100.0 * total_records as f64 / limit_val as f64,
                       current_batch_size)
            } else {
                format!("Processed {} batches, {} total records - Batch size: {}", 
                       batch_count, total_records, current_batch_size)
            };
            eprintln!("Progress: {}", progress_msg);
            
            // Force garbage collection every 400 batches  
            if batch_count % 400 == 0 {
                Python::with_gil(|py| {
                    let gc = py.import_bound("gc").unwrap();
                    let _ = gc.call_method0("collect");
                });
            }
        }
    }

    parquet_writer.close()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Parquet writer: {}", e)))?;

    let completion_msg = if limit.is_some() {
        format!("Conversion complete: {} records (limited from potentially more) written to {}", 
               total_records, parquet_path)
    } else {
        format!("Conversion complete: {} records written to {}", total_records, parquet_path)
    };
    eprintln!("{}", completion_msg);
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    arrow_ipc_path, 
    batch_size = 50000,
    include_sequence = true,
    include_quality = true,
    limit = None
))]
pub fn bam_to_arrow_ipc(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    limit: Option<usize>,
) -> PyResult<()> {
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
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    let mut file = File::open(input_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file '{}': {}", bam_path, e)))?;
    
    let mut bam_reader = bam::io::Reader::new(&mut file);
    
    let header = bam_reader.read_header()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to read BAM header: {}", e)))?;

    let schema = create_bam_schema(include_sequence, include_quality);
    
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file '{}': {}", arrow_ipc_path, e)))?;

    let mut arrow_writer = ArrowIpcWriter::try_new(output_file, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow IPC writer: {}", e)))?;

    let mut reusable_buffers = ReusableBuffers::new(effective_batch_size);
    
    let mut batch_count = 0;
    let mut total_records = 0;
    let target_records = limit.unwrap_or(usize::MAX);

    loop {
        if batch_count % 10 == 0 {
            Python::with_gil(|py| {
                py.check_signals().map_err(|e| {
                    eprintln!("Conversion interrupted by user");
                    e
                })
            })?;
        }

        let remaining_records = target_records.saturating_sub(total_records);
        if remaining_records == 0 {
            eprintln!("Reached limit of {} records", target_records);
            break;
        }

        let current_batch_size = effective_batch_size.min(remaining_records);
        
        let batch = read_bam_batch_enhanced(
            &mut bam_reader, 
            &header, 
            current_batch_size, 
            include_sequence, 
            include_quality,
            &mut reusable_buffers
        )?;
        
        if batch.num_rows() == 0 {
            break;
        }

        arrow_writer.write(&batch)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to write batch {}: {}", batch_count, e)))?;

        total_records += batch.num_rows();
        batch_count += 1;
        
        // Reduce frequency: report every 200 batches instead of 50
        if batch_count % 200 == 0 {
            let progress_msg = if let Some(limit_val) = limit {
                format!("Processed {} batches, {} / {} records ({:.1}%) - Batch size: {}", 
                       batch_count, total_records, limit_val, 
                       100.0 * total_records as f64 / limit_val as f64,
                       current_batch_size)
            } else {
                format!("Processed {} batches, {} total records - Batch size: {}", 
                       batch_count, total_records, current_batch_size)
            };
            eprintln!("Progress: {}", progress_msg);
            
            // Force garbage collection every 400 batches
            if batch_count % 400 == 0 {
                Python::with_gil(|py| {
                    let gc = py.import_bound("gc").unwrap();
                    let _ = gc.call_method0("collect");
                });
            }
        }
    }

    arrow_writer.finish()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Arrow IPC writer: {}", e)))?;

    let completion_msg = if limit.is_some() {
        format!("Conversion complete: {} records (limited from potentially more) written to {}", 
               total_records, arrow_ipc_path)
    } else {
        format!("Conversion complete: {} records written to {}", total_records, arrow_ipc_path)
    };
    eprintln!("{}", completion_msg);
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    arrow_ipc_path, 
    batch_size = 50000,
    include_sequence = true,
    include_quality = true,
    num_threads = 4,
    preserve_order = false,
    limit = None
))]
pub fn bam_to_arrow_ipc_parallel(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    num_threads: usize,
    preserve_order: bool,
    limit: Option<usize>,
) -> PyResult<()> {
    let start_time = Instant::now();
    
    // Thread validation with nice user messages
    let effective_threads = if num_threads == 0 {
        eprintln!("Warning: num_threads cannot be 0, using default of 4 threads");
        4
    } else if num_threads > 8 {
        eprintln!("Notice: For optimal performance, we recommend using 8 threads or fewer based on our benchmarking.");
        eprintln!("Your request for {} threads has been capped to 8 threads for best results.", num_threads);
        8
    } else {
        num_threads
    };
    
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
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    eprintln!("Starting parallel BAM to IPC conversion:");
    eprintln!("  Input: {}", bam_path);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Batch size: {}", effective_batch_size);
    eprintln!("  Threads: {}", effective_threads);
    eprintln!("  Preserve order: {}", preserve_order);
    eprintln!("  Include sequence: {}", include_sequence);
    eprintln!("  Include quality: {}", include_quality);

    let mut file = File::open(input_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file '{}': {}", bam_path, e)))?;
    
    let mut bam_reader = bam::io::Reader::new(&mut file);
    
    let header = bam_reader.read_header()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to read BAM header: {}", e)))?;

    let schema = create_bam_schema(include_sequence, include_quality);
    
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file '{}': {}", arrow_ipc_path, e)))?;

    let mut arrow_writer = ArrowIpcWriter::try_new(output_file, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow IPC writer: {}", e)))?;

    // Setup parallel processing architecture
    let channel_size = (effective_threads * 4).max(16);
    let (batch_sender, batch_receiver): (Sender<(usize, Vec<bam::Record>)>, Receiver<(usize, Vec<bam::Record>)>) = bounded(channel_size);
    let (result_sender, result_receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    let header_arc = Arc::new(header.clone());
    
    // Configure thread pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(effective_threads)
        .build()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create thread pool: {}", e)))?;

    // Spawn processing workers
    let workers_header = header_arc.clone();
    let workers_result_sender = result_sender.clone();
    let workers_batch_receiver = batch_receiver.clone();
    
    let _workers_handle = std::thread::spawn(move || {
        pool.install(|| {
            while let Ok((batch_id, raw_records)) = workers_batch_receiver.recv() {
                let batch_start = Instant::now();
                
                // Process raw BAM records into Arrow RecordBatch
                let batch_result = process_bam_records_to_batch(
                    &raw_records,
                    &workers_header,
                    include_sequence,
                    include_quality,
                    batch_id,
                );
                
                match batch_result {
                    Ok(batch) => {
                        let batch_duration = batch_start.elapsed();
                        debug!("Thread processed batch {} ({} records) in {:?}", 
                                batch_id, batch.num_rows(), batch_duration);
                        
                        if let Err(_) = workers_result_sender.send((batch_id, batch)) {
                            debug!("Failed to send processed batch {}", batch_id);
                            break;
                        }
                    }
                    Err(e) => {
                        debug!("Failed to process batch {}: {}", batch_id, e);
                        break;
                    }
                }
            }
        });
    });

    // Writer thread with conditional ordering
    let writer_handle = std::thread::spawn(move || -> PyResult<usize> {
        let mut total_records = 0;
        
        if preserve_order {
            // Order preservation: Complex writer thread with HashMap buffering
            let mut received_batches: std::collections::HashMap<usize, RecordBatch> = std::collections::HashMap::new();
            let mut next_batch_id = 0;
            
            // Must wait for sequential batches
            while let Ok((batch_id, batch)) = result_receiver.recv() {
                received_batches.insert(batch_id, batch);
                
                // Write batches in order
                while let Some(batch) = received_batches.remove(&next_batch_id) {
                    arrow_writer.write(&batch)
                        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                            format!("Failed to write batch {}: {}", next_batch_id, e)
                        ))?;
                    
                    total_records += batch.num_rows();
                    next_batch_id += 1;
                    
                    // Reduce frequency: report every 200 batches instead of 50
                    if next_batch_id % 200 == 0 {
                        eprintln!("Progress: {} batches written, {} total records", next_batch_id, total_records);
                    }
                }
            }
            
            // Write any remaining batches
            for batch_id in next_batch_id.. {
                if let Some(batch) = received_batches.remove(&batch_id) {
                    arrow_writer.write(&batch)
                        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                            format!("Failed to write final batch {}: {}", batch_id, e)
                        ))?;
                    total_records += batch.num_rows();
                } else {
                    break;
                }
            }
        } else {
            // No order preservation: Direct write as batches arrive - no buffering needed!
            let mut batch_count = 0;
            
            while let Ok((_, batch)) = result_receiver.recv() {
                arrow_writer.write(&batch)
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                        format!("Failed to write batch: {}", e)
                    ))?;
                
                total_records += batch.num_rows();
                batch_count += 1;
                
                // Reduce frequency: report every 200 batches instead of 50
                if batch_count % 200 == 0 {
                    eprintln!("Progress: {} batches written, {} total records", batch_count, total_records);
                }
            }
        }
        
        arrow_writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Arrow IPC writer: {}", e)))?;
        
        Ok(total_records)
    });

    // Main thread: read BAM sequentially and distribute to workers
    let mut batch_id = 0;
    let mut current_batch: Vec<bam::Record> = Vec::with_capacity(effective_batch_size);
    let mut total_records_read = 0;
    let target_records = limit.unwrap_or(usize::MAX);
    
    loop {
        // Check for Python interrupts periodically
        if batch_id % 10 == 0 {
            Python::with_gil(|py| {
                py.check_signals().map_err(|e| {
                    eprintln!("Conversion interrupted by user");
                    e
                })
            })?;
        }

        let remaining_records = target_records.saturating_sub(total_records_read);
        if remaining_records == 0 {
            debug!("Reached limit of {} records", target_records);
            break;
        }

        let current_batch_size = effective_batch_size.min(remaining_records);
        
        // Read a batch of BAM records
        current_batch.clear();
        let mut record = bam::Record::default();
        
        for _ in 0..current_batch_size {
            match bam_reader.read_record(&mut record) {
                Ok(0) => {
                    // EOF reached
                    if !current_batch.is_empty() {
                        // Send the final partial batch
                        batch_sender.send((batch_id, current_batch.clone()))
                            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                                format!("Failed to send final batch: {}", e)
                            ))?;
                        batch_id += 1;
                    }
                    break;
                }
                Ok(_) => {
                    current_batch.push(record.clone());
                    total_records_read += 1;
                }
                Err(e) => {
                    return Err(PyErr::new::<PyRuntimeError, _>(
                        format!("Failed to read BAM record at position {}: {}", total_records_read, e)
                    ));
                }
            }
        }
        
        if current_batch.is_empty() {
            break; // EOF reached
        }
        
        debug!("Read batch {} with {} records", batch_id, current_batch.len());
        
        // Send batch to workers
        batch_sender.send((batch_id, current_batch.clone()))
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to send batch {}: {}", batch_id, e)
            ))?;
        
        batch_id += 1;
    }
    
    // Signal end of input
    drop(batch_sender);
    drop(result_sender);
    
    // Wait for writer to complete
    let final_record_count = writer_handle.join()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Writer thread failed: {:?}", e)))??;
    
    let total_duration = start_time.elapsed();
    let throughput = final_record_count as f64 / total_duration.as_secs_f64();
    
    let completion_msg = format!(
        "Parallel conversion complete: {} records written to {} in {:?} ({:.0} records/sec using {} threads)", 
        final_record_count, arrow_ipc_path, total_duration, throughput, effective_threads
    );
    eprintln!("{}", completion_msg);
    
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    arrow_ipc_path, 
    batch_size = 50000,
    include_sequence = true,
    include_quality = true,
    decompression_threads = 4,
    processing_threads = 2,
    preserve_order = false,
    limit = None
))]
pub fn bam_to_arrow_ipc_gzp_parallel(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    decompression_threads: usize,
    processing_threads: usize,
    preserve_order: bool,
    limit: Option<usize>,
) -> PyResult<()> {
    let start_time = Instant::now();
    
    // Validate thread parameters
    let effective_decompression_threads = if decompression_threads == 0 {
        eprintln!("Warning: decompression_threads cannot be 0, using default of 4 threads");
        4
    } else if decompression_threads > 16 {
        eprintln!("Notice: Capping decompression threads to 16 (requested: {})", decompression_threads);
        16
    } else {
        decompression_threads
    };
    
    let effective_processing_threads = if processing_threads == 0 {
        eprintln!("Warning: processing_threads cannot be 0, using default of 2 threads");
        2
    } else if processing_threads > 8 {
        eprintln!("Notice: Capping processing threads to 8 (requested: {})", processing_threads);
        8
    } else {
        processing_threads
    };
    
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
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    eprintln!("Starting enhanced parallel BAM to IPC conversion (gzp alternative):");
    eprintln!("  Input: {}", bam_path);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Batch size: {}", effective_batch_size);
    eprintln!("  Decompression threads: {} (applied as I/O optimization)", effective_decompression_threads);
    eprintln!("  Processing threads: {}", effective_processing_threads);
    eprintln!("  Preserve order: {}", preserve_order);
    eprintln!("  Include sequence: {}", include_sequence);
    eprintln!("  Include quality: {}", include_quality);

    // WORKAROUND: Since gzp doesn't support BGZF format, use enhanced buffering + parallel processing
    // This maintains the dual-thread architecture but optimizes I/O differently
    
    let file = File::open(input_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file '{}': {}", bam_path, e)))?;
    
    // Create a large buffer based on decompression thread count for better I/O
    let buffer_size = (effective_decompression_threads * 256 * 1024).max(1024 * 1024); // Min 1MB, scale with threads
    let buffered_file = std::io::BufReader::with_capacity(buffer_size, file);
    
    let mut bam_reader = bam::io::Reader::new(buffered_file);
    
    let header = bam_reader.read_header()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to read BAM header: {}", e)))?;

    let schema = create_bam_schema(include_sequence, include_quality);
    
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file '{}': {}", arrow_ipc_path, e)))?;

    let mut arrow_writer = ArrowIpcWriter::try_new(output_file, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow IPC writer: {}", e)))?;

    // Setup parallel processing architecture - similar to existing but specialized for gzp
    let channel_size = (effective_processing_threads * 4).max(16);
    let (batch_sender, batch_receiver): (Sender<(usize, Vec<bam::Record>)>, Receiver<(usize, Vec<bam::Record>)>) = bounded(channel_size);
    let (result_sender, result_receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    let header_arc = Arc::new(header.clone());
    
    // Configure processing thread pool (separate from decompression threads)
    let pool = ThreadPoolBuilder::new()
        .num_threads(effective_processing_threads)
        .build()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create processing thread pool: {}", e)))?;

    // Spawn processing workers (same logic as existing parallel implementation)
    let workers_header = header_arc.clone();
    let workers_result_sender = result_sender.clone();
    let workers_batch_receiver = batch_receiver.clone();
    
    let _workers_handle = std::thread::spawn(move || {
        pool.install(|| {
            while let Ok((batch_id, raw_records)) = workers_batch_receiver.recv() {
                let batch_start = Instant::now();
                
                // Process raw BAM records into Arrow RecordBatch
                let batch_result = process_bam_records_to_batch(
                    &raw_records,
                    &workers_header,
                    include_sequence,
                    include_quality,
                    batch_id,
                );
                
                match batch_result {
                    Ok(batch) => {
                        let batch_duration = batch_start.elapsed();
                        debug!("Thread processed batch {} ({} records) in {:?} with gzp decompression", 
                                batch_id, batch.num_rows(), batch_duration);
                        
                        if let Err(_) = workers_result_sender.send((batch_id, batch)) {
                            debug!("Failed to send processed batch {}", batch_id);
                            break;
                        }
                    }
                    Err(e) => {
                        debug!("Failed to process batch {}: {}", batch_id, e);
                        break;
                    }
                }
            }
        });
    });

    // Writer thread with conditional ordering (reuse existing logic)
    let writer_handle = std::thread::spawn(move || -> PyResult<usize> {
        let mut total_records = 0;
        
        if preserve_order {
            // Order preservation: Complex writer thread with HashMap buffering
            let mut received_batches: std::collections::HashMap<usize, RecordBatch> = std::collections::HashMap::new();
            let mut next_batch_id = 0;
            
            // Must wait for sequential batches
            while let Ok((batch_id, batch)) = result_receiver.recv() {
                received_batches.insert(batch_id, batch);
                
                // Write batches in order
                while let Some(batch) = received_batches.remove(&next_batch_id) {
                    arrow_writer.write(&batch)
                        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                            format!("Failed to write batch {}: {}", next_batch_id, e)
                        ))?;
                    
                    total_records += batch.num_rows();
                    next_batch_id += 1;
                    
                    // Reduce frequency: report every 200 batches instead of 50
                    if next_batch_id % 200 == 0 {
                        eprintln!("Progress: {} batches written, {} total records (gzp parallel)", next_batch_id, total_records);
                    }
                }
            }
            
            // Write any remaining batches
            for batch_id in next_batch_id.. {
                if let Some(batch) = received_batches.remove(&batch_id) {
                    arrow_writer.write(&batch)
                        .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                            format!("Failed to write final batch {}: {}", batch_id, e)
                        ))?;
                    total_records += batch.num_rows();
                } else {
                    break;
                }
            }
        } else {
            // No order preservation: Direct write as batches arrive - no buffering needed!
            let mut batch_count = 0;
            
            while let Ok((_, batch)) = result_receiver.recv() {
                arrow_writer.write(&batch)
                    .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                        format!("Failed to write batch: {}", e)
                    ))?;
                
                total_records += batch.num_rows();
                batch_count += 1;
                
                // Reduce frequency: report every 200 batches instead of 50
                if batch_count % 200 == 0 {
                    eprintln!("Progress: {} batches written, {} total records (gzp parallel)", batch_count, total_records);
                }
            }
        }
        
        arrow_writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Arrow IPC writer: {}", e)))?;
        
        Ok(total_records)
    });

    // Main thread: read BAM sequentially from gzp-decompressed stream and distribute to workers
    let mut batch_id = 0;
    let mut current_batch: Vec<bam::Record> = Vec::with_capacity(effective_batch_size);
    let mut total_records_read = 0;
    let target_records = limit.unwrap_or(usize::MAX);
    
    loop {
        // Check for Python interrupts periodically
        if batch_id % 10 == 0 {
            Python::with_gil(|py| {
                py.check_signals().map_err(|e| {
                    eprintln!("Conversion interrupted by user");
                    e
                })
            })?;
        }

        let remaining_records = target_records.saturating_sub(total_records_read);
        if remaining_records == 0 {
            debug!("Reached limit of {} records", target_records);
            break;
        }

        let current_batch_size = effective_batch_size.min(remaining_records);
        
        // Read a batch of BAM records from gzp-decompressed stream
        current_batch.clear();
        let mut record = bam::Record::default();
        
        for _ in 0..current_batch_size {
            match bam_reader.read_record(&mut record) {
                Ok(0) => {
                    // EOF reached
                    if !current_batch.is_empty() {
                        // Send the final partial batch
                        batch_sender.send((batch_id, current_batch.clone()))
                            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                                format!("Failed to send final batch: {}", e)
                            ))?;
                        batch_id += 1;
                    }
                    break;
                }
                Ok(_) => {
                    current_batch.push(record.clone());
                    total_records_read += 1;
                }
                Err(e) => {
                    return Err(PyErr::new::<PyRuntimeError, _>(
                        format!("Failed to read BAM record at position {} (gzp decompression): {}", total_records_read, e)
                    ));
                }
            }
        }
        
        if current_batch.is_empty() {
            break; // EOF reached
        }
        
        debug!("Read batch {} with {} records (gzp parallel)", batch_id, current_batch.len());
        
        // Send batch to workers
        batch_sender.send((batch_id, current_batch.clone()))
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to send batch {}: {}", batch_id, e)
            ))?;
        
        batch_id += 1;
    }
    
    // Signal end of input
    drop(batch_sender);
    drop(result_sender);
    
    // Wait for writer to complete
    let final_record_count = writer_handle.join()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Writer thread failed: {:?}", e)))??;
    
    let total_duration = start_time.elapsed();
    let throughput = final_record_count as f64 / total_duration.as_secs_f64();
    
    let buffer_mb = (effective_decompression_threads * 256).max(1024); // Recalculate for display
    let completion_msg = format!(
        "Enhanced parallel conversion complete: {} records written to {} in {:?} ({:.0} records/sec using {}MB buffer + {}+{} threads)", 
        final_record_count, arrow_ipc_path, total_duration, throughput, 
        buffer_mb / 1024, effective_processing_threads, 1
    );
    eprintln!("{}", completion_msg);
    
    Ok(())
}

// TEMPORARILY DISABLED: HTSlib parallel implementation - requires libclang system dependency
// TODO: Re-enable once libclang is properly installed in environment
#[cfg(feature = "htslib")]
#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    arrow_ipc_path, 
    batch_size = 15000,
    include_sequence = true,
    include_quality = true,
    bgzf_threads = 8,
    writing_threads = 8,
    read_buffer_mb = Some(1024),
    write_buffer_mb = Some(256),
    limit = None
))]
pub fn bam_to_arrow_ipc_htslib_parallel(
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
    let effective_bgzf_threads = if bgzf_threads == 0 { 4 } else { bgzf_threads.min(16) };
    let effective_writing_threads = writing_threads.max(1).min(16);
    
    // Independent buffer sizing (decoupled from thread count)
    let read_buffer_size = read_buffer_mb
        .unwrap_or(64)  // Default 64MB regardless of threads
        * 1024 * 1024;
        
    let write_buffer_size = write_buffer_mb
        .unwrap_or(32)  // Default 32MB per writer
        * 1024 * 1024;
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    eprintln!("Starting HTSlib parallel BGZF decompression BAM to IPC conversion:");
    eprintln!("  Input: {}", bam_path);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Batch size: {}", effective_batch_size);
    eprintln!("  BGZF threads: {}", effective_bgzf_threads);
    eprintln!("  Writing threads: {}", effective_writing_threads);
    eprintln!("  Read buffer: {}MB", read_buffer_size / (1024 * 1024));
    eprintln!("  Write buffer: {}MB per writer", write_buffer_size / (1024 * 1024));
    eprintln!("  Include sequence: {}", include_sequence);
    eprintln!("  Include quality: {}", include_quality);

    // Setup parallel processing architecture
    let schema = create_bam_schema(include_sequence, include_quality);
    let channel_size = (effective_writing_threads * 4).max(16);
    let (batch_sender, batch_receiver): (Sender<(usize, Vec<hts_bam::Record>)>, Receiver<(usize, Vec<hts_bam::Record>)>) = bounded(channel_size);
    let (result_sender, result_receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    // Configure processing thread pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(effective_writing_threads)
        .build()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create thread pool: {}", e)))?;

    // Writer thread with single buffered output stream
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file: {}", e)))?;
    let buffered_output = std::io::BufWriter::with_capacity(write_buffer_size, output_file);
    let mut arrow_writer = ArrowIpcWriter::try_new(buffered_output, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow writer: {}", e)))?;

    let writer_handle = std::thread::spawn(move || -> PyResult<usize> {
        let mut total_records = 0;
        let mut batch_count = 0;
        
        // Direct write as batches arrive (no order preservation for maximum speed)
        while let Ok((_, batch)) = result_receiver.recv() {
            arrow_writer.write(&batch)
                .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                    format!("Failed to write batch: {}", e)
                ))?;
            
            total_records += batch.num_rows();
            batch_count += 1;
            
            // Report progress every 200 batches
            if batch_count % 200 == 0 {
                eprintln!("Progress: {} batches written, {} total records (HTSlib parallel)", batch_count, total_records);
            }
        }
        
        arrow_writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Arrow writer: {}", e)))?;
        
        Ok(total_records)
    });

    // Main thread: read HTSlib records and distribute to workers
    let mut hts_reader = hts_bam::Reader::from_path(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file: {}", e)))?;
    
    // Enable parallel BGZF decompression
    hts_reader.set_threads(effective_bgzf_threads)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to set threads: {}", e)))?;

    // Build chromosome lookup table from header
    let header_view = hts_reader.header().clone();
    let chrom_lookup = build_chromosome_lookup(&header_view);

    // Spawn processing workers (moved here after chrom_lookup is available)
    let workers_result_sender = result_sender.clone();
    let workers_batch_receiver = batch_receiver.clone();
    let workers_chrom_lookup = chrom_lookup.clone();
    
    let _workers_handle = std::thread::spawn(move || {
        pool.install(|| {
            while let Ok((batch_id, hts_records)) = workers_batch_receiver.recv() {
                let batch_start = Instant::now();
                
                // Process HTSlib records into Arrow RecordBatch
                let batch_result = process_htslib_records_to_batch(
                    &hts_records,
                    &workers_chrom_lookup,
                    include_sequence,
                    include_quality,
                    batch_id,
                );
                
                match batch_result {
                    Ok(batch) => {
                        let batch_duration = batch_start.elapsed();
                        debug!("Thread processed HTSlib batch {} ({} records) in {:?}", 
                                batch_id, batch.num_rows(), batch_duration);
                        
                        if let Err(_) = workers_result_sender.send((batch_id, batch)) {
                            debug!("Failed to send processed HTSlib batch {}", batch_id);
                            break;
                        }
                    }
                    Err(e) => {
                        debug!("Failed to process HTSlib batch {}: {}", batch_id, e);
                        break;
                    }
                }
            }
        });
    });

    let target_records = limit.unwrap_or(usize::MAX);
    
    // Main reading loop: collect HTSlib records into batches and send to workers
    let mut batch_id = 0;
    let mut current_batch: Vec<hts_bam::Record> = Vec::with_capacity(effective_batch_size);
    let mut total_records_read = 0;
    let mut hts_record = hts_bam::Record::new();
    
    loop {
        // Check for Python interrupts periodically
        if batch_id % 10 == 0 {
            Python::with_gil(|py| {
                py.check_signals().map_err(|e| {
                    eprintln!("Conversion interrupted by user");
                    e
                })
            })?;
        }

        let remaining_records = target_records.saturating_sub(total_records_read);
        if remaining_records == 0 {
            debug!("Reached limit of {} records", target_records);
            break;
        }

        let current_batch_size = effective_batch_size.min(remaining_records);
        
        // Read a batch of HTSlib records
        current_batch.clear();
        
        for _ in 0..current_batch_size {
            match hts_reader.read(&mut hts_record) {
                Some(Ok(())) => {
                    current_batch.push(hts_record.clone());
                    total_records_read += 1;
                }
                None => {
                    // EOF reached
                    if !current_batch.is_empty() {
                        // Send the final partial batch
                        batch_sender.send((batch_id, current_batch.clone()))
                            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                                format!("Failed to send final batch: {}", e)
                            ))?;
                        batch_id += 1;
                    }
                    break;
                }
                Some(Err(e)) => {
                    return Err(PyErr::new::<PyRuntimeError, _>(
                        format!("Failed to read BAM record at position {}: {}", total_records_read, e)
                    ));
                }
            }
        }
        
        if current_batch.is_empty() {
            break; // EOF reached
        }
        
        debug!("Read HTSlib batch {} with {} records", batch_id, current_batch.len());
        
        // Send batch to workers
        batch_sender.send((batch_id, current_batch.clone()))
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to send batch {}: {}", batch_id, e)
            ))?;
        
        batch_id += 1;
    }
    
    // Signal end of input
    drop(batch_sender);
    drop(result_sender);
    
    // Wait for writer to complete
    let final_record_count = writer_handle.join()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Writer thread failed: {:?}", e)))??;
    
    let total_duration = start_time.elapsed();
    let throughput = final_record_count as f64 / total_duration.as_secs_f64();
    
    let completion_msg = format!(
        "HTSlib parallel conversion complete: {} records written to {} in {:?} ({:.0} records/sec using {} BGZF + {} processing threads)", 
        final_record_count, arrow_ipc_path, total_duration, throughput, effective_bgzf_threads, effective_writing_threads
    );
    eprintln!("{}", completion_msg);
    
    Ok(())
}

#[cfg(feature = "htslib")]
#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    arrow_ipc_path, 
    batch_size = 15000,
    include_sequence = true,
    include_quality = true,
    max_bgzf_threads = 16,
    writing_threads = 6,
    read_buffer_mb = Some(2048),
    write_buffer_mb = Some(512),
    limit = None
))]
pub fn bam_to_arrow_ipc_htslib_optimized(
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
    let effective_bgzf_threads = max_bgzf_threads.min(32); // Allow high thread counts
    let effective_writing_threads = writing_threads.max(1).min(16);
    
    // Much larger buffers for high-throughput
    let read_buffer_size = read_buffer_mb
        .unwrap_or(2048)  // Default 2GB read buffer
        * 1024 * 1024;
        
    let write_buffer_size = write_buffer_mb
        .unwrap_or(512)   // Default 512MB write buffer
        * 1024 * 1024;
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    eprintln!("Starting HTSlib OPTIMIZED high-throughput BAM to IPC conversion:");
    eprintln!("  Input: {}", bam_path);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Batch size: {}", effective_batch_size);
    eprintln!("  BGZF threads: {} (maximized)", effective_bgzf_threads);
    eprintln!("  Writing threads: {}", effective_writing_threads);
    eprintln!("  Read buffer: {}MB", read_buffer_size / (1024 * 1024));
    eprintln!("  Write buffer: {}MB", write_buffer_size / (1024 * 1024));
    eprintln!("  Include sequence: {}", include_sequence);
    eprintln!("  Include quality: {}", include_quality);

    // Setup parallel processing architecture (streamlined)
    let schema = create_bam_schema(include_sequence, include_quality);
    let channel_size = (effective_writing_threads * 8).max(32); // Larger queues
    let (batch_sender, batch_receiver): (Sender<(usize, Vec<hts_bam::Record>)>, Receiver<(usize, Vec<hts_bam::Record>)>) = bounded(channel_size);
    let (result_sender, result_receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    // Configure processing thread pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(effective_writing_threads)
        .build()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create thread pool: {}", e)))?;

    // Main thread: read HTSlib records with MAXIMUM threading
    let mut hts_reader = hts_bam::Reader::from_path(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file: {}", e)))?;
    
    // MAXIMIZE BGZF decompression threads - this is key!
    eprintln!("Setting BGZF threads to maximum: {}", effective_bgzf_threads);
    hts_reader.set_threads(effective_bgzf_threads)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to set threads: {}", e)))?;

    // Build chromosome lookup table from header
    let header_view = hts_reader.header().clone();
    let chrom_lookup = build_chromosome_lookup(&header_view);

    // Spawn processing workers with optimized chromosome lookup
    let workers_result_sender = result_sender.clone();
    let workers_batch_receiver = batch_receiver.clone();
    let workers_chrom_lookup = chrom_lookup.clone();
    
    let _workers_handle = std::thread::spawn(move || {
        pool.install(|| {
            while let Ok((batch_id, hts_records)) = workers_batch_receiver.recv() {
                let batch_start = Instant::now();
                
                // Process HTSlib records into Arrow RecordBatch
                let batch_result = process_htslib_records_to_batch(
                    &hts_records,
                    &workers_chrom_lookup,
                    include_sequence,
                    include_quality,
                    batch_id,
                );
                
                match batch_result {
                    Ok(batch) => {
                        let batch_duration = batch_start.elapsed();
                        debug!("Thread processed optimized batch {} ({} records) in {:?}", 
                                batch_id, batch.num_rows(), batch_duration);
                        
                        if let Err(_) = workers_result_sender.send((batch_id, batch)) {
                            debug!("Failed to send processed optimized batch {}", batch_id);
                            break;
                        }
                    }
                    Err(e) => {
                        debug!("Failed to process optimized batch {}: {}", batch_id, e);
                        break;
                    }
                }
            }
        });
    });

    // Writer thread with larger buffer
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file: {}", e)))?;
    let buffered_output = std::io::BufWriter::with_capacity(write_buffer_size, output_file);
    let mut arrow_writer = ArrowIpcWriter::try_new(buffered_output, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow writer: {}", e)))?;

    let writer_handle = std::thread::spawn(move || -> PyResult<usize> {
        let mut total_records = 0;
        let mut batch_count = 0;
        
        // Direct write as batches arrive (no order preservation for maximum speed)
        while let Ok((_, batch)) = result_receiver.recv() {
            arrow_writer.write(&batch)
                .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                    format!("Failed to write batch: {}", e)
                ))?;
            
            total_records += batch.num_rows();
            batch_count += 1;
            
            // Report progress every 100 batches
            if batch_count % 100 == 0 {
                eprintln!("Progress: {} batches written, {} total records (optimized)", batch_count, total_records);
            }
        }
        
        arrow_writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Arrow writer: {}", e)))?;
        
        Ok(total_records)
    });

    let target_records = limit.unwrap_or(usize::MAX);
    
    // Main reading loop: focus on feeding HTSlib's optimized decompression
    let mut batch_id = 0;
    let mut current_batch: Vec<hts_bam::Record> = Vec::with_capacity(effective_batch_size);
    let mut total_records_read = 0;
    let mut hts_record = hts_bam::Record::new();
    
    // Pre-allocate for efficiency
    current_batch.reserve(effective_batch_size);
    
    loop {
        // Check for Python interrupts periodically
        if batch_id % 10 == 0 {
            Python::with_gil(|py| {
                py.check_signals().map_err(|e| {
                    eprintln!("Conversion interrupted by user");
                    e
                })
            })?;
        }

        let remaining_records = target_records.saturating_sub(total_records_read);
        if remaining_records == 0 {
            debug!("Reached limit of {} records", target_records);
            break;
        }

        let current_batch_size = effective_batch_size.min(remaining_records);
        
        // Read a batch of HTSlib records - let HTSlib's threading do the work!
        current_batch.clear();
        
        for _ in 0..current_batch_size {
            match hts_reader.read(&mut hts_record) {
                Some(Ok(())) => {
                    current_batch.push(hts_record.clone());
                    total_records_read += 1;
                }
                None => {
                    // EOF reached
                    if !current_batch.is_empty() {
                        // Send the final partial batch
                        batch_sender.send((batch_id, current_batch.clone()))
                            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                                format!("Failed to send final batch: {}", e)
                            ))?;
                        batch_id += 1;
                    }
                    break;
                }
                Some(Err(e)) => {
                    return Err(PyErr::new::<PyRuntimeError, _>(
                        format!("Failed to read BAM record at position {}: {}", total_records_read, e)
                    ));
                }
            }
        }
        
        if current_batch.is_empty() {
            break; // EOF reached
        }
        
        debug!("Read optimized batch {} with {} records", batch_id, current_batch.len());
        
        // Send batch to workers
        batch_sender.send((batch_id, current_batch.clone()))
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to send batch {}: {}", batch_id, e)
            ))?;
        
        batch_id += 1;
    }
    
    // Signal end of input
    drop(batch_sender);
    drop(result_sender);
    
    // Wait for writer to complete
    let final_record_count = writer_handle.join()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Writer thread failed: {:?}", e)))??;
    
    let total_duration = start_time.elapsed();
    let throughput = final_record_count as f64 / total_duration.as_secs_f64();
    
    let completion_msg = format!(
        "HTSlib OPTIMIZED conversion complete: {} records written to {} in {:?} ({:.0} records/sec using {} BGZF + {} processing threads)", 
        final_record_count, arrow_ipc_path, total_duration, throughput, effective_bgzf_threads, effective_writing_threads
    );
    eprintln!("{}", completion_msg);
    
    Ok(())
}

#[cfg(feature = "htslib")]
#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    arrow_ipc_path, 
    batch_size = 15000,
    include_sequence = true,
    include_quality = true,
    num_workers = 4,
    chunk_size_mb = 64,
    bgzf_threads_per_worker = 2,
    limit = None
))]
pub fn bam_to_arrow_ipc_htslib_mmap_parallel(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    num_workers: usize,
    chunk_size_mb: usize,
    bgzf_threads_per_worker: usize,
    limit: Option<usize>,
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
    let effective_num_workers = num_workers.max(1).min(16);
    let _effective_chunk_size = (chunk_size_mb * 1024 * 1024).max(1024 * 1024); // Min 1MB
    let effective_bgzf_threads = bgzf_threads_per_worker.max(1).min(8);
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    eprintln!("Starting HTSlib MEMORY-MAPPED parallel BAM to IPC conversion:");
    eprintln!("  Input: {}", bam_path);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Batch size: {}", effective_batch_size);
    eprintln!("  Number of workers: {}", effective_num_workers);
    eprintln!("  Chunk size: {}MB", chunk_size_mb);
    eprintln!("  BGZF threads per worker: {}", effective_bgzf_threads);
    eprintln!("  Include sequence: {}", include_sequence);
    eprintln!("  Include quality: {}", include_quality);

    // Get file size for chunking strategy
    let file_size = std::fs::metadata(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to get file size: {}", e)))?
        .len();
    
    eprintln!("  File size: {:.1}MB", file_size as f64 / (1024.0 * 1024.0));

    // Create header reader to get chromosome lookup
    let temp_reader = hts_bam::Reader::from_path(bam_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to open BAM file: {}", e)))?;
    let header_view = temp_reader.header().clone();
    let chrom_lookup = build_chromosome_lookup(&header_view);

    // For now, implement a simplified approach:
    // Create multiple readers that read different portions of the file
    // This is a step toward true memory mapping
    
    // Setup channels
    let schema = create_bam_schema(include_sequence, include_quality);
    let channel_size = (effective_num_workers * 4).max(16);
    let (result_sender, result_receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    // Create worker threads that each process a different portion of the file
    let mut worker_handles = Vec::new();
    let target_records = limit.unwrap_or(usize::MAX);
    let records_per_worker = target_records / effective_num_workers;
    
    for worker_id in 0..effective_num_workers {
        let bam_path_clone = bam_path.to_string();
        let result_sender_clone = result_sender.clone();
        let chrom_lookup_clone = chrom_lookup.clone();
        let start_record = worker_id * records_per_worker;
        let end_record = if worker_id == effective_num_workers - 1 {
            target_records // Last worker takes remainder
        } else {
            (worker_id + 1) * records_per_worker
        };
        
        let handle = std::thread::spawn(move || {
            eprintln!("Worker {} processing records {} to {}", worker_id, start_record, end_record);
            
            // Each worker creates its own reader
            match hts_bam::Reader::from_path(&bam_path_clone) {
                Ok(mut worker_reader) => {
                    // Set threads for this worker
                    if let Err(e) = worker_reader.set_threads(effective_bgzf_threads) {
                        eprintln!("Worker {} failed to set threads: {}", worker_id, e);
                        return;
                    }
                    
                    // Simple approach: skip to start position
                    // TODO: Replace with true memory mapping + BGZF parsing
                    let mut skipped = 0;
                    let mut worker_record = hts_bam::Record::new();
                    
                    // Skip to start position
                    while skipped < start_record {
                        match worker_reader.read(&mut worker_record) {
                            Some(Ok(())) => skipped += 1,
                            Some(Err(e)) => {
                                eprintln!("Worker {} skip error: {}", worker_id, e);
                                return;
                            }
                            None => {
                                eprintln!("Worker {} reached EOF during skip", worker_id);
                                return;
                            }
                        }
                    }
                    
                    // Process records for this worker
                    let mut processed = 0;
                    let mut current_batch: Vec<hts_bam::Record> = Vec::with_capacity(effective_batch_size);
                    let mut batch_id = worker_id * 1000;
                    
                    while processed < (end_record - start_record) {
                        match worker_reader.read(&mut worker_record) {
                            Some(Ok(())) => {
                                current_batch.push(worker_record.clone());
                                processed += 1;
                                
                                if current_batch.len() >= effective_batch_size {
                                    // Process batch
                                    match process_htslib_records_to_batch(
                                        &current_batch,
                                        &chrom_lookup_clone,
                                        include_sequence,
                                        include_quality,
                                        batch_id,
                                    ) {
                                        Ok(batch) => {
                                            if let Err(_) = result_sender_clone.send((batch_id, batch)) {
                                                eprintln!("Worker {} failed to send batch {}", worker_id, batch_id);
                                                return;
                                            }
                                        }
                                        Err(e) => {
                                            eprintln!("Worker {} batch processing error: {}", worker_id, e);
                                            return;
                                        }
                                    }
                                    current_batch.clear();
                                    batch_id += 1;
                                }
                            }
                            Some(Err(e)) => {
                                eprintln!("Worker {} read error: {}", worker_id, e);
                                break;
                            }
                            None => {
                                eprintln!("Worker {} reached EOF", worker_id);
                                break;
                            }
                        }
                    }
                    
                    // Process final batch
                    if !current_batch.is_empty() {
                        match process_htslib_records_to_batch(
                            &current_batch,
                            &chrom_lookup_clone,
                            include_sequence,
                            include_quality,
                            batch_id,
                        ) {
                            Ok(batch) => {
                                let _ = result_sender_clone.send((batch_id, batch));
                            }
                            Err(e) => {
                                eprintln!("Worker {} final batch error: {}", worker_id, e);
                            }
                        }
                    }
                    
                    eprintln!("Worker {} completed: {} records processed", worker_id, processed);
                }
                Err(e) => {
                    eprintln!("Worker {} failed to open BAM file: {}", worker_id, e);
                }
            }
        });
        
        worker_handles.push(handle);
    }
    
    // Writer thread
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file: {}", e)))?;
    let mut arrow_writer = ArrowIpcWriter::try_new(output_file, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow writer: {}", e)))?;

    let writer_handle = std::thread::spawn(move || -> PyResult<usize> {
        let mut total_records = 0;
        let mut batch_count = 0;
        
        while let Ok((_, batch)) = result_receiver.recv() {
            arrow_writer.write(&batch)
                .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                    format!("Failed to write batch: {}", e)
                ))?;
            
            total_records += batch.num_rows();
            batch_count += 1;
            
            if batch_count % 50 == 0 {
                eprintln!("Progress: {} batches written, {} total records (mmap)", batch_count, total_records);
            }
        }
        
        arrow_writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Arrow writer: {}", e)))?;
        
        Ok(total_records)
    });
    
    // Wait for all workers to complete
    for (i, handle) in worker_handles.into_iter().enumerate() {
        if let Err(e) = handle.join() {
            eprintln!("Worker {} thread panicked: {:?}", i, e);
        }
    }
    
    // Signal end of input
    drop(result_sender);
    
    // Wait for writer to complete
    let final_record_count = writer_handle.join()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Writer thread failed: {:?}", e)))??;
    
    let total_duration = start_time.elapsed();
    let throughput = final_record_count as f64 / total_duration.as_secs_f64();
    
    let completion_msg = format!(
        "Memory-mapped parallel conversion complete: {} records written to {} in {:?} ({:.0} records/sec using {} workers)", 
        final_record_count, arrow_ipc_path, total_duration, throughput, effective_num_workers
    );
    eprintln!("{}", completion_msg);
    
    Ok(())
}

// Helper function to build chromosome lookup table
#[cfg(feature = "htslib")]
fn build_chromosome_lookup(header: &hts_bam::HeaderView) -> Arc<Vec<String>> {
    let mut lookup = Vec::new();
    for tid in 0..header.target_count() {
        let chrom_name = String::from_utf8_lossy(header.tid2name(tid as u32)).into_owned();
        lookup.push(chrom_name);
    }
    Arc::new(lookup)
}

// Zero-copy quality score processing helper
#[cfg(feature = "htslib")]
fn quality_to_string_zero_copy(qual: &[u8]) -> String {
    if qual.is_empty() {
        return String::new();
    }
    
    let mut result = String::with_capacity(qual.len());
    unsafe {
        let bytes = result.as_mut_vec();
        bytes.extend_from_slice(qual);           // Single memcpy operation
        for byte in bytes {
            *byte += b'!';                       // In-place addition (PHRED+33)
        }
    }
    result
}

// Multi-reader seek-based functionality
#[cfg(feature = "htslib")]
fn test_multi_position_reading(bam_path: &str) -> Result<bool, Box<dyn std::error::Error + Send + Sync>> {
    use std::fs;
    
    let file_size = fs::metadata(bam_path)?.len();
    let _mid_point = file_size / 2;
    
    // Test if we can create multiple readers and seek to different positions
    let mut reader1 = hts_bam::Reader::from_path(bam_path)?;
    let mut reader2 = hts_bam::Reader::from_path(bam_path)?;
    
    // Try to seek reader2 to middle of file
    let _seek_result = reader2.set_threads(1); // Ensure single-threaded for seeking
    
    // For now, let's test basic dual-reader creation
    // More sophisticated seeking will be implemented next
    let mut record1 = hts_bam::Record::new();
    let mut record2 = hts_bam::Record::new();
    
    // Test that both readers can read independently
    let read1_ok = reader1.read(&mut record1).is_some();
    let read2_ok = reader2.read(&mut record2).is_some();
    
    Ok(read1_ok && read2_ok)
}

#[cfg(feature = "htslib")]
fn discover_safe_seek_points(bam_path: &str, num_segments: usize) -> Result<Vec<u64>, Box<dyn std::error::Error + Send + Sync>> {
    use std::fs;
    
    let file_size = fs::metadata(bam_path)?.len();
    let segment_size = file_size / num_segments as u64;
    
    let mut seek_points = Vec::new();
    seek_points.push(0); // Always start from beginning
    
    for i in 1..num_segments {
        let approximate_pos = i as u64 * segment_size;
        
        // For initial implementation, use approximate positions
        // TODO: Implement proper record boundary detection
        seek_points.push(approximate_pos);
    }
    
    Ok(seek_points)
}

// Multi-reader pipeline architecture
#[cfg(feature = "htslib")]
struct MultiReaderPipeline {
    bam_path: String,
    num_readers: usize,
    seek_points: Vec<u64>,
    chrom_lookup: Arc<Vec<String>>,
}

#[cfg(feature = "htslib")]
impl MultiReaderPipeline {
    fn new(bam_path: &str, num_readers: usize) -> Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        // Create initial reader to get header and discover seek points
        let temp_reader = hts_bam::Reader::from_path(bam_path)?;
        let header_view = temp_reader.header().clone();
        let chrom_lookup = build_chromosome_lookup(&header_view);
        
        let seek_points = discover_safe_seek_points(bam_path, num_readers)?;
        
        Ok(MultiReaderPipeline {
            bam_path: bam_path.to_string(),
            num_readers,
            seek_points,
            chrom_lookup,
        })
    }
    
    fn spawn_reader_threads(
        &self,
        batch_size: usize,
        batch_sender: Sender<(usize, Vec<hts_bam::Record>)>,
        target_records: usize,
    ) -> Vec<std::thread::JoinHandle<()>> {
        let mut handles = Vec::new();
        
        for reader_id in 0..self.num_readers {
            let bam_path = self.bam_path.clone();
            let seek_start = self.seek_points[reader_id];
            let seek_end = self.seek_points.get(reader_id + 1).copied().unwrap_or(u64::MAX);
            let batch_sender_clone = batch_sender.clone();
            let records_per_reader = target_records / self.num_readers;
            
            let handle = std::thread::spawn(move || {
                debug!("Reader {} starting: seek_start={}, seek_end={}, target_records={}", 
                       reader_id, seek_start, seek_end, records_per_reader);
                
                if let Err(e) = read_segment(
                    &bam_path,
                    reader_id,
                    seek_start,
                    seek_end,
                    batch_size,
                    records_per_reader,
                    batch_sender_clone,
                ) {
                    debug!("Reader {} failed: {}", reader_id, e);
                }
            });
            
            handles.push(handle);
        }
        
        handles
    }
}

#[cfg(feature = "htslib")]
fn read_segment(
    bam_path: &str,
    reader_id: usize,
    seek_start: u64,
    _seek_end: u64,
    batch_size: usize,
    max_records: usize,
    batch_sender: Sender<(usize, Vec<hts_bam::Record>)>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut reader = hts_bam::Reader::from_path(bam_path)?;
    reader.set_threads(1)?; // Single-threaded for seeking
    
    // For initial implementation, only seek if not the first reader
    if reader_id > 0 && seek_start > 0 {
        // TODO: Implement proper BGZF seeking
        // For now, skip seeking and let all readers read from start (testing concurrent access)
        debug!("Reader {} would seek to position {}", reader_id, seek_start);
    }
    
    let mut current_batch: Vec<hts_bam::Record> = Vec::with_capacity(batch_size);
    let mut total_records_read = 0;
    let mut batch_id = reader_id * 1000; // Offset batch IDs by reader
    let mut record = hts_bam::Record::new();
    
    // Skip records for readers other than the first (simple load distribution)
    let skip_count = reader_id * (max_records / 4); // Simple skip strategy
    for _ in 0..skip_count {
        if reader.read(&mut record).is_none() {
            break;
        }
    }
    
    loop {
        if total_records_read >= max_records {
            break;
        }
        
        match reader.read(&mut record) {
            Some(Ok(())) => {
                current_batch.push(record.clone());
                total_records_read += 1;
                
                if current_batch.len() >= batch_size {
                    if let Err(_) = batch_sender.send((batch_id, current_batch.clone())) {
                        debug!("Reader {} failed to send batch {}", reader_id, batch_id);
                        break;
                    }
                    current_batch.clear();
                    batch_id += 1;
                }
            }
            Some(Err(e)) => {
                debug!("Reader {} read error: {}", reader_id, e);
                break;
            }
            None => {
                debug!("Reader {} reached EOF", reader_id);
                break;
            }
        }
    }
    
    // Send final partial batch
    if !current_batch.is_empty() {
        let _ = batch_sender.send((batch_id, current_batch));
    }
    
    debug!("Reader {} completed: {} records read", reader_id, total_records_read);
    Ok(())
}

#[cfg(feature = "htslib")]
#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    arrow_ipc_path, 
    batch_size = 15000,
    include_sequence = true,
    include_quality = true,
    num_readers = 2,
    bgzf_threads = 4,
    writing_threads = 10,
    read_buffer_mb = Some(1024),
    write_buffer_mb = Some(256),
    limit = None
))]
pub fn bam_to_arrow_ipc_htslib_multi_reader_parallel(
    bam_path: &str,
    arrow_ipc_path: &str,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
    num_readers: usize,
    bgzf_threads: usize,
    writing_threads: usize,
    read_buffer_mb: Option<usize>,
    write_buffer_mb: Option<usize>,
    limit: Option<usize>,
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
    let effective_num_readers = num_readers.max(1).min(8); // Limit to reasonable range
    let effective_bgzf_threads = if bgzf_threads == 0 { 4 } else { bgzf_threads.min(16) };
    let effective_writing_threads = writing_threads.max(1).min(16);
    
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                format!("Failed to create output directory: {}", e)
            ))?;
    }

    eprintln!("Starting HTSlib MULTI-READER parallel BAM to IPC conversion:");
    eprintln!("  Input: {}", bam_path);
    eprintln!("  Output: {}", arrow_ipc_path);
    eprintln!("  Batch size: {}", effective_batch_size);
    eprintln!("  Number of readers: {}", effective_num_readers);
    eprintln!("  BGZF threads per reader: {}", effective_bgzf_threads);
    eprintln!("  Writing threads: {}", effective_writing_threads);
    eprintln!("  Include sequence: {}", include_sequence);
    eprintln!("  Include quality: {}", include_quality);

    // Test feasibility first
    match test_multi_position_reading(bam_path) {
        Ok(true) => eprintln!("✅ Multi-reader feasibility test passed"),
        Ok(false) => {
            eprintln!("⚠️  Multi-reader test failed, falling back to single reader");
            return bam_to_arrow_ipc_htslib_parallel(
                bam_path, arrow_ipc_path, batch_size, include_sequence, include_quality,
                bgzf_threads, writing_threads, read_buffer_mb, write_buffer_mb, limit
            );
        }
        Err(e) => {
            eprintln!("❌ Multi-reader test error: {}, falling back to single reader", e);
            return bam_to_arrow_ipc_htslib_parallel(
                bam_path, arrow_ipc_path, batch_size, include_sequence, include_quality,
                bgzf_threads, writing_threads, read_buffer_mb, write_buffer_mb, limit
            );
        }
    }

    // Setup multi-reader pipeline
    let pipeline = MultiReaderPipeline::new(bam_path, effective_num_readers)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create multi-reader pipeline: {}", e)))?;

    // Setup channels and processing architecture
    let schema = create_bam_schema(include_sequence, include_quality);
    let channel_size = (effective_writing_threads * 4).max(16);
    let (batch_sender, batch_receiver): (Sender<(usize, Vec<hts_bam::Record>)>, Receiver<(usize, Vec<hts_bam::Record>)>) = bounded(channel_size);
    let (result_sender, result_receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    // Configure processing thread pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(effective_writing_threads)
        .build()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create thread pool: {}", e)))?;

    // Spawn processing workers
    let workers_result_sender = result_sender.clone();
    let workers_batch_receiver = batch_receiver.clone();
    let workers_chrom_lookup = pipeline.chrom_lookup.clone();
    
    let _workers_handle = std::thread::spawn(move || {
        pool.install(|| {
            while let Ok((batch_id, hts_records)) = workers_batch_receiver.recv() {
                let batch_start = Instant::now();
                
                let batch_result = process_htslib_records_to_batch(
                    &hts_records,
                    &workers_chrom_lookup,
                    include_sequence,
                    include_quality,
                    batch_id,
                );
                
                match batch_result {
                    Ok(batch) => {
                        let batch_duration = batch_start.elapsed();
                        debug!("Thread processed multi-reader batch {} ({} records) in {:?}", 
                                batch_id, batch.num_rows(), batch_duration);
                        
                        if let Err(_) = workers_result_sender.send((batch_id, batch)) {
                            debug!("Failed to send processed multi-reader batch {}", batch_id);
                            break;
                        }
                    }
                    Err(e) => {
                        debug!("Failed to process multi-reader batch {}: {}", batch_id, e);
                        break;
                    }
                }
            }
        });
    });

    // Writer thread
    let write_buffer_size = write_buffer_mb.unwrap_or(256) * 1024 * 1024;
    let output_file = File::create(output_path)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create output file: {}", e)))?;
    let buffered_output = std::io::BufWriter::with_capacity(write_buffer_size, output_file);
    let mut arrow_writer = ArrowIpcWriter::try_new(buffered_output, &schema)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create Arrow writer: {}", e)))?;

    let writer_handle = std::thread::spawn(move || -> PyResult<usize> {
        let mut total_records = 0;
        let mut batch_count = 0;
        
        while let Ok((_, batch)) = result_receiver.recv() {
            arrow_writer.write(&batch)
                .map_err(|e| PyErr::new::<PyRuntimeError, _>(
                    format!("Failed to write batch: {}", e)
                ))?;
            
            total_records += batch.num_rows();
            batch_count += 1;
            
            if batch_count % 100 == 0 {
                eprintln!("Progress: {} batches written, {} total records (multi-reader)", batch_count, total_records);
            }
        }
        
        arrow_writer.finish()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to close Arrow writer: {}", e)))?;
        
        Ok(total_records)
    });

    // Start multi-reader threads
    let target_records = limit.unwrap_or(usize::MAX);
    let reader_handles = pipeline.spawn_reader_threads(effective_batch_size, batch_sender, target_records);
    
    // Wait for all readers to complete
    for (i, handle) in reader_handles.into_iter().enumerate() {
        if let Err(e) = handle.join() {
            eprintln!("Reader {} thread panicked: {:?}", i, e);
        }
    }
    
    // Signal end of input
    drop(result_sender);
    
    // Wait for writer to complete
    let final_record_count = writer_handle.join()
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Writer thread failed: {:?}", e)))??;
    
    let total_duration = start_time.elapsed();
    let throughput = final_record_count as f64 / total_duration.as_secs_f64();
    
    let completion_msg = format!(
        "Multi-reader parallel conversion complete: {} records written to {} in {:?} ({:.0} records/sec using {} readers + {} processing threads)", 
        final_record_count, arrow_ipc_path, total_duration, throughput, effective_num_readers, effective_writing_threads
    );
    eprintln!("{}", completion_msg);
    
    Ok(())
}

// Helper function to process HTSlib records into Arrow RecordBatch
#[cfg(feature = "htslib")]
fn process_htslib_records_to_batch(
    hts_records: &[hts_bam::Record],
    chrom_lookup: &Arc<Vec<String>>,
    include_sequence: bool,
    include_quality: bool,
    batch_id: usize,
) -> Result<RecordBatch, Box<dyn std::error::Error + Send + Sync>> {
    // Pre-allocate vectors for Arrow arrays
    let mut names = Vec::with_capacity(hts_records.len());
    let mut chroms = Vec::with_capacity(hts_records.len());
    let mut starts = Vec::with_capacity(hts_records.len());
    let mut ends = Vec::with_capacity(hts_records.len());
    let mut flags = Vec::with_capacity(hts_records.len());
    let mut sequences: Option<Vec<Option<String>>> = if include_sequence { Some(Vec::with_capacity(hts_records.len())) } else { None };
    let mut qualities: Option<Vec<Option<String>>> = if include_quality { Some(Vec::with_capacity(hts_records.len())) } else { None };
    
    for hts_record in hts_records {
        // 1. Read name
        let name = String::from_utf8_lossy(hts_record.qname()).to_string();
        names.push(name);
        
        // 2. Chromosome - use tid to resolve reference name
        let chrom = if hts_record.tid() >= 0 {
            let tid = hts_record.tid() as usize;
            if tid < chrom_lookup.len() {
                Some(chrom_lookup[tid].clone())
            } else {
                None  // Invalid tid
            }
        } else {
            None  // Unmapped read
        };
        chroms.push(chrom);
        
        // 3. Start/End positions (convert from 0-based to 1-based)
        let start_pos = if hts_record.pos() >= 0 {
            Some((hts_record.pos() + 1) as u32)
        } else {
            None
        };
        starts.push(start_pos);
        
        // Calculate end position from start + read length
        let end_pos = start_pos.map(|start| start + hts_record.seq_len() as u32 - 1);
        ends.push(end_pos);
        
        // 4. Flags
        flags.push(hts_record.flags() as u32);
        
        // 5. Sequence (if requested)
        if let Some(ref mut seq_vec) = sequences {
            let seq = hts_record.seq();
            if seq.len() > 0 {
                let mut sequence_string = String::with_capacity(seq.len());
                for i in 0..seq.len() {
                    let base = match seq[i] {
                        1 => 'A',  // HTSlib encoding: A=1, C=2, G=4, T=8
                        2 => 'C',
                        4 => 'G', 
                        8 => 'T',
                        15 => 'N', // 15 = all bits set
                        _ => 'N',
                    };
                    sequence_string.push(base);
                }
                seq_vec.push(Some(sequence_string));
            } else {
                seq_vec.push(None);
            }
        }
        
        // 6. Quality scores (if requested) - ZERO-COPY OPTIMIZED
        if let Some(ref mut qual_vec) = qualities {
            let qual = hts_record.qual();
            if !qual.is_empty() && qual[0] != 255 {
                let quality_string = quality_to_string_zero_copy(qual);
                qual_vec.push(Some(quality_string));
            } else {
                qual_vec.push(None);
            }
        }
    }
    
    // Build Arrow arrays
    let name_array = Arc::new(StringArray::from(names));
    let chrom_array = Arc::new(StringArray::from(chroms));
    let start_array = Arc::new(UInt32Array::from(starts));
    let end_array = Arc::new(UInt32Array::from(ends));
    let flag_array = Arc::new(UInt32Array::from(flags));

    // Build column arrays
    let mut columns: Vec<Arc<dyn Array>> = vec![
        name_array,
        chrom_array,
        start_array,
        end_array,
        flag_array,
    ];

    // Add sequence if requested
    if include_sequence {
        if let Some(seq_vec) = sequences {
            let seq_array = Arc::new(StringArray::from(seq_vec));
            columns.push(seq_array);
        }
    }

    // Add quality if requested
    if include_quality {
        if let Some(qual_vec) = qualities {
            let qual_array = Arc::new(StringArray::from(qual_vec));
            columns.push(qual_array);
        }
    }

    // Create RecordBatch
    let schema = create_bam_schema(include_sequence, include_quality);
    RecordBatch::try_new(schema, columns)
        .map_err(|e| format!("Failed to create RecordBatch for HTSlib batch {}: {}", batch_id, e).into())
}


// Helper function to process raw BAM records into Arrow RecordBatch
fn process_bam_records_to_batch(
    records: &[bam::Record],
    header: &sam::Header,
    include_sequence: bool,
    include_quality: bool,
    batch_id: usize,
) -> Result<RecordBatch, Box<dyn std::error::Error + Send + Sync>> {
    let mut buffers = ReusableBuffers::new(records.len());
    
    for record in records {
        extract_record_data_enhanced(
            record,
            header,
            &mut buffers,
            include_sequence,
            include_quality,
        ).map_err(|e| format!("Failed to extract record data in batch {}: {}", batch_id, e))?;
    }
    
    // Create Arrow arrays from buffers
    let name_array = Arc::new(StringArray::from(buffers.names.clone()));
    let chrom_array = Arc::new(StringArray::from(buffers.chroms.clone()));
    let start_array = Arc::new(UInt32Array::from(buffers.starts.clone()));
    let end_array = Arc::new(UInt32Array::from(buffers.ends.clone()));
    let flags_array = Arc::new(UInt32Array::from(buffers.flags.clone()));
    
    let mut arrays: Vec<ArrayRef> = vec![
        name_array,
        chrom_array, 
        start_array,
        end_array,
        flags_array,
    ];
    
    if include_sequence {
        if let Some(ref sequences) = buffers.sequences {
            arrays.push(Arc::new(StringArray::from(sequences.clone())));
        }
    }
    
    if include_quality {
        if let Some(ref quality_scores) = buffers.quality_scores {
            arrays.push(Arc::new(StringArray::from(quality_scores.clone())));
        }
    }

    let schema = create_bam_schema(include_sequence, include_quality);
    RecordBatch::try_new(schema, arrays)
        .map_err(|e| format!("Failed to create record batch for batch {}: {}", batch_id, e).into())
}


fn create_bam_schema(include_sequence: bool, include_quality: bool) -> Arc<Schema> {
    let mut fields = vec![
        Field::new("name", DataType::Utf8, false),        // Read ID - never null
        Field::new("chrom", DataType::Utf8, true),        // Chromosome - can be null for unmapped
        Field::new("start", DataType::UInt32, true),      // 1-based start position
        Field::new("end", DataType::UInt32, true),        // 1-based end position  
        Field::new("flags", DataType::UInt32, false),     // SAM flags - never null
    ];

    if include_sequence {
        fields.push(Field::new("sequence", DataType::Utf8, true));
    }
    
    if include_quality {
        fields.push(Field::new("quality_scores", DataType::Utf8, true));
    }

    Arc::new(Schema::new(fields))
}





fn decode_base(encoded: u8) -> char {
    match encoded {
        1 => 'A',
        2 => 'C',
        4 => 'G',
        8 => 'T',
        15 => 'N', 
        _ => 'N',  
    }
}

fn calculate_bam_alignment_length(cigar: &bam::record::Cigar) -> u32 {
    let mut length = 0u32;
    for op_result in cigar.iter() {
        if let Ok(op) = op_result {
            match op.kind() {
                CigarOpKind::Match |
                CigarOpKind::SequenceMatch |
                CigarOpKind::SequenceMismatch |
                CigarOpKind::Deletion |
                CigarOpKind::Skip => {
                    length += op.len() as u32;
                }
                _ => {} // Insertions, soft clips, etc. don't contribute to reference length
            }
        }
        // Skip invalid CIGAR operations
    }
    length
}



// TODO: Add these columns later if needed:
// - cigar: String              // Alignment details  
// - mapping_quality: UInt8     // Alignment confidence (0-255)
// - mate_chrom: String         // Mate chromosome for paired reads
// - mate_start: UInt32         // Mate start position
// - template_length: Int32     // Insert size for paired reads

fn _future_expanded_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        // Current fields
        Field::new("name", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, true),
        Field::new("start", DataType::UInt32, true),
        Field::new("end", DataType::UInt32, true),
        Field::new("flags", DataType::UInt32, false),
        Field::new("sequence", DataType::Utf8, true),
        Field::new("quality_scores", DataType::Utf8, true),
        
        // Future fields ??
        // Field::new("cigar", DataType::Utf8, true),
        // Field::new("mapping_quality", DataType::UInt8, true),
        // Field::new("mate_chrom", DataType::Utf8, true),
        // Field::new("mate_start", DataType::UInt32, true),
        // Field::new("template_length", DataType::Int32, true),
    ]))
}

fn parse_compression(compression: &str) -> Compression {
    match compression.to_lowercase().as_str() {
        "snappy" => Compression::SNAPPY,
        "gzip" => Compression::GZIP(Default::default()),   
        "lz4" => Compression::LZ4,
        "brotli" => Compression::BROTLI(Default::default()), 
        "zstd" => Compression::ZSTD(Default::default()),   
        "uncompressed" | "none" => Compression::UNCOMPRESSED,
        _ => {
            eprintln!("Unknown compression '{}', defaulting to snappy", compression);
            Compression::SNAPPY
        }
    }
}
