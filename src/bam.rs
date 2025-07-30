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
    batch_size = 20000,
    include_sequence = true,
    include_quality = true,
    bgzf_threads = 4,
    writing_threads = 12,
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

    // Spawn processing workers
    let workers_result_sender = result_sender.clone();
    let workers_batch_receiver = batch_receiver.clone();
    
    let _workers_handle = std::thread::spawn(move || {
        pool.install(|| {
            while let Ok((batch_id, hts_records)) = workers_batch_receiver.recv() {
                let batch_start = Instant::now();
                
                // Process HTSlib records into Arrow RecordBatch
                let batch_result = process_htslib_records_to_batch(
                    &hts_records,
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


// Helper function to process HTSlib records into Arrow RecordBatch
#[cfg(feature = "htslib")]
fn process_htslib_records_to_batch(
    hts_records: &[hts_bam::Record],
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
            // We need the header to resolve chromosome names, but we can't access it here
            // For now, use the tid directly as string - this is a limitation we need to fix
            Some(format!("chr_{}", hts_record.tid()))
        } else {
            None
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
        
        // 6. Quality scores (if requested)
        if let Some(ref mut qual_vec) = qualities {
            let qual = hts_record.qual();
            if !qual.is_empty() && qual[0] != 255 {
                let quality_string: String = qual.iter()
                    .map(|&q| char::from(q + b'!'))
                    .collect();
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
