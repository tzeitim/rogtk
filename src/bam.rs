use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use pyo3::Python;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use noodles::{bam, sam, bgzf};
use sam::alignment::record::cigar::op::Kind as CigarOpKind;
use sam::alignment::record::QualityScores; 
use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, Encoding};
use parquet::file::properties::WriterProperties;


#[pyfunction]
#[pyo3(signature = (
    bam_path, 
    parquet_path, 
    batch_size = 100000,  // Increased default batch size for better throughput
    include_sequence = true,
    include_quality = true,
    compression = "snappy",
    limit = None
))]
pub fn bam_to_parquet_improved(
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
        
        let batch = read_bam_batch_improved(
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
        
        // More frequent progress updates with memory info
        if batch_count % 50 == 0 {
            let progress_msg = if let Some(limit_val) = limit {
                format!("Processed {} batches, {} / {} records ({:.1}%) - Batch size: {}", 
                       batch_count, total_records, limit_val, 
                       100.0 * total_records as f64 / limit_val as f64,
                       current_batch_size)
            } else {
                format!("Processed {} batches, {} total records - Batch size: {}", 
                       batch_count, total_records, current_batch_size)
            };
            eprintln!("{}", progress_msg);
            
            // Force garbage collection every 100 batches
            if batch_count % 100 == 0 {
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

fn read_bam_batch_improved(
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
                extract_record_data_improved(
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

fn extract_record_data_improved(
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
    batch_size = 10000,
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
        let current_batch_size = batch_size.min(remaining_records);
        
        let batch = read_bam_batch(
            &mut bam_reader, 
            &header, 
            current_batch_size, 
            include_sequence, 
            include_quality
        )?;
        
        if batch.num_rows() == 0 {
            // no more batches to go; we are done
            break;
        }

        parquet_writer.write(&batch)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to write batch {}: {}", batch_count, e)))?;

        total_records += batch.num_rows();
        batch_count += 1;
        
        if batch_count % 100 == 0 {
            let progress_msg = if let Some(limit_val) = limit {
                format!("Processed {} batches, {} / {} records ({:.1}%)", 
                       batch_count, total_records, limit_val, 
                       100.0 * total_records as f64 / limit_val as f64)
            } else {
                format!("Processed {} batches, {} total records", batch_count, total_records)
            };
            eprintln!("{}", progress_msg);
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


fn read_bam_batch(
    reader: &mut bam::io::Reader<bgzf::Reader<&mut File>>,
    header: &sam::Header,
    batch_size: usize,
    include_sequence: bool,
    include_quality: bool,
) -> PyResult<RecordBatch> {
    let mut names = Vec::with_capacity(batch_size);
    let mut chroms = Vec::with_capacity(batch_size);
    let mut starts = Vec::with_capacity(batch_size);
    let mut ends = Vec::with_capacity(batch_size);
    let mut flags = Vec::with_capacity(batch_size);
    
    let mut sequences = if include_sequence { 
        Some(Vec::with_capacity(batch_size)) 
    } else { 
        None 
    };
    let mut quality_scores = if include_quality { 
        Some(Vec::with_capacity(batch_size)) 
    } else { 
        None 
    };

    let mut count = 0;
    let mut record = bam::Record::default();

    while count < batch_size {
        match reader.read_record(&mut record) {
            Ok(0) => break, // EOF - no more records
            Ok(_) => {
                extract_record_data(
                    &record, 
                    header, 
                    &mut names, 
                    &mut chroms, 
                    &mut starts, 
                    &mut ends, 
                    &mut flags,
                    &mut sequences,
                    &mut quality_scores
                )?;
                count += 1;
            }
            Err(e) => {
                return Err(PyErr::new::<PyRuntimeError, _>(
                    format!("Error reading BAM record {}: {}", count, e)
                ));
            }
        }
    }

    create_record_batch(names, chroms, starts, ends, flags, sequences, quality_scores, include_sequence, include_quality)
}


fn extract_record_data(
    record: &bam::Record,
    header: &sam::Header,
    names: &mut Vec<String>,
    chroms: &mut Vec<Option<String>>,
    starts: &mut Vec<Option<u32>>,
    ends: &mut Vec<Option<u32>>,
    flags: &mut Vec<u32>,
    sequences: &mut Option<Vec<Option<String>>>,
    quality_scores: &mut Option<Vec<Option<String>>>,
) -> PyResult<()> {
    // 1. Extract read name
    let name = record.name()
        .map(|bstr| bstr.to_string())
        .unwrap_or_else(|| "unknown".to_string());
    names.push(name);

    // 2. Extract chromosome/reference name
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
    chroms.push(chrom);

    // 3. Extract coordinates (convert to 1-based to match polars-bio)
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
    starts.push(start_pos);
    ends.push(end_pos);

    // 4. Extract flags
    flags.push(record.flags().bits() as u32); // Convert u16 to u32

    // 5. Extract sequence if requested
    if let Some(ref mut seq_vec) = sequences {
        // Get the sequence object and keep it alive
        let sequence_obj = record.sequence();
        let seq_len = sequence_obj.len();
        let sequence_bytes = sequence_obj.as_ref();
        let mut sequence = String::with_capacity(seq_len);
        
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
            sequence.push(decode_base(base_encoded));
        }
        
        seq_vec.push(if sequence.is_empty() { None } else { Some(sequence) });
    }

    // 6. Extract quality scores if requested  
    if let Some(ref mut qual_vec) = quality_scores {
        let quality = record.quality_scores()
            .iter()
            .map(|q| char::from(q + b'!')) // Convert to ASCII PHRED+33
            .collect::<String>();
        qual_vec.push(if quality.is_empty() { None } else { Some(quality) });
    }

    Ok(())
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

fn create_record_batch(
    names: Vec<String>,
    chroms: Vec<Option<String>>,
    starts: Vec<Option<u32>>,
    ends: Vec<Option<u32>>,
    flags: Vec<u32>,
    sequences: Option<Vec<Option<String>>>,
    quality_scores: Option<Vec<Option<String>>>,
    include_sequence: bool,
    include_quality: bool,
) -> PyResult<RecordBatch> {
    // Create Arrow arrays from our data
    let mut arrays: Vec<Arc<dyn Array>> = vec![
        Arc::new(StringArray::from(names)),
        Arc::new(StringArray::from(chroms)),
        Arc::new(UInt32Array::from(starts)),
        Arc::new(UInt32Array::from(ends)),
        Arc::new(UInt32Array::from(flags)),
    ];

    if include_sequence {
        if let Some(seq_vec) = sequences {
            arrays.push(Arc::new(StringArray::from(seq_vec)));
        }
    }

    if include_quality {
        if let Some(qual_vec) = quality_scores {
            arrays.push(Arc::new(StringArray::from(qual_vec)));
        }
    }

    let schema = create_bam_schema(include_sequence, include_quality);
    RecordBatch::try_new(schema, arrays)
        .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Failed to create record batch: {}", e)))
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
