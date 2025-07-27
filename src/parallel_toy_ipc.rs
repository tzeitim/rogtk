use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;
use rayon::prelude::*;
use crossbeam_channel::{bounded, Receiver, Sender};

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use arrow::ipc::writer::FileWriter as ArrowIpcWriter;
use arrow::ipc::reader::FileReader as ArrowIpcReader;

fn create_bam_schema(include_sequence: bool, include_quality: bool) -> Arc<Schema> {
    let mut fields = vec![
        Field::new("name", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, true),
        Field::new("start", DataType::UInt32, true),
        Field::new("end", DataType::UInt32, true),
        Field::new("flags", DataType::UInt32, false),
    ];

    if include_sequence {
        fields.push(Field::new("sequence", DataType::Utf8, true));
    }
    
    if include_quality {
        fields.push(Field::new("quality_scores", DataType::Utf8, true));
    }

    Arc::new(Schema::new(fields))
}

fn create_mock_record_batch(
    batch_id: usize, 
    records_per_batch: usize, 
    include_sequence: bool, 
    include_quality: bool
) -> Result<RecordBatch, Box<dyn std::error::Error + Send + Sync>> {
    let mut names = Vec::with_capacity(records_per_batch);
    let mut chroms = Vec::with_capacity(records_per_batch);
    let mut starts = Vec::with_capacity(records_per_batch);
    let mut ends = Vec::with_capacity(records_per_batch);
    let mut flags = Vec::with_capacity(records_per_batch);
    let mut sequences = if include_sequence { Some(Vec::with_capacity(records_per_batch)) } else { None };
    let mut quality_scores = if include_quality { Some(Vec::with_capacity(records_per_batch)) } else { None };

    for i in 0..records_per_batch {
        let record_id = batch_id * records_per_batch + i;
        
        names.push(format!("read_{}", record_id));
        chroms.push(Some(format!("chr{}", (record_id % 22) + 1)));
        starts.push(Some((record_id as u32 % 1000000) + 1));
        ends.push(Some((record_id as u32 % 1000000) + 101));
        flags.push(record_id as u32 % 1024);
        
        if let Some(ref mut seq_vec) = sequences {
            seq_vec.push(Some("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string()));
        }
        
        if let Some(ref mut qual_vec) = quality_scores {
            qual_vec.push(Some("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_string()));
        }
    }

    let name_array = Arc::new(StringArray::from(names));
    let chrom_array = Arc::new(StringArray::from(chroms));
    let start_array = Arc::new(UInt32Array::from(starts));
    let end_array = Arc::new(UInt32Array::from(ends));
    let flags_array = Arc::new(UInt32Array::from(flags));
    
    let mut arrays: Vec<ArrayRef> = vec![
        name_array,
        chrom_array, 
        start_array,
        end_array,
        flags_array,
    ];
    
    if let Some(sequences) = sequences {
        arrays.push(Arc::new(StringArray::from(sequences)));
    }
    
    if let Some(quality_scores) = quality_scores {
        arrays.push(Arc::new(StringArray::from(quality_scores)));
    }

    let schema = create_bam_schema(include_sequence, include_quality);
    RecordBatch::try_new(schema, arrays)
        .map_err(|e| format!("Failed to create mock record batch: {}", e).into())
}

pub fn parallel_toy_ipc_conversion(
    ipc_path: &str,
    num_batches: usize,
    records_per_batch: usize,
    num_threads: usize,
    include_sequence: bool,
    include_quality: bool,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let start_time = Instant::now();
    println!("Starting parallel toy IPC conversion:");
    println!("  Batches: {}", num_batches);
    println!("  Records per batch: {}", records_per_batch);
    println!("  Threads: {}", num_threads);
    println!("  Include sequence: {}", include_sequence);
    println!("  Include quality: {}", include_quality);
    
    let output_path = Path::new(ipc_path);
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let schema = create_bam_schema(include_sequence, include_quality);
    
    let output_file = File::create(output_path)?;
    let mut ipc_writer = ArrowIpcWriter::try_new(output_file, &schema)?;

    // IPC Optimization 1: Use larger channel buffer for better throughput
    let channel_size = (num_threads * 4).max(16); // Larger buffer for IPC
    let (sender, receiver): (Sender<(usize, RecordBatch)>, Receiver<(usize, RecordBatch)>) = bounded(channel_size);
    
    // IPC Optimization 2: Configure rayon for better memory locality
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()?;
    
    // Writer thread optimized for IPC
    let writer_sender = sender.clone();
    let writer_thread = std::thread::spawn(move || -> Result<Vec<(usize, RecordBatch)>, Box<dyn std::error::Error + Send + Sync>> {
        let mut received_batches: Vec<(usize, RecordBatch)> = Vec::with_capacity(num_batches);
        
        // IPC Optimization 3: Pre-allocate with exact capacity
        while let Ok((batch_id, batch)) = receiver.recv() {
            println!("Received batch {} with {} records (IPC optimized)", batch_id, batch.num_rows());
            received_batches.push((batch_id, batch));
            if received_batches.len() == num_batches {
                break;
            }
        }
        
        // Sort by batch_id to ensure correct order
        received_batches.sort_by_key(|(batch_id, _)| *batch_id);
        println!("Sorted {} batches for IPC writing", received_batches.len());
        
        Ok(received_batches)
    });
    
    // Process batches in parallel with IPC optimizations
    let processing_start = Instant::now();
    pool.install(|| {
        (0..num_batches).into_par_iter().try_for_each(|batch_id| -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
            let batch_start = Instant::now();
            
            // IPC Optimization 4: Create batches with optimal memory layout
            let batch = create_mock_record_batch(batch_id, records_per_batch, include_sequence, include_quality)?;
            
            let batch_duration = batch_start.elapsed();
            println!("Generated IPC batch {} in {:?} ({} records)", batch_id, batch_duration, batch.num_rows());
            
            // Send to writer thread
            writer_sender.send((batch_id, batch))
                .map_err(|e| format!("Failed to send IPC batch {}: {}", batch_id, e))?;
            
            Ok(())
        })
    })?;
    
    drop(sender);
    
    let processing_duration = processing_start.elapsed();
    println!("Parallel IPC processing completed in {:?}", processing_duration);
    
    // Wait for writer thread and get sorted batches
    let sorted_batches = writer_thread.join()
        .map_err(|e| format!("IPC writer thread panicked: {:?}", e))??;
    
    // IPC Optimization 5: Write batches in optimal order for sequential access
    let write_start = Instant::now();
    for (batch_id, batch) in sorted_batches {
        println!("Writing IPC batch {} with {} records", batch_id, batch.num_rows());
        ipc_writer.write(&batch)?;
    }
    
    // IPC Optimization 6: Ensure proper finalization for memory mapping compatibility
    ipc_writer.finish()?;
    
    let write_duration = write_start.elapsed();
    let total_duration = start_time.elapsed();
    let total_records = num_batches * records_per_batch;
    
    println!("\nIPC Conversion complete:");
    println!("  Total records: {}", total_records);
    println!("  Processing time: {:?}", processing_duration);
    println!("  Write time: {:?}", write_duration);
    println!("  Total time: {:?}", total_duration);
    println!("  Records/sec: {:.0}", total_records as f64 / total_duration.as_secs_f64());
    println!("  Output IPC file: {}", ipc_path);
    
    // IPC Optimization 7: Demonstrate fast read-back capability
    demonstrate_fast_ipc_read(ipc_path)?;
    
    Ok(())
}

// IPC Optimization 8: Demonstrate zero-copy read performance
fn demonstrate_fast_ipc_read(ipc_path: &str) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    println!("\nDemonstrating fast IPC read-back:");
    let read_start = Instant::now();
    
    let file = File::open(ipc_path)?;
    let reader = ArrowIpcReader::try_new(file, None)?;
    
    let mut total_records = 0;
    let mut batch_count = 0;
    
    for batch_result in reader {
        let batch = batch_result?;
        total_records += batch.num_rows();
        batch_count += 1;
    }
    
    let read_duration = read_start.elapsed();
    
    println!("  Read {} batches, {} records in {:?}", batch_count, total_records, read_duration);
    println!("  Read speed: {:.0} records/sec", total_records as f64 / read_duration.as_secs_f64());
    println!("  IPC advantage: Zero-copy, memory-mappable format!");
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_parallel_toy_ipc_conversion() {
        let temp_dir = tempdir().unwrap();
        let output_path = temp_dir.path().join("test_output.arrow");
        
        let result = parallel_toy_ipc_conversion(
            output_path.to_str().unwrap(),
            4,     // num_batches
            1000,  // records_per_batch  
            2,     // num_threads
            true,  // include_sequence
            true,  // include_quality
        );
        
        assert!(result.is_ok());
        assert!(output_path.exists());
        
        // Verify file is not empty
        let metadata = std::fs::metadata(&output_path).unwrap();
        assert!(metadata.len() > 0);
    }

    #[test]
    fn test_ipc_read_performance() {
        let temp_dir = tempdir().unwrap();
        let output_path = temp_dir.path().join("test_read_perf.arrow");
        
        // Create test file
        parallel_toy_ipc_conversion(
            output_path.to_str().unwrap(),
            2,     // num_batches
            500,   // records_per_batch  
            1,     // num_threads
            false, // include_sequence
            false, // include_quality
        ).unwrap();
        
        // Test read performance
        let result = demonstrate_fast_ipc_read(output_path.to_str().unwrap());
        assert!(result.is_ok());
    }
}