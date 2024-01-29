use flate2::read::MultiGzDecoder;
use itertools::Itertools;
use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::fs::File;
use std::io::{self, BufRead};

use std::sync::Arc;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_to_parquet_schema;
use parquet::arrow::ArrowWriter;
use arrow::array::{Array, StringArray};
use parquet::file::properties::WriterProperties;
use parquet::basic::Compression;

use arrow::datatypes::{Schema, Field, DataType};
use parquet::file::writer::SerializedFileWriter;

//extern crate fasten;
/*use seal::pair::{
    Alignment, AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, SmithWaterman, Step,Strategy,
};
*/

#[pyfunction]
fn process_files(in_fn1: String, in_fn2: String, cbc_len: usize, umi_len: usize, limit: Option<usize>)
    -> PyResult<Vec<(String, String, String, String)>> {
    let gz_file1 = File::open(in_fn1)?;
    let gz_file2 = File::open(in_fn2)?;

    let decoder1 = MultiGzDecoder::new(gz_file1);
    let decoder2 = MultiGzDecoder::new(gz_file2);

    let reader1 = io::BufReader::new(decoder1).lines().filter_map(Result::ok);
    let reader2 = io::BufReader::new(decoder2).lines().filter_map(Result::ok);

    // Define schema
    let arrow_schema = Arc::new(Schema::new(vec![
        Field::new("cbc_str", DataType::Utf8, false),
        Field::new("umi_str", DataType::Utf8, false),
        Field::new("cbc_qual", DataType::Utf8, false),
        Field::new("umi_qual", DataType::Utf8, false),
    ]));

    let parquet_schema = arrow_to_parquet_schema(&arrow_schema).unwrap();

    // Initialize Parquet writer
    let file = File::create("output.parquet")?;
    // WriterProperties can be used to set Parquet file options
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();

    let mut writer = match ArrowWriter::try_new(file, arrow_schema.clone(), Some(props)){
    Ok(writer) => writer,
    Err(e) => {
        // Create a PyErr from the Parquet error
        let err_msg = format!("Parquet error: {}", e);
        return Err(PyErr::new::<PyRuntimeError, _>(err_msg));
    }
    };


    let iter1: Box<dyn Iterator<Item = String>> = match limit {
        Some(l) => Box::new(reader1.take(l)),
        None => Box::new(reader1),
    };

    let iter2: Box<dyn Iterator<Item = String>> = match limit {
        Some(l) => Box::new(reader2.take(l)),
        None => Box::new(reader2),
    };

    let mut cbc_str_buffer = Vec::new(); // Buffer to accumulate record batches
    let mut cbc_qual_buffer = Vec::new(); // Buffer to accumulate record batches
    let mut umi_str_buffer = Vec::new(); // Buffer to accumulate record batches
    let mut umi_qual_buffer = Vec::new(); // Buffer to accumulate record batches
                                          //
    let mut chunk_count = 0; // Counter for chunks
                             //
    let mut results = Vec::new();
    for (chunk1, chunk2) in iter1.chunks(4).into_iter().zip(iter2.chunks(4).into_iter()) {
        let chunk1: Vec<_> = chunk1.collect();
        let chunk2: Vec<_> = chunk2.collect();

        let (read_id1, seq1, _plus1, qual1) = (&chunk1[0], &chunk1[1], &chunk1[2], &chunk1[3]);
        let (read_id2, seq2, _plus2, qual2) = (&chunk2[0], &chunk2[1], &chunk2[2], &chunk2[3]);

        let cbc_str = seq1.get(0..cbc_len).expect("invalid range of string").to_string();
        let umi_str = seq1.get(cbc_len..cbc_len+umi_len).expect("invalid range of string").to_string();

        let cbc_qual = qual1.get(0..cbc_len).expect("invalid range of string").to_string();
        let umi_qual = qual1.get(cbc_len..cbc_len+umi_len).expect("invalid range of string").to_string();
        
        cbc_str_buffer.push(cbc_str);
        cbc_qual_buffer.push(cbc_qual);
        umi_str_buffer.push(umi_str);
        umi_qual_buffer.push(umi_qual);

        chunk_count += 1;
        
        // Check if chunk_count reached a million or it's the last iteration
        if chunk_count == 100_000   {


            // Create a record batch
            let record_batch = RecordBatch::try_new(
                arrow_schema.clone(),
                vec![
                    Arc::new(StringArray::from(cbc_str_buffer.clone())),
                    Arc::new(StringArray::from(umi_str_buffer.clone())),
                    Arc::new(StringArray::from(cbc_qual_buffer.clone())),
                    Arc::new(StringArray::from(umi_qual_buffer.clone())),
                ],
            ).map_err(|e| { // Manually handle ArrowError
                let err_msg = format!("Arrow error: {}", e);
                PyErr::new::<PyRuntimeError, _>(err_msg)
            })?;
                // Write batch to Parquet
                
            writer.write(&record_batch).map_err(|e| {
                let err_msg = format!("Parquet write error: {}", e);
                PyErr::new::<PyRuntimeError, _>(err_msg)
            })?;

            // Reset buffer and counter
            cbc_str_buffer.clear();
            cbc_qual_buffer.clear();
            umi_str_buffer.clear();
            umi_qual_buffer.clear();

            chunk_count = 0;
        }
    }

    // Create a record batch
    let record_batch = RecordBatch::try_new(
        arrow_schema.clone(),
        vec![
            Arc::new(StringArray::from(cbc_str_buffer.clone())),
            Arc::new(StringArray::from(umi_str_buffer.clone())),
            Arc::new(StringArray::from(cbc_qual_buffer.clone())),
            Arc::new(StringArray::from(umi_qual_buffer.clone())),
        ],
    ).map_err(|e| { // Manually handle ArrowError
        let err_msg = format!("Arrow error: {}", e);
        PyErr::new::<PyRuntimeError, _>(err_msg)
    })?;
        // Write batch to Parquet
        
    writer.write(&record_batch).map_err(|e| {
        let err_msg = format!("Parquet write error: {}", e);
        PyErr::new::<PyRuntimeError, _>(err_msg)
    })?;

    // Reset buffer and counter
    cbc_str_buffer.clear();
    cbc_qual_buffer.clear();
    umi_str_buffer.clear();
    umi_qual_buffer.clear();

    chunk_count = 0;

    writer.close().unwrap();
    Ok(results)
}


/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
fn porsa_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
fn parse_cigar(cigar: &str) -> PyResult<Vec<(char, String, String)>>{ 
    let mut result = Vec::new();
    let mut num_buf = String::new();
    let mut ref_pos = 0;

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buf.push(c);
        } else {
            let len = num_buf.parse::<usize>().unwrap();
            match c {
                'D' => {
                    // Decompose deletions into individual positions
                    for del_pos in ref_pos..(ref_pos + len) {
                        result.push(('D', del_pos.to_string(), 1.to_string()));
                    }
                    ref_pos += len;
                }
                'I' => {
                    // Add insertions as a single triplet
                    result.push(('I', ref_pos.to_string(), len.to_string()));
                }
                _ => {
                    // For other operations, just move the reference position
                    ref_pos += len;
                }
            }
            num_buf.clear();
        }
    }

    Ok(result)
}


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]

fn rogtk(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(porsa_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(parse_cigar, m)?)?;
    m.add_function(wrap_pyfunction!(process_files, m)?)?;
    Ok(())
}

