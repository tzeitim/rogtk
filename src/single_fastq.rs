use flate2::read::MultiGzDecoder;
use itertools::Itertools;
use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
//use pyo3::exceptions::PyValueError;
use std::fs::File;
use std::io::{self, BufRead};

use std::sync::Arc;
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use arrow::array::StringArray;
use parquet::file::properties::WriterProperties;
use parquet::basic::{Compression, ZstdLevel};

use arrow::datatypes::{Schema, Field, DataType};

// zstd-3: ~half the size of snappy at essentially snappy's write speed;
// matches the Arrow/DuckDB default. ponytail: fixed to fastq_to_parquet,
// apply the same default to bam_to_parquet/bams_to_parquet if that matters too.
fn parse_compression(name: &str, level: i32) -> PyResult<Compression> {
    match name.to_lowercase().as_str() {
        "zstd" => ZstdLevel::try_new(level)
            .map(Compression::ZSTD)
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(format!("Invalid zstd level {}: {}", level, e))),
        "snappy" => Ok(Compression::SNAPPY),
        "gzip" => Ok(Compression::GZIP(Default::default())),
        "lz4" => Ok(Compression::LZ4),
        "uncompressed" | "none" => Ok(Compression::UNCOMPRESSED),
        other => Err(PyErr::new::<PyRuntimeError, _>(format!(
            "Unknown compression '{}'. Expected one of: zstd, snappy, gzip, lz4, uncompressed", other
        ))),
    }
}

#[pyfunction]
#[pyo3(signature = (in_fn1, out_fn, limit=None, compression="zstd", compression_level=3))]
pub fn fastq_to_parquet(
       in_fn1: String,
       out_fn: String,
       limit: Option<usize>,
       compression: &str,
       compression_level: i32,
       )
    -> PyResult<()> {

    let gz_file1 = File::open(in_fn1)?;

    let decoder1 = MultiGzDecoder::new(gz_file1);

    let reader1 = io::BufReader::new(decoder1).lines().filter_map(Result::ok);

    let arrow_schema = Arc::new(Schema::new(vec![
        Field::new("read_id", DataType::Utf8, false),
        Field::new("r1_seq", DataType::Utf8, false),
        Field::new("r1_qual", DataType::Utf8, false),
    ]));

    // Initialize Parquet writer
    let file = File::create(out_fn)?;
    // WriterProperties can be used to set Parquet file options
    let props = WriterProperties::builder()
        .set_compression(parse_compression(compression, compression_level)?)
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

    let mut read_id_buffer = Vec::new(); 
    let mut read1_seq_buffer= Vec::new();
    let mut read1_qual_buffer = Vec::new();

    let mut chunk_count = 0; 
                             
    for chunk1 in iter1.chunks(4).into_iter(){
        let chunk1: Vec<_> = chunk1.collect();

        let (read_id1, seq1, _plus1, qual1) = (&chunk1[0], &chunk1[1], &chunk1[2], &chunk1[3]);

        let read_id1 = read_id1.trim_start_matches('@').trim_end().to_string();

        let read1_seq = seq1.trim_end().to_string();
        let read1_qual = qual1.trim_end().to_string();

        read_id_buffer.push(read_id1);
        read1_seq_buffer.push(read1_seq);
        read1_qual_buffer.push(read1_qual);

        chunk_count += 1;
        
        // against my expectations a relatively small buffer is faster than a larger (10M) one
        if chunk_count == 10_000{


            // Create a record batch
            let record_batch = RecordBatch::try_new(
                arrow_schema.clone(),
                vec![
                    Arc::new(StringArray::from(read_id_buffer.clone())),
                    Arc::new(StringArray::from(read1_seq_buffer.clone())),
                    Arc::new(StringArray::from(read1_qual_buffer.clone())),
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
            read_id_buffer.clear();
            read1_seq_buffer.clear();
            read1_qual_buffer.clear();

            chunk_count = 0;
        }
    }

    // Create last record batch
    let record_batch = RecordBatch::try_new(
        arrow_schema.clone(),
        vec![
            Arc::new(StringArray::from(read_id_buffer.clone())),
            Arc::new(StringArray::from(read1_seq_buffer.clone())),
            Arc::new(StringArray::from(read1_qual_buffer.clone())),
        ],
    ).map_err(|e| { 
        let err_msg = format!("Arrow error: {}", e);
        PyErr::new::<PyRuntimeError, _>(err_msg)
    })?;
    // Write last batch to Parquet
        
    writer.write(&record_batch).map_err(|e| {
        let err_msg = format!("Parquet write error: {}", e);
        PyErr::new::<PyRuntimeError, _>(err_msg)
    })?;

    writer.close().unwrap();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_compression_covers_known_and_unknown_names() {
        assert!(matches!(parse_compression("zstd", 3), Ok(Compression::ZSTD(_))));
        assert!(matches!(parse_compression("snappy", 3), Ok(Compression::SNAPPY)));
        assert!(matches!(parse_compression("uncompressed", 3), Ok(Compression::UNCOMPRESSED)));
        assert!(parse_compression("zstd", 99).is_err(), "zstd level out of range should error");
        assert!(parse_compression("bogus", 3).is_err(), "unknown codec name should error");
    }
}
