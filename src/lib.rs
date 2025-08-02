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
use parquet::basic::Compression;

use arrow::datatypes::{Schema, Field, DataType};

mod expressions;
mod single_fastq;
mod fracture;
mod graph_viz;
mod fracture_opt;
mod djfind;
mod bam;
mod bam_htslib;
// DEPRECATED: Hybrid approaches preserved as code relics (see PERFORMANCE_ROADMAP.md Phase 7)
// Performance: ~109k rec/sec vs 205k+ rec/sec optimized single-reader
// Issue: BGZF I/O serialization bottleneck with multiple concurrent readers
// mod bam_htslib_hybrid;
// mod bam_htslib_hybrid_optimized;
// mod bam_htslib_hybrid_minimal;
mod umi_score;

use crate::single_fastq::{fastq_to_parquet};
use crate::fracture::{fracture_fasta, fracture_sequences};
use crate::bam::{bam_to_parquet, bam_to_arrow_ipc, bam_to_arrow_ipc_parallel, bam_to_arrow_ipc_gzp_parallel};
#[cfg(feature = "htslib")]
use crate::bam::{bam_to_arrow_ipc_htslib_parallel, bam_to_arrow_ipc_htslib_multi_reader_parallel, bam_to_arrow_ipc_htslib_optimized, bam_to_arrow_ipc_htslib_mmap_parallel};
#[cfg(feature = "htslib")]
use crate::bam_htslib::{bam_to_arrow_ipc_htslib_bgzf_blocks};
// DEPRECATED: Hybrid function imports removed from package
// Files preserved as code relics for research/educational purposes
// #[cfg(feature = "htslib")]
// use crate::bam_htslib_hybrid::{bam_to_arrow_ipc_htslib_hybrid_segments};
// #[cfg(feature = "htslib")]
// use crate::bam_htslib_hybrid_optimized::{bam_to_arrow_ipc_htslib_hybrid_optimized};
// #[cfg(feature = "htslib")]
// use crate::bam_htslib_hybrid_minimal::{bam_to_arrow_ipc_htslib_hybrid_minimal_fix};


//extern crate fasten;
/*use seal::pair::{
    Alignment, AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, SmithWaterman, Step,Strategy,
};
*/
fn reverse_complement(dna: &str) -> String {
    dna.chars()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => c, // Handle non-DNA characters or add error handling
        })
        .rev()
        .collect::<String>()
}


#[pyfunction]
#[pyo3(signature = (in_fn1, in_fn2, out_fn, limit=None, do_rev_comp=None))]
fn merge_paired_fastqs(
       in_fn1: String,
       in_fn2: String,
       out_fn: String,
       limit: Option<usize>,
       do_rev_comp: Option<bool>,
       )
    -> PyResult<()> {

    let do_rev_comp = do_rev_comp.unwrap_or(false);

    let gz_file1 = File::open(in_fn1)?;
    let gz_file2 = File::open(in_fn2)?;

    let decoder1 = MultiGzDecoder::new(gz_file1);
    let decoder2 = MultiGzDecoder::new(gz_file2);

    let reader1 = io::BufReader::new(decoder1).lines().filter_map(Result::ok);
    let reader2 = io::BufReader::new(decoder2).lines().filter_map(Result::ok);

    let arrow_schema = Arc::new(Schema::new(vec![
        Field::new("read_id", DataType::Utf8, false),
        Field::new("r1_seq", DataType::Utf8, false),
        Field::new("r1_qual", DataType::Utf8, false),
        Field::new("r2_seq", DataType::Utf8, false),
        Field::new("r2_qual", DataType::Utf8, false),
    ]));

    // Initialize Parquet writer
    let file = File::create(out_fn)?;
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

    let mut read_id_buffer = Vec::new(); 
    let mut read1_seq_buffer= Vec::new();
    let mut read1_qual_buffer = Vec::new();
    let mut read2_seq_buffer= Vec::new();
    let mut read2_qual_buffer = Vec::new();

    let mut chunk_count = 0; 
                             
    for (chunk1, chunk2) in iter1.chunks(4).into_iter().zip(iter2.chunks(4).into_iter()) {
        let chunk1: Vec<_> = chunk1.collect();
        let chunk2: Vec<_> = chunk2.collect();

        let (read_id1, seq1, _plus1, qual1) = (&chunk1[0], &chunk1[1], &chunk1[2], &chunk1[3]);
        let (_read_id2, seq2, _plus2, qual2) = (&chunk2[0], &chunk2[1], &chunk2[2], &chunk2[3]);

        let read_id1 = read_id1.trim_start_matches('@').trim_end().to_string();

        let read1_seq = seq1.trim_end().to_string();
        let read1_qual = qual1.trim_end().to_string();

        let read2_seq = if do_rev_comp{
            reverse_complement(seq2.trim_end())
        }else{
            seq2.trim_end().to_string()
        };

        let read2_qual = if do_rev_comp {
            qual2.trim_end().chars().rev().collect::<String>()
        }else{
            qual2.trim_end().to_string()
        };
            
        read_id_buffer.push(read_id1);
        read1_seq_buffer.push(read1_seq);
        read1_qual_buffer.push(read1_qual);
        read2_seq_buffer.push(read2_seq);
        read2_qual_buffer.push(read2_qual);


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
                    Arc::new(StringArray::from(read2_seq_buffer.clone())),
                    Arc::new(StringArray::from(read2_qual_buffer.clone())),
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
            read2_seq_buffer.clear();
            read2_qual_buffer.clear();

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
            Arc::new(StringArray::from(read2_seq_buffer.clone())),
            Arc::new(StringArray::from(read2_qual_buffer.clone())),
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

#[pyfunction]
#[pyo3(signature = (in_fn1, in_fn2, cbc_len, umi_len, out_fn, limit=None, do_rev_comp=None))]
fn parse_paired_fastqs(
    // Parse pair of fastqs according to a given chemistry (10x in this case)
       in_fn1: String,
       in_fn2: String,
       cbc_len: usize,
       umi_len: usize,
       out_fn: String,
       limit: Option<usize>,
       do_rev_comp: Option<bool>)
    -> PyResult<()> {

    let do_rev_comp = do_rev_comp.unwrap_or(false);

    let gz_file1 = File::open(in_fn1)?;
    let gz_file2 = File::open(in_fn2)?;

    let decoder1 = MultiGzDecoder::new(gz_file1);
    let decoder2 = MultiGzDecoder::new(gz_file2);

    let reader1 = io::BufReader::new(decoder1).lines().filter_map(Result::ok);
    let reader2 = io::BufReader::new(decoder2).lines().filter_map(Result::ok);

    // Define schema
    let arrow_schema = Arc::new(Schema::new(vec![
        Field::new("read_id", DataType::Utf8, false),
        Field::new("start", DataType::Utf8, false),
        Field::new("end", DataType::Utf8, false),
        Field::new("cbc", DataType::Utf8, false),
        Field::new("umi", DataType::Utf8, false),
        Field::new("cbc_qual", DataType::Utf8, false),
        Field::new("umi_qual", DataType::Utf8, false),
        Field::new("seq", DataType::Utf8, false),
        Field::new("qual", DataType::Utf8, false),
    ]));

    // Initialize Parquet writer
    let file = File::create(out_fn)?;
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

    let mut read_id_buffer = Vec::new(); 
    let mut start_buffer = Vec::new(); 
    let mut end_buffer = Vec::new();
    let mut cbc_str_buffer = Vec::new(); 
    let mut cbc_qual_buffer = Vec::new();
    let mut umi_str_buffer = Vec::new(); 
    let mut umi_qual_buffer = Vec::new();
    let mut read2_seq_buffer= Vec::new();
    let mut read2_qual_buffer = Vec::new();

    let mut chunk_count = 0; 
                             
    for (chunk1, chunk2) in iter1.chunks(4).into_iter().zip(iter2.chunks(4).into_iter()) {
        let chunk1: Vec<_> = chunk1.collect();
        let chunk2: Vec<_> = chunk2.collect();

        let (read_id1, seq1, _plus1, qual1) = (&chunk1[0], &chunk1[1], &chunk1[2], &chunk1[3]);
        let (_read_id2, seq2, _plus2, qual2) = (&chunk2[0], &chunk2[1], &chunk2[2], &chunk2[3]);

        let read_id1 = read_id1.trim_start_matches('@').trim_end().to_string();
        let cbc_str = seq1.get(0..cbc_len).expect("invalid range of string").to_string();
        let umi_str = seq1.get(cbc_len..cbc_len+umi_len).expect("invalid range of string").to_string();

        let cbc_qual = qual1.get(0..cbc_len).expect("invalid range of string").to_string();
        let umi_qual = qual1.get(cbc_len..cbc_len+umi_len).expect("invalid range of string").to_string();

        let read2_seq = if do_rev_comp{
            reverse_complement(seq2.trim_end())
        }else{
            seq2.trim_end().to_string()
        };

        let read2_qual = if do_rev_comp {
            qual2.trim_end().chars().rev().collect::<String>()
        }else{
            qual2.trim_end().to_string()
        };
            
        read_id_buffer.push(read_id1);
        start_buffer.push(0.to_string());
        end_buffer.push(1.to_string());

        cbc_str_buffer.push(cbc_str);
        cbc_qual_buffer.push(cbc_qual);
        umi_str_buffer.push(umi_str);
        umi_qual_buffer.push(umi_qual);
        read2_seq_buffer.push(read2_seq);
        read2_qual_buffer.push(read2_qual);


        chunk_count += 1;
        
        // Check if chunk_count reached a million or it's the last iteration
        if chunk_count == 10_000_000   {


            // Create a record batch
            let record_batch = RecordBatch::try_new(
                arrow_schema.clone(),
                vec![
                    Arc::new(StringArray::from(read_id_buffer.clone())),
                    Arc::new(StringArray::from(start_buffer.clone())),
                    Arc::new(StringArray::from(end_buffer.clone())),
                    Arc::new(StringArray::from(cbc_str_buffer.clone())),
                    Arc::new(StringArray::from(umi_str_buffer.clone())),
                    Arc::new(StringArray::from(cbc_qual_buffer.clone())),
                    Arc::new(StringArray::from(umi_qual_buffer.clone())),
                    Arc::new(StringArray::from(read2_seq_buffer.clone())),
                    Arc::new(StringArray::from(read2_qual_buffer.clone())),
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
            start_buffer.clear();
            end_buffer.clear();

            cbc_str_buffer.clear();
            cbc_qual_buffer.clear();
            umi_str_buffer.clear();
            umi_qual_buffer.clear();
            read2_seq_buffer.clear();
            read2_qual_buffer.clear();

            chunk_count = 0;
        }
    }

    // Create a record batch
    let record_batch = RecordBatch::try_new(
        arrow_schema.clone(),
        vec![
            Arc::new(StringArray::from(read_id_buffer.clone())),
            Arc::new(StringArray::from(start_buffer.clone())),
            Arc::new(StringArray::from(end_buffer.clone())),
            Arc::new(StringArray::from(cbc_str_buffer.clone())),
            Arc::new(StringArray::from(umi_str_buffer.clone())),
            Arc::new(StringArray::from(cbc_qual_buffer.clone())),
            Arc::new(StringArray::from(umi_qual_buffer.clone())),
            Arc::new(StringArray::from(read2_seq_buffer.clone())),
            Arc::new(StringArray::from(read2_qual_buffer.clone())),
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
    start_buffer.clear();
    end_buffer.clear();
    cbc_str_buffer.clear();
    cbc_qual_buffer.clear();
    umi_str_buffer.clear();
    umi_qual_buffer.clear();


    writer.close().unwrap();
    Ok(())
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
fn oparse_cigar(cigar: &str) -> PyResult<Vec<(char, String, String)>>{ 
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
    //Some(result.into_iter().map(|(c, p, l)| format!("({}, {}, {})", c, p, l)).collect::<String>())

    Ok(result)
}



#[pymodule]
fn _internal(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

#[pymodule]
fn rogtk(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(oparse_cigar, m)?)?;
    m.add_function(wrap_pyfunction!(merge_paired_fastqs, m)?)?;
    m.add_function(wrap_pyfunction!(parse_paired_fastqs, m)?)?;
    m.add_function(wrap_pyfunction!(fastq_to_parquet, m)?)?;
    m.add_function(wrap_pyfunction!(fracture_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(fracture_sequences, m)?)?;
    m.add_function(wrap_pyfunction!(bam_to_parquet, m)?)?;
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc, m)?)?;
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_parallel, m)?)?;
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_gzp_parallel, m)?)?;
    #[cfg(feature = "htslib")]
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_parallel, m)?)?;
    #[cfg(feature = "htslib")]
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_multi_reader_parallel, m)?)?;
    #[cfg(feature = "htslib")]
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_optimized, m)?)?;
    #[cfg(feature = "htslib")]
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_mmap_parallel, m)?)?;
    #[cfg(feature = "htslib")]
    m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_bgzf_blocks, m)?)?;
    // DEPRECATED: Hybrid functions removed from Python API
    // Preserved as code relics - see PERFORMANCE_ROADMAP.md Phase 7 analysis
    // Use bam_to_arrow_ipc_htslib_optimized() for best performance (205k+ rec/sec)
    // #[cfg(feature = "htslib")]
    // m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_hybrid_segments, m)?)?;
    // #[cfg(feature = "htslib")]
    // m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_hybrid_optimized, m)?)?;
    // #[cfg(feature = "htslib")]
    // m.add_function(wrap_pyfunction!(bam_to_arrow_ipc_htslib_hybrid_minimal_fix, m)?)?;

    //m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
