use flate2::read::MultiGzDecoder;
use itertools::Itertools;
use pyo3::prelude::*;
use std::fs::File;
use std::io::{self, BufRead};


//extern crate fasten;
/*use seal::pair::{
    Alignment, AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, SmithWaterman, Step,Strategy,
};
*/

#[pyfunction]
fn process_files(in_fn1: String, in_fn2: String, limit: Option<usize>) -> PyResult<Vec<(Vec<String>, Vec<String>)>> {
    let gz_file1 = File::open(in_fn1)?;
    let gz_file2 = File::open(in_fn2)?;

    let decoder1 = MultiGzDecoder::new(gz_file1);
    let decoder2 = MultiGzDecoder::new(gz_file2);

    let reader1 = io::BufReader::new(decoder1).lines().filter_map(Result::ok);
    let reader2 = io::BufReader::new(decoder2).lines().filter_map(Result::ok);

    let iter1: Box<dyn Iterator<Item = String>> = match limit {
        Some(l) => Box::new(reader1.take(l)),
        None => Box::new(reader1),
    };

    let iter2: Box<dyn Iterator<Item = String>> = match limit {
        Some(l) => Box::new(reader2.take(l)),
        None => Box::new(reader2),
    };

    let mut results = Vec::new();
    for (chunk1, chunk2) in iter1.chunks(4).into_iter().zip(iter2.chunks(4).into_iter()) {
        let chunk1: Vec<_> = chunk1.collect();
        let chunk2: Vec<_> = chunk2.collect();
        results.push((chunk1, chunk2));
    }

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

