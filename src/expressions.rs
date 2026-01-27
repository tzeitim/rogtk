use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
use serde::Deserialize;

use crate::fracture::assemble_sequences;
use crate::djfind::AssemblyMethod;
use crate::umi_score::calculate_umi_complexity;

use log::debug;
use std::collections::HashMap;

// ============================================================================
// CIGAR Insertion Extraction Functions
// ============================================================================

/// Represents an insertion extracted from CIGAR alignment
#[derive(Debug, Clone)]
struct Insertion {
    ref_pos: usize,
    length: usize,
    sequence: String,
}

/// Parse CIGAR string and extract insertions with their actual sequences from the query.
///
/// Returns a HashMap mapping reference position -> inserted sequence
fn extract_insertions_from_cigar(seq: &str, cigar: &str) -> HashMap<usize, String> {
    let mut insertions = HashMap::new();
    let mut num_buf = String::new();
    let mut seq_pos: usize = 0;  // position in query sequence
    let mut ref_pos: usize = 0;  // position in reference

    let seq_bytes = seq.as_bytes();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num_buf.push(c);
        } else {
            if let Ok(len) = num_buf.parse::<usize>() {
                match c {
                    'M' | '=' | 'X' => {
                        // Match/mismatch: advances both
                        seq_pos += len;
                        ref_pos += len;
                    }
                    'I' => {
                        // Insertion: extract the actual inserted bases
                        if seq_pos + len <= seq_bytes.len() {
                            let inserted = String::from_utf8_lossy(
                                &seq_bytes[seq_pos..seq_pos + len]
                            ).to_string();
                            insertions.insert(ref_pos, inserted);
                        }
                        seq_pos += len;
                        // ref_pos doesn't advance for insertions
                    }
                    'D' | 'N' => {
                        // Deletion/skip: only advances reference
                        ref_pos += len;
                    }
                    'S' => {
                        // Soft clip: only advances query
                        seq_pos += len;
                    }
                    'H' | 'P' => {
                        // Hard clip/padding: doesn't advance either
                    }
                    _ => {}
                }
            }
            num_buf.clear();
        }
    }

    insertions
}

/// Enrich an allele string by replacing [pos:NI] with [pos:NI:ACTG]
///
/// Input: "TAGTCATTAC[78:5I]ACTTAGACAGGTG"
/// Output: "TAGTCATTAC[78:5I:GCTAG]ACTTAGACAGGTG"
fn enrich_allele_with_insertions(
    allele: &str,
    insertions: &HashMap<usize, String>,
    output: &mut String
) {
    let mut chars = allele.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '[' {
            // Start of CIGAR notation - parse it
            let mut bracket_content = String::new();
            let mut found_close = false;

            for inner in chars.by_ref() {
                if inner == ']' {
                    found_close = true;
                    break;
                }
                bracket_content.push(inner);
            }

            if found_close {
                // Parse content like "78:5I" or "20:432D" or "None"
                if bracket_content == "None" {
                    output.push('[');
                    output.push_str(&bracket_content);
                    output.push(']');
                } else if let Some((pos_str, rest)) = bracket_content.split_once(':') {
                    if let Ok(pos) = pos_str.parse::<usize>() {
                        // Check if it's an insertion (ends with I)
                        if rest.ends_with('I') {
                            // It's an insertion - try to enrich with actual sequence
                            // Allele notation uses 1-based positions, CIGAR uses 0-based
                            // Try both pos-1 (convert to 0-based) and pos (in case already 0-based)
                            let ins_seq = if pos > 0 {
                                insertions.get(&(pos - 1)).or_else(|| insertions.get(&pos))
                            } else {
                                insertions.get(&pos)
                            };
                            if let Some(seq) = ins_seq {
                                output.push('[');
                                output.push_str(&bracket_content);
                                output.push(':');
                                output.push_str(seq);
                                output.push(']');
                            } else {
                                // No insertion found at this position, keep original
                                output.push('[');
                                output.push_str(&bracket_content);
                                output.push(']');
                            }
                        } else {
                            // Not an insertion (deletion, etc.) - keep as-is
                            output.push('[');
                            output.push_str(&bracket_content);
                            output.push(']');
                        }
                    } else {
                        // Couldn't parse position - keep original
                        output.push('[');
                        output.push_str(&bracket_content);
                        output.push(']');
                    }
                } else {
                    // Invalid format - keep original
                    output.push('[');
                    output.push_str(&bracket_content);
                    output.push(']');
                }
            } else {
                // No closing bracket found - just output what we have
                output.push('[');
                output.push_str(&bracket_content);
            }
        } else {
            output.push(c);
        }
    }
}

/// Polars expression: Enrich allele strings with actual insertion sequences
///
/// Takes 3 inputs:
/// - inputs[0]: allele string column (e.g., "TAGTCATTAC[78:5I]ACTTAGACAGGTG")
/// - inputs[1]: sequence column (the full aligned sequence)
/// - inputs[2]: CIGAR string column
///
/// Returns: enriched allele string with insertion sequences
#[polars_expr(output_type=String)]
fn enrich_allele_insertions_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let allele_ca: &StringChunked = inputs[0].str()?;
    let seq_ca: &StringChunked = inputs[1].str()?;
    let cigar_ca: &StringChunked = inputs[2].str()?;

    // Process each row
    let results: StringChunked = allele_ca
        .into_iter()
        .zip(seq_ca.into_iter())
        .zip(cigar_ca.into_iter())
        .map(|((allele_opt, seq_opt), cigar_opt)| {
            match (allele_opt, seq_opt, cigar_opt) {
                (Some(allele), Some(seq), Some(cigar)) => {
                    let insertions = extract_insertions_from_cigar(seq, cigar);
                    let mut output = String::with_capacity(allele.len() + 50);
                    enrich_allele_with_insertions(allele, &insertions, &mut output);
                    Some(output)
                }
                (Some(allele), _, _) => Some(allele.to_string()), // Return original if missing seq/cigar
                _ => None,
            }
        })
        .collect();

    Ok(results.into_series())
}

/// Polars expression: Extract all insertions from CIGAR as a string map
///
/// Takes 2 inputs:
/// - inputs[0]: sequence column
/// - inputs[1]: CIGAR string column
///
/// Returns: String in format "pos1:seq1|pos2:seq2|..."
#[polars_expr(output_type=String)]
fn extract_cigar_insertions_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let seq_ca: &StringChunked = inputs[0].str()?;
    let cigar_ca: &StringChunked = inputs[1].str()?;

    let results: StringChunked = seq_ca
        .into_iter()
        .zip(cigar_ca.into_iter())
        .map(|(seq_opt, cigar_opt)| {
            match (seq_opt, cigar_opt) {
                (Some(seq), Some(cigar)) => {
                    let insertions = extract_insertions_from_cigar(seq, cigar);
                    if insertions.is_empty() {
                        Some(String::new())
                    } else {
                        let mut pairs: Vec<_> = insertions.into_iter().collect();
                        pairs.sort_by_key(|(pos, _)| *pos);
                        let result = pairs
                            .into_iter()
                            .map(|(pos, seq)| format!("{}:{}", pos, seq))
                            .collect::<Vec<_>>()
                            .join("|");
                        Some(result)
                    }
                }
                _ => None,
            }
        })
        .collect();

    Ok(results.into_series())
}

// ============================================================================
// Original CIGAR parsing (indel positions only, no sequences)
// ============================================================================

fn parse_cigar_str(cigar: &str, output: &mut String, block_dels: bool) {
    let mut num_buf = String::new();
    let mut ref_pos = 0;

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buf.push(c);
        } else {
            if let Ok(len) = num_buf.parse::<usize>() {
                match c {
                    'D' => {
                        if block_dels {
                            output.push_str(&format!("D,{},{}|", ref_pos, len));
                        } else{
                            for del_pos in ref_pos..(ref_pos + len) {
                                output.push_str(&format!("D,{},1|", del_pos));
                            }
                        }
                        ref_pos += len;
                    }
                    'I' => {
                        output.push_str(&format!("I,{},{}|", ref_pos, len));
                    }
                    _ => {
                        ref_pos += len;
                    }
                }
            }
            num_buf.clear();
        }
    }

    if output.ends_with('|') {
        output.pop();
    }
}

#[derive(Deserialize)]
struct ParseCigarKwargs {
    block_dels: bool,
}

#[derive(Deserialize)]
struct Phred2NmKwargs {
    base: u8,
}

#[polars_expr(output_type=String)]
fn parse_cigar_series(inputs: &[Series], kwargs: ParseCigarKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    //let out: StringChunked = ca.apply_into_string_amortized(parse_cigar_str);
    let out: StringChunked = ca.apply_into_string_amortized(|value, output|{
            parse_cigar_str(value, output, kwargs.block_dels);
    });
    Ok(out.into_series())
}

fn __phred_to_numeric_str(phred_char: &str, output: &mut String, base: u8){
    for phred_char in phred_char.chars() {
    // PHRED+33 format (common in FASTQ files)
        output.push_str(&format!("{}|",phred_char as u8 - base));
    }
    
    if output.ends_with('|') {
        output.pop();
    }
}
//fn phred_to_numeric(phred_char: &str, output: &mut u8, base: u8){
//    for phred_char in phred_char.chars() {
//    // PHRED+33 format (common in FASTQ files)
//        output.push(phred_char as u8 - base);
//    }
//}

fn list_u8_dtype(input_fields: &[Field]) -> PolarsResult<Field> {
    let field = Field::new(input_fields[0].name.clone(), DataType::List(Box::new(DataType::UInt8)));
    Ok(field.clone())
}

fn _list_i32_dtype(input_fields: &[Field]) -> PolarsResult<Field> {
    let field = Field::new(input_fields[0].name.clone(), DataType::List(Box::new(DataType::Int32)));
    Ok(field.clone())
}

fn _list_str_dtype(input_fields: &[Field]) -> PolarsResult<Field> {
    let field = Field::new(input_fields[0].name.clone(), DataType::List(Box::new(DataType::String)));
    Ok(field.clone())
}

#[polars_expr(output_type_func=list_u8_dtype)]
fn panic_phred_to_numeric_series(inputs: &[Series], kwargs: Phred2NmKwargs) -> PolarsResult<Series> {
    //let values = inputs[0].list()?;
    let ca: &StringChunked = inputs[0].str()?;

    let name = PlSmallStr::from_str("numeric_phred");
    let mut output: ListPrimitiveChunkedBuilder<UInt8Type> = ListPrimitiveChunkedBuilder::new(
        name,
        ca.len(), 
        ca.len(), 
        DataType::UInt8,
    );

    ca.for_each(|phred_str_opt: Option<&str>| {
        if let Some(phred_str) = phred_str_opt {
            let mut byte_array: Vec<u8> = Vec::with_capacity(phred_str.len());

            for phred_char in phred_str.chars() {
                // PHRED+33 format (common in FASTQ files)
                let phred_num: u8 = phred_char as u8 - kwargs.base;
                byte_array.push(phred_num);
            }
            output.append_slice(&byte_array);
            dbg!(&byte_array); // Debug print the byte array
            byte_array.clear();
        }
    });
    let result = output.finish();
    println!("Hello, rusty world! 113");
    println!("DataType of result: {:?}", result.dtype());

    Ok(result.into_series())
}
//#[polars_expr(output_type_func=list_i32_dtype)]
//fn nn(inputs: &[Series]) -> PolarsResult<Series> {
//    let s = &inputs[0];
//    let ca = s.str()?;
//    //polars_ensure!(
//    //    ca.dtype() == &DataType::List(Box::new(DataType::String)),
//    //    ComputeError: "Expected `List(String)`, got: {}", ca.dtype()
//    //);
//
//    let out: ListChunked = ca.apply(
//        |opt_s: Option<&str>| {
//            // Modified core logic
//            let name = "my_chunked_array".into();
//            let vec = vec![1, 2, 3, 4, 5];
//            let int_arr = ChunkedArray::<Int32Type>::from_vec(name, vec);
//            int_arr
//        }
//    ).to_list();
//
//    eprintln!("oilo");
//    out.into_series();
//    let series: Series = out.into_series();
//        //println!("{}", out);
//    Ok(series)
//
//}
#[polars_expr(output_type_func=list_u8_dtype)]
fn phred_to_numeric_series(inputs: &[Series], kwargs: Phred2NmKwargs) -> PolarsResult<Series> {
    //let values = inputs[0].list()?;
    let ca: &StringChunked = inputs[0].str()?;

    let name = PlSmallStr::from_str("numeric_phred");
    let mut output: ListPrimitiveChunkedBuilder<UInt8Type> = ListPrimitiveChunkedBuilder::new(
        name,
        ca.len(), 
        ca.len(), 
        DataType::UInt8,
    );

    ca.for_each(|phred_str_opt: Option<&str>| {
        if let Some(phred_str) = phred_str_opt {
            let mut byte_array: Vec<u8> = Vec::with_capacity(phred_str.len());

            for phred_char in phred_str.chars() {
                // PHRED+33 format (common in FASTQ files)
                let phred_num: u8 = phred_char as u8 - kwargs.base;
                byte_array.push(phred_num);
            }
            output.append_slice(&byte_array);
            dbg!(&byte_array); // Debug print the byte array
            byte_array.clear();
        }
    });
    let result = output.finish();
    println!("Hello, rusty world! 113");
    println!("DataType of result: {:?}", result.dtype());

    Ok(result.into_series())
}

#[polars_expr(output_type=String)]
fn phred_to_numeric_series_str(inputs: &[Series], kwargs: Phred2NmKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    //let out: StringChunked = ca.apply_into_string_amortized(parse_cigar_str);
    let out: StringChunked = ca.apply_into_string_amortized(|value, output|{
            split_string(value, output, kwargs.base);
    });
    Ok(out.into_series())
}

//#[polars_expr(output_type=list_str_dtype)]
//fn phred_to_numeric_series_str(inputs: &[Series], kwargs: Phred2NmKwargs) -> PolarsResult<Series> {
//    let ca: &StringChunked = inputs[0].str()?;
//    //let out: StringChunked = ca.apply_into_string_amortized(parse_cigar_str);
//    let out: StringChunked = ca.apply_into_string_amortized(|value, output|{
//            split_string(value, output, kwargs.base);
//    });
//    
//    let result = out.into_series().str()?.split(lit("|")).list().head(lit(1));
//
//    Ok(result)
//}

fn split_string(phred_str: &str, output: &mut String, base: u8) {

    for phred_char in phred_str.chars() {
        //output.push_str(&format!("{}|", phred_char));
        output.push_str(&format!("{}|",phred_char as u8 - base));
    }

    if output.ends_with('|') {
        output.pop();
    }
}


// fracture
#[derive(Deserialize, Debug)]
#[allow(dead_code)]
struct AssemblyKwargs {
    k: usize,
    min_coverage: usize,
    method: String,
    // optional parameters
    start_anchor: Option<String>,
    end_anchor: Option<String>,
    min_length: Option<usize>,
    export_graphs: Option<bool>,
    only_largest: Option<bool>,
    auto_k: Option<bool>,
    prefix: Option<String>,
    // optimization related
    max_iterations: Option<bool>,
    explore_k: Option<bool>,
    prioritize_lenth: Option<bool>,
}

// Default string output type for the expression
fn output_string_type(input_fields: &[Field]) -> PolarsResult<Field> {
    let field = Field::new(input_fields[0].name().clone(), DataType::String);
    Ok(field)
}

#[polars_expr(output_type_func=output_string_type)]
fn assemble_sequences_expr(inputs: &[Series], kwargs: AssemblyKwargs) -> PolarsResult<Series> {
    debug!("Received kwargs: {:?}", kwargs);
    

    let method = match kwargs.method.as_str() {
        "compression" => {
            if kwargs.start_anchor.is_some() || kwargs.end_anchor.is_some() {
                return Err(PolarsError::ComputeError(
                        "Anchor sequences should not be provided for compression method".into()
                ));
            }
            AssemblyMethod::Compression
        },
        "shortest_path" => {
            match (&kwargs.start_anchor, &kwargs.end_anchor) {
                (Some(start), Some(end)) => AssemblyMethod::ShortestPath {
                    start_anchor: start.clone(),
                    end_anchor: end.clone(),
                },
                _ => return Err(PolarsError::ComputeError(
                        "Both start_anchor and end_anchor are required for shortest_path method".into()
                )),
            }
        },
        "shortest_path_auto" => {
            if kwargs.start_anchor.is_some() || kwargs.end_anchor.is_some() {
                return Err(PolarsError::ComputeError(
                        "Anchor sequences should not be provided for shortest_path_auto method".into()
                ));
            }
            AssemblyMethod::ShortestPathAuto
        },
        _ => return Err(PolarsError::ComputeError(
                "Invalid assembly method. Must be 'compression', 'shortest_path', or 'shortest_path_auto'".into()
        )),
    };
    
    // Extract sequences from input series
    let ca = inputs[0].str()?;
    
    // Convert string chunk to Vec<String>
    let sequences: Vec<String> = ca.into_iter()
        .flatten()
        .map(|s| s.to_string())
        .collect();
    
    // Call assembly function with only_largest always set to true
    let contigs = assemble_sequences(
        sequences,
        kwargs.k,
        kwargs.min_coverage,
        method,
        kwargs.export_graphs,
        Some(true), // Hardcoded to true
        kwargs.min_length,
        kwargs.auto_k,
        kwargs.prefix,
    ).map_err(|e| PolarsError::ComputeError(
        format!("Assembly failed: {}", e).into()
    ))?;

    let result = contigs.join("\n");
    Ok(StringChunked::from_slice(
        PlSmallStr::from_str("assembled_sequences"),
        &[result.as_str()]
    ).into_series())
}

/// Assembly expression that takes anchor sequences from input columns instead of kwargs.
/// This enables per-group dynamic anchors in group_by operations.
///
/// inputs[0] = sequences (String column)
/// inputs[1] = start_anchor (String column - first value used)
/// inputs[2] = end_anchor (String column - first value used)
#[polars_expr(output_type_func=output_string_type)]
fn assemble_sequences_with_anchors_expr(inputs: &[Series], kwargs: AssemblyKwargs) -> PolarsResult<Series> {
    debug!("Received kwargs for dynamic anchors: {:?}", kwargs);

    // Validate we have at least 3 inputs
    if inputs.len() < 3 {
        return Err(PolarsError::ComputeError(
            "assemble_sequences_with_anchors requires 3 inputs: sequences, start_anchor, end_anchor".into()
        ));
    }

    // Extract anchor strings from the first row of each anchor column
    let start_anchor_series = inputs[1].str()?;
    let end_anchor_series = inputs[2].str()?;

    let start_anchor = start_anchor_series.get(0)
        .ok_or_else(|| PolarsError::ComputeError(
            "start_anchor column is empty".into()
        ))?
        .to_string();

    let end_anchor = end_anchor_series.get(0)
        .ok_or_else(|| PolarsError::ComputeError(
            "end_anchor column is empty".into()
        ))?
        .to_string();

    debug!("Dynamic anchors - start: {}, end: {}", start_anchor, end_anchor);

    // Build assembly method - for dynamic anchors, we only support shortest_path
    let method = match kwargs.method.as_str() {
        "compression" => {
            return Err(PolarsError::ComputeError(
                "compression method is not supported with dynamic anchors; use shortest_path".into()
            ));
        },
        "shortest_path" => AssemblyMethod::ShortestPath {
            start_anchor: start_anchor.clone(),
            end_anchor: end_anchor.clone(),
        },
        "shortest_path_auto" => {
            return Err(PolarsError::ComputeError(
                "shortest_path_auto method is not supported with dynamic anchors; use shortest_path".into()
            ));
        },
        _ => return Err(PolarsError::ComputeError(
            "Invalid assembly method for dynamic anchors. Must be 'shortest_path'".into()
        )),
    };

    // Extract sequences from input series
    let ca = inputs[0].str()?;

    // Convert string chunk to Vec<String>
    let sequences: Vec<String> = ca.into_iter()
        .flatten()
        .map(|s| s.to_string())
        .collect();

    // Call assembly function
    let contigs = assemble_sequences(
        sequences,
        kwargs.k,
        kwargs.min_coverage,
        method,
        kwargs.export_graphs,
        Some(true), // Hardcoded to true
        kwargs.min_length,
        kwargs.auto_k,
        kwargs.prefix,
    ).map_err(|e| PolarsError::ComputeError(
        format!("Assembly failed: {}", e).into()
    ))?;

    let result = contigs.join("\n");
    Ok(StringChunked::from_slice(
        PlSmallStr::from_str("assembled_sequences"),
        &[result.as_str()]
    ).into_series())
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct SweepParams {
    k_start: usize,
    k_end: usize,
    k_step: usize,
    cov_start: usize,
    cov_end: usize, 
    cov_step: usize,
    method: String,
    start_anchor: Option<String>,
    end_anchor: Option<String>,
    min_length: Option<usize>,
    export_graphs: Option<bool>,
    prefix: Option<String>,
    auto_k: Option<bool>,
}

fn struct_output_type(input_fields: &[Field]) -> PolarsResult<Field> {
    let fields = vec![
        Field::new("k".into(), DataType::Int64),
        Field::new("min_coverage".into(), DataType::Int64),
        Field::new("contig_length".into(), DataType::Int64),
    ];
    let struct_type = DataType::Struct(fields);
    let field = Field::new(input_fields[0].name().clone(), struct_type);
    Ok(field)
}

#[polars_expr(output_type_func=struct_output_type)]
fn sweep_assembly_params_expr(inputs: &[Series], kwargs: SweepParams) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    
    let sequences: Vec<String> = ca.into_iter()
        .flatten()
        .map(|s| s.to_string())
        .collect();
    
    let method = AssemblyMethod::from_str(
        &kwargs.method,
        kwargs.start_anchor.clone(),
        kwargs.end_anchor.clone(),
    ).map_err(|e| PolarsError::ComputeError(
        format!("Invalid assembly method: {}", e).into()
    ))?;

    let mut k_values = Vec::new();
    let mut cov_values = Vec::new();
    let mut len_values = Vec::new();

    for k in (kwargs.k_start..=kwargs.k_end).step_by(kwargs.k_step) {
        for min_cov in (kwargs.cov_start..=kwargs.cov_end).step_by(kwargs.cov_step) {
            // Run assembly with current parameters
            //pub fn assemble_sequences(
            //    sequences: Vec<String>, 
            //    k: usize, 
            //    min_coverage: usize, 
            //    method: AssemblyMethod, 
            //    export_graphs: Option<bool>,
            //    only_largest: Option<bool>,
            //    min_length: Option<usize>,
            //    auto_k: Option<bool>,
            //    prefix: Option<String>,
            //
            match assemble_sequences(
                sequences.clone(),
                k,
                min_cov,
                // TODO: Optimize by using references to AssemblyMethod instead of cloning.
                // This would require updating function signatures across the codebase 
                // to accept &AssemblyMethod instead of AssemblyMethod.
                method.clone(),
                kwargs.export_graphs,
                Some(true),
                None,
                Some(false),
                kwargs.prefix.clone(),
            ) {
                Ok(contigs) => {
                    let length = if contigs.is_empty() { 0 } else { contigs[0].len() };
                    k_values.push(k as i64);
                    cov_values.push(min_cov as i64);
                    len_values.push(length as i64);
                }
                Err(_) => {
                    k_values.push(k as i64);
                    cov_values.push(min_cov as i64);
                    len_values.push(0i64);
                }
            }
        }
    }

    // Create series for each field
    let k_series = Series::new("k".into(), k_values);
    let cov_series = Series::new("min_coverage".into(), cov_values);
    let len_series = Series::new("contig_length".into(), len_values);

    // Create struct from series
    let fields = vec![k_series, cov_series, len_series];
    let df = DataFrame::new(fields)?;
    
    // Use the input field name for the output struct - no need for extra .into()
    Ok(df.into_struct(inputs[0].name().clone()).into())
}

#[polars_expr(output_type=String)]
fn reverse_complement_series(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    let out: StringChunked = ca.apply_into_string_amortized(|value, output| {
        reverse_complement_to_output(value, output);
    });
    Ok(out.into_series())
}

fn reverse_complement_to_output(dna: &str, output: &mut String) {
    for c in dna.chars().rev() {
        match c {
            'A' => output.push('T'),
            'T' => output.push('A'),
            'C' => output.push('G'),
            'G' => output.push('C'),
            'N' => output.push('N'),
            _ => output.push(c), 
        }
    }
}

/////////////////////////
// Helper function to generate fuzzy patterns in Rust
//
//
fn generate_fuzzy_pattern(string: &str, wildcard: &str, include_original: bool, max_length: usize) -> String {
    if string.is_empty() {
        return string.to_string();
    }
    
    let mut fuzz = Vec::new();
    
    if include_original {
        fuzz.push(string.to_string());
    }
    
    // Skip fuzzy generation for very long strings to avoid performance issues
    if string.len() <= max_length {
        for i in 0..string.len() {
            let variant: String = string.chars()
                .enumerate()
                .map(|(j, c)| if i == j { wildcard.to_string() } else { c.to_string() })
                .collect::<Vec<String>>()
                .join("");
            fuzz.push(variant);
        }
        
        // Add end-substitution pattern: replace last character with .
        if string.len() > 0 {
            let end_substitution = format!("{}.", &string[..string.len()-1]);
            fuzz.push(end_substitution);
        }
    }
    
    fuzz.join("|")
}


#[derive(Deserialize, Debug)]
struct HammingKwargs {
    target: String,
    max_distance: Option<u32>,
}

#[derive(Deserialize, Debug)]
struct FuzzyReplaceKwargs {
    pattern: String,
    replacement: String,
    literal: Option<bool>,
}

#[derive(Deserialize, Debug)]
struct FuzzyMatchKwargs {
    target: String,
    wildcard: Option<String>,
    include_original: Option<bool>,
    max_length: Option<usize>,
}

#[derive(Deserialize, Debug)]
struct FuzzyReplaceNativeKwargs {
    target: String,
    replacement: String,
    wildcard: Option<String>,
    include_original: Option<bool>,
    max_length: Option<usize>,
    replace_all: Option<bool>,
}

// Hamming distance functions
#[polars_expr(output_type=UInt32)]
fn hamming_distance_expr(inputs: &[Series], kwargs: HammingKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    let target = &kwargs.target;
    let target_chars: Vec<char> = target.chars().collect();
    
    let results: Vec<Option<u32>> = ca.iter().map(|seq_opt| {
        match seq_opt {
            Some(seq) => {
                if seq.len() != target.len() { 
                    Some(u32::MAX) // Different lengths = invalid
                } else {
                    let distance: u32 = seq.chars()
                       .zip(target_chars.iter())
                       .map(|(a, b)| if a != *b { 1u32 } else { 0u32 })
                       .sum();
                    Some(distance)
                }
            }
            None => None,
        }
    }).collect();
    
    let out = UInt32Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=Boolean)]
fn hamming_within_expr(inputs: &[Series], kwargs: HammingKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    let target = &kwargs.target;
    let max_dist = kwargs.max_distance.unwrap_or(1);
    let target_chars: Vec<char> = target.chars().collect();
    
    let results: Vec<Option<bool>> = ca.iter().map(|seq_opt| {
        match seq_opt {
            Some(seq) => {
                if seq.len() != target.len() { 
                    Some(false) 
                } else {
                    let distance: u32 = seq.chars()
                       .zip(target_chars.iter())
                       .map(|(a, b)| if a != *b { 1u32 } else { 0u32 })
                       .sum();
                    Some(distance <= max_dist)
                }
            }
            None => None,
        }
    }).collect();
    
    let out = BooleanChunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

// Fuzzy pattern functions with pre-generated patterns
#[polars_expr(output_type=String)]
fn fuzzy_replace_expr(inputs: &[Series], kwargs: FuzzyReplaceKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    let pattern = &kwargs.pattern;
    let replacement = &kwargs.replacement;
    let literal = kwargs.literal.unwrap_or(false);
    
    if literal {
        // Simple string replacement
        let out: StringChunked = ca.apply_into_string_amortized(|value, output| {
            output.push_str(&value.replace(pattern, replacement));
        });
        Ok(out.into_series())
    } else {
        // Regex replacement - compile once for efficiency
        let re = regex::Regex::new(pattern).map_err(|e| {
            PolarsError::ComputeError(format!("Invalid regex pattern: {}", e).into())
        })?;
        
        let out: StringChunked = ca.apply_into_string_amortized(|value, output| {
            let result = re.replace_all(value, replacement);
            output.push_str(&result);
        });
        Ok(out.into_series())
    }
}

#[polars_expr(output_type=Boolean)]
fn fuzzy_contains_expr(inputs: &[Series], kwargs: FuzzyReplaceKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    let pattern = &kwargs.pattern;
    let literal = kwargs.literal.unwrap_or(false);
    
    if literal {
        let results: Vec<Option<bool>> = ca.iter().map(|value_opt| {
            match value_opt {
                Some(value) => Some(value.contains(pattern)),
                None => None,
            }
        }).collect();
        
        let out = BooleanChunked::from_iter_options(ca.name().clone(), results.into_iter());
        Ok(out.into_series())
    } else {
        let re = regex::Regex::new(pattern).map_err(|e| {
            PolarsError::ComputeError(format!("Invalid regex pattern: {}", e).into())
        })?;
        
        let results: Vec<Option<bool>> = ca.iter().map(|value_opt| {
            match value_opt {
                Some(value) => Some(re.is_match(value)),
                None => None,
            }
        }).collect();
        
        let out = BooleanChunked::from_iter_options(ca.name().clone(), results.into_iter());
        Ok(out.into_series())
    }
}

#[polars_expr(output_type=Boolean)]
fn fuzzy_contains_native_expr(inputs: &[Series], kwargs: FuzzyMatchKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    let target = &kwargs.target;
    let wildcard = kwargs.wildcard.as_deref().unwrap_or(".{0,1}");
    let include_original = kwargs.include_original.unwrap_or(true);
    let max_length = kwargs.max_length.unwrap_or(100);

    // Generate fuzzy pattern in Rust
    let pattern = generate_fuzzy_pattern(target, wildcard, include_original, max_length);

    let re = regex::Regex::new(&pattern).map_err(|e| {
        PolarsError::ComputeError(format!("Invalid regex pattern: {}", e).into())
    })?;

    let results: Vec<Option<bool>> = ca.iter().map(|value_opt| {
        match value_opt {
            Some(value) => Some(re.is_match(value)),
            None => None,
        }
    }).collect();

    let out = BooleanChunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=String)]
fn fuzzy_replace_native_expr(inputs: &[Series], kwargs: FuzzyReplaceNativeKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    let target = &kwargs.target;
    let replacement = &kwargs.replacement;
    let wildcard = kwargs.wildcard.as_deref().unwrap_or(".{0,1}");
    let include_original = kwargs.include_original.unwrap_or(true);
    let max_length = kwargs.max_length.unwrap_or(100);
    let replace_all = kwargs.replace_all.unwrap_or(false);    

    let pattern = generate_fuzzy_pattern(target, wildcard, include_original, max_length);

    let re = regex::Regex::new(&pattern).map_err(|e| {
        PolarsError::ComputeError(format!("Invalid regex pattern: {}", e).into())
    })?;

    let out: StringChunked = ca.apply_into_string_amortized(|value, output| {
        let result = if replace_all {
            re.replace_all(value, replacement)
        } else {
            re.replace(value, replacement)  
        };
        output.push_str(&result);
    });

    Ok(out.into_series())
}

// UMI Complexity Scoring Functions
fn umi_complexity_struct_output_type(input_fields: &[Field]) -> PolarsResult<Field> {
    let fields = vec![
        Field::new("shannon_entropy".into(), DataType::Float64),
        Field::new("linguistic_complexity".into(), DataType::Float64),
        Field::new("homopolymer_fraction".into(), DataType::Float64),
        Field::new("dinucleotide_entropy".into(), DataType::Float64),
        Field::new("longest_homopolymer_run".into(), DataType::UInt32),
        Field::new("dust_score".into(), DataType::Float64),
        Field::new("combined_score".into(), DataType::Float64),
    ];
    let struct_type = DataType::Struct(fields);
    let field = Field::new(input_fields[0].name().clone(), struct_type);
    Ok(field)
}

#[polars_expr(output_type_func=umi_complexity_struct_output_type)]
fn umi_complexity_all_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let mut shannon_values = Vec::new();
    let mut linguistic_values = Vec::new();
    let mut homopolymer_values = Vec::new();
    let mut dinucleotide_values = Vec::new();
    let mut longest_homo_values = Vec::new();
    let mut dust_values = Vec::new();
    let mut combined_values = Vec::new();
    
    for umi_opt in ca.iter() {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                shannon_values.push(Some(scores.shannon_entropy));
                linguistic_values.push(Some(scores.linguistic_complexity));
                homopolymer_values.push(Some(scores.homopolymer_fraction));
                dinucleotide_values.push(Some(scores.dinucleotide_entropy));
                longest_homo_values.push(Some(scores.longest_homopolymer_run as u32));
                dust_values.push(Some(scores.dust_score));
                combined_values.push(Some(scores.combined_score));
            }
            None => {
                shannon_values.push(None);
                linguistic_values.push(None);
                homopolymer_values.push(None);
                dinucleotide_values.push(None);
                longest_homo_values.push(None);
                dust_values.push(None);
                combined_values.push(None);
            }
        }
    }
    
    // Create series for each field
    let shannon_series = Series::new("shannon_entropy".into(), shannon_values);
    let linguistic_series = Series::new("linguistic_complexity".into(), linguistic_values);
    let homopolymer_series = Series::new("homopolymer_fraction".into(), homopolymer_values);
    let dinucleotide_series = Series::new("dinucleotide_entropy".into(), dinucleotide_values);
    let longest_homo_series = Series::new("longest_homopolymer_run".into(), longest_homo_values);
    let dust_series = Series::new("dust_score".into(), dust_values);
    let combined_series = Series::new("combined_score".into(), combined_values);
    
    // Create struct from series
    let fields = vec![shannon_series, linguistic_series, homopolymer_series, dinucleotide_series, longest_homo_series, dust_series, combined_series];
    let df = DataFrame::new(fields)?;
    
    Ok(df.into_struct(inputs[0].name().clone()).into())
}

#[polars_expr(output_type=Float64)]
fn umi_shannon_entropy_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let results: Vec<Option<f64>> = ca.iter().map(|umi_opt| {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                Some(scores.shannon_entropy)
            }
            None => None,
        }
    }).collect();
    
    let out = Float64Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=Float64)]
fn umi_linguistic_complexity_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let results: Vec<Option<f64>> = ca.iter().map(|umi_opt| {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                Some(scores.linguistic_complexity)
            }
            None => None,
        }
    }).collect();
    
    let out = Float64Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=Float64)]
fn umi_homopolymer_fraction_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let results: Vec<Option<f64>> = ca.iter().map(|umi_opt| {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                Some(scores.homopolymer_fraction)
            }
            None => None,
        }
    }).collect();
    
    let out = Float64Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=Float64)]
fn umi_dinucleotide_entropy_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let results: Vec<Option<f64>> = ca.iter().map(|umi_opt| {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                Some(scores.dinucleotide_entropy)
            }
            None => None,
        }
    }).collect();
    
    let out = Float64Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=Float64)]
fn umi_combined_score_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let results: Vec<Option<f64>> = ca.iter().map(|umi_opt| {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                Some(scores.combined_score)
            }
            None => None,
        }
    }).collect();
    
    let out = Float64Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=UInt32)]
fn umi_longest_homopolymer_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let results: Vec<Option<u32>> = ca.iter().map(|umi_opt| {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                Some(scores.longest_homopolymer_run as u32)
            }
            None => None,
        }
    }).collect();
    
    let out = UInt32Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}

#[polars_expr(output_type=Float64)]
fn umi_dust_score_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    
    let results: Vec<Option<f64>> = ca.iter().map(|umi_opt| {
        match umi_opt {
            Some(umi) => {
                let scores = calculate_umi_complexity(umi);
                Some(scores.dust_score)
            }
            None => None,
        }
    }).collect();
    
    let out = Float64Chunked::from_iter_options(ca.name().clone(), results.into_iter());
    Ok(out.into_series())
}
