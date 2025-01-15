use serde::Deserialize;
use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
use crate::fracture::assemble_sequences;
use crate::djfind::AssemblyMethod;

use log::debug;

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

    // Create AssemblyMethod based on method string
    let method = match kwargs.method.as_str() {
        "compression" => Ok(AssemblyMethod::Compression),
        "shortest_path" => AssemblyMethod::from_str(
            &kwargs.method,
            kwargs.start_anchor,
            kwargs.end_anchor,
        ),
        _ => Err("Invalid assembly method. Must be 'compression' or 'shortest_path'".to_string())
    };

    // Handle any errors from method parsing
    let method = match method {
        Ok(m) => m,
        Err(e) => {
            debug!("Method parsing error: {}", e);
            return Ok(StringChunked::from_slice(
                PlSmallStr::from_str("assembled_sequences"),
                &[""]
            ).into_series());
        }
    };

    // Extract sequences from input series
    let ca = inputs[0].str()?;
    
    // Convert string chunk to Vec<String>
    let sequences: Vec<String> = ca.into_iter()
        .flatten()
        .map(|s| s.to_string())
        .collect();
    
    // Call assembly function with only_largest always set to true
    match assemble_sequences(
        sequences,
        kwargs.k,
        kwargs.min_coverage,
        method,
        kwargs.export_graphs,
        Some(true), // Hardcoded to true
        kwargs.min_length,
        kwargs.auto_k,
        kwargs.prefix,
    ) {
        Ok(contigs) => {
            let result = contigs.join("\n");
            Ok(StringChunked::from_slice(
                    PlSmallStr::from_str("assembled_sequences"),
                    &[result.as_str()]
            ).into_series())
        },
        Err(e) => {
            debug!("Assembly failed: {}", e);
            Ok(StringChunked::from_slice(
                    PlSmallStr::from_str("assembled_sequences"),
                    &[""]
            ).into_series())
        }
    }
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


