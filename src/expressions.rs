use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
//use polars_core::series::amortized_iter::AmortSeries;
//use crate::expressions::polars_plan::prelude::lit;
use serde::Deserialize;
//use std::fmt::Write;

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
