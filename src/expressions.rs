use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
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

#[polars_expr(output_type=String)]
fn parse_cigar_series(inputs: &[Series], kwargs: ParseCigarKwargs) -> PolarsResult<Series> {
    let ca: &StringChunked = inputs[0].str()?;
    //let out: StringChunked = ca.apply_to_buffer(parse_cigar_str);
    let out: StringChunked = ca.apply_to_buffer(|value, output|{
            parse_cigar_str(value, output, kwargs.block_dels);
    });
    Ok(out.into_series())
}
