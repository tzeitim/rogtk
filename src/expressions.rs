use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
use std::fmt::Write;

fn pig_latin_str(value: &str, output: &mut String) {
    if let Some(first_char) = value.chars().next() {
        write!(output, "{}{}ay", &value[1..], first_char).unwrap()
    }
}

#[polars_expr(output_type=String)]
fn pig_latinnify(inputs: &[Series]) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let out: StringChunked = ca.apply_to_buffer(pig_latin_str);
    Ok(out.into_series())
}

fn parse_cigar_str(cigar: &str, output: &mut String) {
    let mut num_buf = String::new();
    let mut ref_pos = 0;

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buf.push(c);
        } else {
            if let Ok(len) = num_buf.parse::<usize>() {
                match c {
                    'D' => {
                        for del_pos in ref_pos..(ref_pos + len) {
                            output.push_str(&format!("D,{},1,", del_pos));
                        }
                        ref_pos += len;
                    }
                    'I' => {
                        output.push_str(&format!("I,{},{},", ref_pos, len));
                    }
                    _ => {
                        ref_pos += len;
                    }
                }
            }
            num_buf.clear();
        }
    }

    // Remove the trailing comma, if present
    if output.ends_with(',') {
        output.pop();
    }
}

#[polars_expr(output_type=String)]
fn parse_cigar_series(inputs: &[Series]) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let out: StringChunked = ca.apply_to_buffer(parse_cigar_str);
    Ok(out.into_series())
}

