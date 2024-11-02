use anyhow::Result;
use std::collections::HashSet;
use serde::Deserialize;
use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
use crate::fracture::assemble_sequences;
use log::*;


#[derive(Deserialize)]
pub struct OptimizeParams {
    pub start_k: usize,
    pub start_min_coverage: usize,
    pub start_anchor: String,
    pub end_anchor: String,
    pub max_iterations: Option<usize>,
    pub explore_k: Option<bool>,
}

// Types module
mod types {
    
    #[derive(Debug, Clone, Copy, Hash, Eq, PartialEq)]
    pub struct ParamPoint {
        pub k: usize,
        pub min_coverage: usize,
    }

    #[derive(Debug, Clone)]
    pub struct AssemblyResult {
        pub contig: String,
        pub params: ParamPoint,
        pub length: usize,
        pub has_anchors: bool,
        pub input_sequences: usize,  // Added field for input sequence count
    }


    impl AssemblyResult {
        pub fn new(contig: String, k: usize, min_coverage: usize, start_anchor: &str, end_anchor: &str, input_sequences: usize) -> Self {
            let length = contig.len();
            let has_anchors = contig.contains(start_anchor) && contig.contains(end_anchor);
            Self {
                contig,
                params: ParamPoint { k, min_coverage },
                length,
                has_anchors,
                input_sequences,
            }
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub enum Direction {
        West,
        East,
        North,
        South,
    }

    impl Direction {
        pub fn apply(&self, point: ParamPoint) -> Option<ParamPoint> {
            match self {
                Direction::West => {
                    if point.min_coverage > 1 {
                        Some(ParamPoint { k: point.k, min_coverage: point.min_coverage - 1 })
                    } else {
                        None
                    }
                },
                Direction::East => Some(ParamPoint { 
                    k: point.k, 
                    min_coverage: point.min_coverage + 1 
                }),
                Direction::North => {
                    if point.k > 4 {
                        Some(ParamPoint { k: point.k - 1, min_coverage: point.min_coverage })
                    } else {
                        None
                    }
                },
                Direction::South => {
                    if point.k < 64 {
                        Some(ParamPoint { k: point.k + 1, min_coverage: point.min_coverage })
                    } else {
                        None
                    }
                },
            }
        }
    }
}

use types::*;

pub fn optimize_assembly(
    sequences: &[String],
    params: ParamPoint,
    start_anchor: &str,
    end_anchor: &str,
    max_iterations: usize,
    explore_k: bool,
) -> Result<Option<AssemblyResult>> {
    info!("Starting assembly optimization with:");
    info!("  Initial k: {}, min_coverage: {}", params.k, params.min_coverage);
    info!("  Max iterations: {}, explore_k: {}", max_iterations, explore_k);
    info!("  Start anchor: {}", start_anchor);
    info!("  End anchor: {}", end_anchor);
    info!("  Number of input sequences: {}", sequences.len());

    let mut tested_params = HashSet::new();
    tested_params.insert(params);
    
    let current = assemble_and_check(sequences, params, start_anchor, end_anchor, sequences.len())?;
    debug!("Initial assembly result:");
    debug!("  Contig length: {}", current.length);
    debug!("  Has anchors: {}", current.has_anchors);
    
    if current.has_anchors {
        info!("Found solution on initial attempt!");
        return Ok(Some(current));
    }

    // Track current best length and all parameter points that achieved it
    let mut best_length = current.length;
    let mut current_points = vec![current.params];
    
    for iteration in 0..max_iterations {
        debug!("\nIteration {}/{}", iteration + 1, max_iterations);
        let mut candidates = Vec::new();
        
        let directions = if explore_k {
            vec![Direction::West, Direction::East, Direction::North, Direction::South]
        } else {
            vec![Direction::West, Direction::East]
        };
        
        // Try all directions from all current points
        for &point in &current_points {
            debug!("Exploring from point k={}, min_coverage={}", point.k, point.min_coverage);
            
            for direction in &directions {
                debug!("Trying direction: {:?}", direction);
                if let Some(new_params) = direction.apply(point) {
                    if !tested_params.contains(&new_params) {
                        debug!("Testing new params - k: {}, min_coverage: {}", new_params.k, new_params.min_coverage);
                        tested_params.insert(new_params);
                        
                        let result = assemble_and_check(sequences, new_params, start_anchor, end_anchor, sequences.len())?;
                        debug!("  Result - length: {}, has_anchors: {}", result.length, result.has_anchors);
                        
                        if result.has_anchors {
                            info!("Found solution at iteration {}!", iteration + 1);
                            info!("Final parameters: k={}, min_coverage={}", new_params.k, new_params.min_coverage);
                            return Ok(Some(result));
                        }
                        candidates.push(result);
                    } else {
                        debug!("Parameters already tested, skipping");
                    }
                } else {
                    debug!("Invalid parameters for direction {:?}, skipping", direction);
                }
            }
        }
        
        if candidates.is_empty() {
            info!("No more valid parameter combinations to try");
            break;
        }

        // Find best length among new candidates
        let max_length = candidates.iter()
            .map(|r| r.length)
            .max()
            .unwrap();

        if max_length > best_length {
            // If we found better solutions, only keep those points
            debug!("Found better length: {} > {}", max_length, best_length);
            best_length = max_length;
            current_points = candidates.iter()
                .filter(|r| r.length == max_length)
                .map(|r| r.params)
                .collect();
            
            debug!("New best points:");
            for point in &current_points {
                debug!("  k: {}, min_coverage: {}", point.k, point.min_coverage);
            }
        } else if max_length == best_length {
            // If we found equal solutions, add those points
            debug!("Found equal length solutions, adding to exploration points");
            current_points.extend(
                candidates.iter()
                    .filter(|r| r.length == max_length)
                    .map(|r| r.params)
            );
            debug!("Current exploration points:");
            for point in &current_points {
                debug!("  k: {}, min_coverage: {}", point.k, point.min_coverage);
            }
        }
        // If max_length < best_length, we keep our current points and continue exploring from them
    }
    
    info!("Optimization completed without finding anchors");
    info!("Parameters tested: {}", tested_params.len());
    info!("Best result: length={}", best_length);
    info!("Best points:");
    for point in &current_points {
        info!("  k={}, min_coverage={}", point.k, point.min_coverage);
    }
    
    // Return result with parameters from first current point
    // (all points have same length at this point)
    let final_result = assemble_and_check(sequences, current_points[0], start_anchor, end_anchor, sequences.len())?;
    Ok(Some(final_result))
}

fn assemble_and_check(
    sequences: &[String], 
    params: ParamPoint,
    start_anchor: &str,
    end_anchor: &str,
    input_sequences: usize,
) -> Result<AssemblyResult> {
    debug!("Attempting assembly with k={}, min_coverage={}", params.k, params.min_coverage);
    
    let contigs = assemble_sequences(
        sequences.to_vec(),
        params.k,
        params.min_coverage,
        None,
        Some(true), // only_largest
        None,
        None,
        None,
    )?;
    
    let contig = if contigs.is_empty() {
        debug!("No contigs produced");
        String::new()
    } else {
        debug!("Produced contig of length {}", contigs[0].len());
        contigs[0].clone()
    };
    
    Ok(AssemblyResult::new(
        contig,
        params.k,
        params.min_coverage,
        start_anchor,
        end_anchor,
        input_sequences
    ))
}

#[polars_expr(output_type_func=output_type)]
pub fn optimize_assembly_expr(inputs: &[Series], kwargs: OptimizeParams) -> PolarsResult<Series> {
    let _ = env_logger::try_init();

    let ca = inputs[0].str()?;
    let sequences: Vec<String> = ca.into_iter()
        .flatten()
        .map(|s| s.to_string())
        .collect();
        
    let sequence_count = sequences.len();
    let max_iterations = kwargs.max_iterations.unwrap_or(50);
    let explore_k = kwargs.explore_k.unwrap_or(false);
    
    let start_params = ParamPoint {
        k: kwargs.start_k,
        min_coverage: kwargs.start_min_coverage,
    };
    
    info!("Starting assembly optimization expression");
    info!("Input sequences: {}", sequence_count);
    
    match optimize_assembly(
        &sequences,
        start_params,
        &kwargs.start_anchor,
        &kwargs.end_anchor,
        max_iterations,
        explore_k,
    ) {
        Ok(Some(result)) => {
            info!("Optimization successful");
            info!("Final contig length: {}", result.length);
            
            let df = DataFrame::new(vec![
                Series::new("contig".into(), vec![result.contig]),
                Series::new("k".into(), vec![result.params.k as u32]),
                Series::new("min_coverage".into(), vec![result.params.min_coverage as u32]),
                Series::new("length".into(), vec![result.length as u32]),
                Series::new("input_sequences".into(), vec![result.input_sequences as u32]),
            ])?;
            
            Ok(df.into_struct(inputs[0].name().clone()).into())
        },
        Ok(None) => {
            info!("Optimization completed without finding anchors");
            
            let df = DataFrame::new(vec![
                Series::new("contig".into(), vec!["" as &str]),
                Series::new("k".into(), vec![0u32]),
                Series::new("min_coverage".into(), vec![0u32]),
                Series::new("length".into(), vec![0u32]),
                Series::new("input_sequences".into(), vec![sequence_count as u32]),
            ])?;
            
            Ok(df.into_struct(inputs[0].name().clone()).into())
        },
        Err(e) => {
            error!("Assembly optimization failed: {}", e);
            Err(PolarsError::ComputeError(
                format!("Assembly optimization failed: {}", e).into()
            ))
        },
    }
}

fn output_type(input_fields: &[Field]) -> PolarsResult<Field> {
    let fields = vec![
        Field::new("contig".into(), DataType::String),
        Field::new("k".into(), DataType::UInt32),
        Field::new("min_coverage".into(), DataType::UInt32),
        Field::new("length".into(), DataType::UInt32),
        Field::new("input_sequences".into(), DataType::UInt32),
    ];
    let struct_type = DataType::Struct(fields);
    Ok(Field::new(input_fields[0].name().clone(), struct_type))
}
