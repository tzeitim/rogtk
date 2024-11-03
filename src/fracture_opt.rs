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
    use log::debug;
    
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
        pub input_sequences: usize,
    }

    impl AssemblyResult {
        pub fn new(contig: String, k: usize, min_coverage: usize, start_anchor: &str, end_anchor: &str, input_sequences: usize) -> Self {
            let has_start = contig.contains(start_anchor);
            let has_end = contig.contains(end_anchor);
            let has_anchors = has_start && has_end;

            if !has_start {
                debug!("Missing start anchor sequence");
            }
            if !has_end {
                debug!("Missing end anchor sequence");
            }

            Self {
                contig: contig.clone(),
                params: ParamPoint { k, min_coverage },
                length: contig.len(), // Keep actual length for optimization
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
    
    // Track best result that has both anchors
    let mut best_with_anchors: Option<AssemblyResult> = if current.has_anchors {
        info!("Found initial solution!");
        Some(current.clone())
    } else {
        None
    };

    // Start with current parameters
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
                        
                        // If we found a solution with both anchors
                        if result.has_anchors {
                            match &best_with_anchors {
                                None => {
                                    info!("Found first solution with both anchors at iteration {}!", iteration + 1);
                                    best_with_anchors = Some(result.clone());
                                }
                                Some(best) if result.length > best.length => {
                                    info!("Found better solution with both anchors: {} > {}", result.length, best.length);
                                    best_with_anchors = Some(result.clone());
                                }
                                _ => {}
                            }
                        }
                        
                        // Keep exploring from this point if it produced a contig
                        if !result.contig.is_empty() {
                            candidates.push(result);
                        }
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

        // Update exploration points based on promising directions
        // Prefer points that produced contigs with at least one anchor
        let candidates_with_start: Vec<_> = candidates.iter()
            .filter(|r| r.contig.contains(start_anchor) || r.contig.contains(end_anchor))
            .collect();

        if !candidates_with_start.is_empty() {
            // If we found contigs with anchors, prefer those directions
            current_points = candidates_with_start.iter()
                .map(|r| r.params)
                .collect();
        } else {
            // Otherwise, use length as a heuristic for promising directions
            let max_length = candidates.iter()
                .map(|r| r.length)
                .max()
                .unwrap();

            current_points = candidates.iter()
                .filter(|r| r.length == max_length)
                .map(|r| r.params)
                .collect();
        }
    }
    
    info!("Optimization completed after testing {} parameter combinations", tested_params.len());
    
    // Return the best result that had both anchors, if any
    match best_with_anchors {
        Some(result) => {
            info!("Found valid solution with both anchors");
            info!("Final parameters: k={}, min_coverage={}", result.params.k, result.params.min_coverage);
            info!("Final contig length: {}", result.length);
            Ok(Some(result))
        }
        None => {
            info!("No valid contig found containing both anchor sequences");
            Ok(None)
        }
    }
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
        Some(false), // don't export
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
        Ok(None) | Err(_) => {
            info!("No valid contig found with both anchor sequences");
            
            let df = DataFrame::new(vec![
                Series::new("contig".into(), vec!["" as &str]),
                Series::new("k".into(), vec![0u32]),
                Series::new("min_coverage".into(), vec![0u32]),
                Series::new("length".into(), vec![0u32]),
                Series::new("input_sequences".into(), vec![sequence_count as u32]),
            ])?;
            
            Ok(df.into_struct(inputs[0].name().clone()).into())
        }
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
