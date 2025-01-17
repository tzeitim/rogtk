use anyhow::Result;
use std::collections::HashSet;
use serde::Deserialize;
use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
use crate::fracture::assemble_sequences;
use crate::djfind::AssemblyMethod;
use log::*;

#[derive(Deserialize)]
#[allow(dead_code)]
pub struct OptimizeParams {
    pub start_anchor: Option<String>,
    pub end_anchor: Option<String>,

    pub method: String,
    pub min_length: Option<usize>,
    pub export_graphs: Option<bool>,
    pub prefix: Option<String>,

    pub start_k: usize,
    pub start_min_coverage: usize,
    pub max_iterations: Option<usize>,
    pub explore_k: Option<bool>,
    pub prioritize_length: Option<bool>,
}

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

#[derive(Debug, Clone)]
struct ExplorationPath {
    params: ParamPoint,
    length: usize,
    steps_without_improvement: usize,
    direction_history: Vec<Direction>,
}

pub fn optimize_assembly(
    sequences: &[String],
    params: ParamPoint,
    start_anchor: &str,
    end_anchor: &str,
    max_iterations: usize,
    explore_k: bool,
    prioritize_length: bool,
    method: AssemblyMethod,
) -> Result<Option<AssemblyResult>> {
    let mut tested_params = HashSet::new();
    tested_params.insert(params);

    // Track best results
    let mut best_anchored_result: Option<AssemblyResult> = None;
    let mut best_length_result: Option<AssemblyResult> = None;

    // Get initial result
    let current = assemble_and_check(sequences, params, start_anchor, end_anchor, sequences.len(), method.clone())?;
    
    // Update best results
    if current.has_anchors {
        debug!("Current seq {}", current.contig);
        best_anchored_result = Some(current.clone());
    }
    if best_length_result.as_ref().map_or(true, |r| current.length > r.length) {
        best_length_result = Some(current.clone());
    }

    // Initialize paths for exploration
    let mut paths = vec![ExplorationPath {
        params: current.params,
        length: current.length,
        steps_without_improvement: 0,
        direction_history: Vec::new(),
    }];

    let directions = if explore_k {
        vec![Direction::West, Direction::East, Direction::North, Direction::South]
    } else {
        vec![Direction::West, Direction::East]
    };

    for iteration in 0..max_iterations {
        debug!("\nIteration {}/{}", iteration + 1, max_iterations);
        let mut new_paths = Vec::new();

        for path in &paths {
            debug!("Exploring from point k={}, min_coverage={}", 
                  path.params.k, path.params.min_coverage);

            for direction in &directions {
                if let Some(new_params) = direction.apply(path.params) {
                    if !tested_params.contains(&new_params) {
                        tested_params.insert(new_params);
                        
                        let result = assemble_and_check(
                            sequences, new_params, start_anchor, end_anchor, sequences.len(), method.clone())?;

                        // Update best results
                        if result.has_anchors && best_anchored_result.as_ref().map_or(true, |r| result.length > r.length) {
                            best_anchored_result = Some(result.clone());
                        }
                        if best_length_result.as_ref().map_or(true, |r| result.length > r.length) {
                            best_length_result = Some(result.clone());
                        }

                        // Early return if we find anchored solution and aren't prioritizing length
                        if result.has_anchors && !prioritize_length {
                            info!("Found solution with both anchors at iteration {}!", iteration + 1);
                            return Ok(Some(result));
                        }

                        // Create new path if result is promising
                        if !result.contig.is_empty() {
                            let mut new_direction_history = path.direction_history.clone();
                            new_direction_history.push(*direction);

                            let steps_without_improvement = 
                                if result.length > path.length { 0 } else { path.steps_without_improvement + 1 };

                            new_paths.push(ExplorationPath {
                                params: new_params,
                                length: result.length,
                                steps_without_improvement,
                                direction_history: new_direction_history,
                            });
                        }
                    }
                }
            }
        }

        if new_paths.is_empty() {
            break;
        }

        paths = select_promising_paths(new_paths);
    }

    // Return best result based on priority
    if prioritize_length {
        info!("Optimization completed, returning longest contig");
        Ok(best_length_result)
    } else {
        info!("Optimization completed, returning best anchored contig");
        Ok(best_anchored_result)
    }
}

fn select_promising_paths(mut paths: Vec<ExplorationPath>) -> Vec<ExplorationPath> {
    // Sort paths by length and recency of improvement
    paths.sort_by(|a, b| {
        b.length.cmp(&a.length)
            .then(a.steps_without_improvement.cmp(&b.steps_without_improvement))
    });

    // Keep top N paths (adjust N based on desired exploration breadth)
    const MAX_ACTIVE_PATHS: usize = 4;
    paths.truncate(MAX_ACTIVE_PATHS);
    paths
}

fn assemble_and_check(
    sequences: &[String], 
    params: ParamPoint,
    start_anchor: &str,
    end_anchor: &str,
    input_sequences: usize,
    method: AssemblyMethod,
) -> Result<AssemblyResult> {
    debug!("Attempting assembly with k={}, min_coverage={}", params.k, params.min_coverage);
    
    let contigs = assemble_sequences(
        sequences.to_vec(),
        params.k,
        params.min_coverage,
        method,
        Some(false), // don't export
        Some(true),  // only_largest
        None,
        None,
        None,
    )?;
    
    let contig = if contigs.is_empty() {
        debug!("No contigs produced");
        String::new()
    } else {
        let length = contigs[0].len();
        debug!("Produced contig of length {}", length);
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
    let prioritize_length = kwargs.prioritize_length.unwrap_or(false);
    
    let start_anchor = kwargs.start_anchor.ok_or_else(||
        PolarsError::ComputeError("start_anchor is required".into()))?;

    let end_anchor = kwargs.end_anchor.ok_or_else(||
        PolarsError::ComputeError("end_anchor is required".into()))?;

    let assembly_method = AssemblyMethod::from_str(
        &kwargs.method,
        Some(start_anchor.clone()),  
        Some(end_anchor.clone())
    ).map_err(|e| PolarsError::ComputeError(e.into()))?;

    let start_params = ParamPoint {
        k: kwargs.start_k,
        min_coverage: kwargs.start_min_coverage,
    };
    
    info!("Starting assembly optimization expression");
    info!("Input sequences: {}", sequence_count);
    
    match optimize_assembly(
        &sequences,
        start_params,
        &start_anchor,
        &end_anchor,
        max_iterations,
        explore_k,
        prioritize_length,
        assembly_method,
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
            info!("No valid contig found");
            
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
