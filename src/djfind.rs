use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::algo::dijkstra;
use petgraph::Incoming;
use std::collections::HashMap;
use debruijn::*;
use debruijn::graph::DebruijnGraph;
use log::*;
use std::fmt::Debug;
use serde::Deserialize;

#[allow(dead_code)]
pub struct PathFindingResult {
    pub path: Vec<String>,
    pub total_weight: f64,
    pub mean_coverage: f64,
    pub assembled_sequence: String,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AssemblyMethod {
    Compression,
    ShortestPath {
        start_anchor: String,
        end_anchor: String,
    },
    ShortestPathAuto,

}

impl AssemblyMethod {
    pub fn from_str(method: &str, start_anchor: Option<String>, end_anchor: Option<String>) -> Result<Self, String> {
        match method {
            "compression" => {
                if start_anchor.is_some() || end_anchor.is_some() {
                    return Err("Anchor sequences should not be provided for compression method".to_string());
                }
                Ok(AssemblyMethod::Compression)
            },
            "shortest_path" => {
                match (start_anchor, end_anchor) {
                    (Some(start), Some(end)) => Ok(AssemblyMethod::ShortestPath { 
                        start_anchor: start, 
                        end_anchor: end 
                    }),
                    _ => Err("Both start_anchor and end_anchor are required for shortest_path method".to_string()),
                }
            },
            "shortest_path_auto" => {
                if start_anchor.is_some() || end_anchor.is_some() {
                    return Err("Anchor sequences should not be provided for shortest_path_auto method".to_string());
                }
                Ok(AssemblyMethod::ShortestPathAuto)
            },
            _ => Err(format!("Unknown assembly method: {}", method)),
        }
    }
}

fn concatenate_path_sequences(sequences: &[String], k: usize) -> String {
    if sequences.is_empty() {
        return String::new();
    }
    
    let mut final_sequence = sequences[0].clone();
    
    // For each subsequent sequence
    for next_seq in sequences.iter().skip(1) {
        // Add only the non-overlapping part of the next sequence
        final_sequence.push_str(&next_seq[k - 1..]);
    }
    
    final_sequence
}


// Convert DebruijnGraph to petgraph::Graph for path finding
fn convert_to_petgraph<K: Kmer, D: std::fmt::Debug + Copy>(
    graph: &DebruijnGraph<K, D>,
) -> (DiGraph<String, f64>, HashMap<NodeIndex, usize>) {
    let mut pg = DiGraph::new();
    let mut node_map = HashMap::new();
    
    // First pass: create nodes
    for node_id in 0..graph.len() {
        let node = graph.get_node(node_id);
        let sequence = node.sequence().to_string();
        let idx = pg.add_node(sequence);
        node_map.insert(idx, node_id);
    }
    
    // Second pass: create edges with weights
    let node_indices: Vec<NodeIndex> = pg.node_indices().collect();
    for from_idx in &node_indices {
        let from_id = node_map[from_idx];
        let from_node = graph.get_node(from_id);
        let from_coverage = match from_node.data() {
            d => format!("{:?}", d).parse::<f64>().unwrap_or(1.0)
        };
        
        // Add edges for right connections
        for (to_id, _, _) in from_node.r_edges() {
            let to_node = graph.get_node(to_id);
            let to_coverage = match to_node.data() {
                d => format!("{:?}", d).parse::<f64>().unwrap_or(1.0)
            };
            
            // Find corresponding petgraph NodeIndex for to_id
            if let Some(to_idx) = node_indices.iter().find(|&&idx| node_map[&idx] == to_id) {
                // Weight is inverse of mean coverage between nodes
                //let weight = 2.0 / (from_coverage + to_coverage);
                // negative log
                let mean_coverage = (from_coverage + to_coverage) / 2.0;
                let weight = -mean_coverage.ln();
                pg.add_edge(*from_idx, *to_idx, weight);
            }
        }
    }
    
    (pg, node_map)
}

// Find nodes containing specific sequences
fn find_anchor_nodes(
    graph: &DiGraph<String, f64>,
    start_seq: &str,
    end_seq: &str,
) -> (Vec<NodeIndex>, Vec<NodeIndex>) {
    let mut start_nodes = Vec::new();
    let mut end_nodes = Vec::new();
    
    for node_idx in graph.node_indices() {
        let sequence = &graph[node_idx];
        if sequence.starts_with(start_seq) {
            debug!("Found start node idx {} {}", node_idx.index(), &graph[node_idx]);
            start_nodes.push(node_idx);
        }
        if sequence.ends_with(end_seq) {
            debug!("Found end node idx {} {}", node_idx.index(), &graph[node_idx]);
            end_nodes.push(node_idx);
        }
    }
    
    debug!("Found {} start nodes containing: {}", 
    start_nodes.len(), 
    start_nodes.iter().map(|&n| graph[n].as_str()).collect::<Vec<_>>().join(", ")
);
debug!("Found {} end nodes containing: {}", 
    end_nodes.len(),
    end_nodes.iter().map(|&n| graph[n].as_str()).collect::<Vec<_>>().join(", ")
);

    (start_nodes, end_nodes)
}

// Find shortest path between anchors
pub fn find_shortest_path(
    graph: &DiGraph<String, f64>,
    start_nodes: &[NodeIndex],
    end_nodes: &[NodeIndex],
) -> Option<(Vec<NodeIndex>, f64)> {
    let mut best_path = None;
    let mut min_total_weight = f64::INFINITY;
    const MAX_ITERATIONS: usize = 1000; // Safety limit
    const FLOAT_EPSILON: f64 = 1e-9; // For floating-point comparison in path reconstruction

    for &start in start_nodes {
        let distances = dijkstra(graph, start, None, |e| *e.weight());
        
        for &end in end_nodes {
            if let Some(&total_weight) = distances.get(&end) {
                debug!("Found path to end {} with weight {}", graph[end], total_weight);
                
                if total_weight < min_total_weight {
                    let mut path = vec![end];
                    let mut current = end;
                    let mut path_valid = false;
                    let mut iterations = 0;
                    
                    while current != start {
                        iterations += 1;
                        if iterations > MAX_ITERATIONS {
                            trace!("Hit iteration limit while searching path");
                            break;
                        }

                        trace!("Iteration {}: Current node: {}", iterations, graph[current]);
                        let mut best_prev = None;
                        let mut best_dist = f64::INFINITY;
                        let current_dist = distances[&current];
                        
                        trace!("Current distance: {}", current_dist);
                        
                        for neighbor in graph.neighbors_directed(current, Incoming) {
                            if let Some(&neighbor_dist) = distances.get(&neighbor) {
                                // Get the edge weight from neighbor to current
                                if let Some(edge_idx) = graph.find_edge(neighbor, current) {
                                    let edge_weight = graph[edge_idx];
                                    let expected_dist = neighbor_dist + edge_weight;

                                    trace!("  Checking neighbor: {} dist={}, edge_weight={}, expected_dist={}",
                                           graph[neighbor], neighbor_dist, edge_weight, expected_dist);

                                    // Check if this neighbor is on the shortest path (within float epsilon)
                                    if (expected_dist - current_dist).abs() < FLOAT_EPSILON {
                                        // Use neighbor_dist as tie-breaker for determinism
                                        if neighbor_dist < best_dist {
                                            best_dist = neighbor_dist;
                                            best_prev = Some(neighbor);
                                            trace!("    Valid predecessor on shortest path: {}", graph[neighbor]);
                                        }
                                    }
                                }
                            }
                        }
                        
                        match best_prev {
                            Some(prev) => {
                                trace!("Selected best previous node: {} with dist {}", graph[prev], best_dist);
                                path.push(prev);
                                current = prev;
                                if current == start {
                                    path_valid = true;
                                }
                            }
                            None => {
                                trace!("No valid previous node found, breaking");
                                break;
                            }
                        }
                    }
                    
                    if path_valid {
                        debug!("Found valid path with {} nodes after {} iterations", path.len(), iterations);
                        path.reverse();
                        best_path = Some((path, total_weight));
                        min_total_weight = total_weight;
                    } else {
                        debug!("Path reconstruction failed after {} iterations", iterations);
                    }
                }
            }
        }
    }

    best_path
}
// Extract path sequences
pub fn extract_path_sequences(
    graph: &DiGraph<String, f64>,
    path: &[NodeIndex],
) -> Vec<String> {
    path.iter().map(|&idx| graph[idx].clone()).collect()
}

// Main function to perform path-based assembly
pub fn assemble_with_path_finding<K: Kmer + Send + Sync + Debug + 'static>(
    preliminary_graph: &DebruijnGraph<K, u16>,
    start_anchor: &str,
    end_anchor: &str,
) -> Result<PathFindingResult, String> {
    let _ = env_logger::try_init();
    info!("Converting de Bruijn graph to weighted digraph for path finding");

    // Convert to petgraph
    let (pg, _node_map) = convert_to_petgraph(preliminary_graph);

    // Find anchor nodes
    let (start_nodes, end_nodes) = find_anchor_nodes(&pg, start_anchor, end_anchor);

    if start_nodes.is_empty() {
        return Err(format!("No nodes containing start anchor '{}' found", start_anchor));
    }
    if end_nodes.is_empty() {
        return Err(format!("No nodes containing end anchor '{}' found", end_anchor));
    }

    info!("Found {} start nodes and {} end nodes", start_nodes.len(), end_nodes.len());

    // Find shortest path
    if let Some((path, total_weight)) = find_shortest_path(&pg, &start_nodes, &end_nodes) {
        let sequences = extract_path_sequences(&pg, &path);
        let mean_coverage = 1.0 / (total_weight / path.len() as f64);

        let assembled_sequence = concatenate_path_sequences(&sequences, K::k());

        info!("Found path with total weight {} and mean coverage {}", 
            total_weight, mean_coverage);

        info!("Assembled sequence length: {}", assembled_sequence.len());
        warn!("Assembled sequence : {}", assembled_sequence);


        Ok(PathFindingResult {
            path: sequences,
            total_weight,
            mean_coverage,
            assembled_sequence,
        })

    } else {
        Err("No valid path found between anchors".to_string())
    }
}
// Add this function to djfind.rs for automatic endpoint detection and path finding

/// Find candidate endpoints based on in/out degree analysis
/// Returns (start_candidates, end_candidates)
fn find_endpoint_candidates<K: Kmer>(
    graph: &DebruijnGraph<K, u16>,
) -> (Vec<usize>, Vec<usize>) {
    let mut start_candidates = Vec::new();
    let mut end_candidates = Vec::new();
    
    // Calculate average coverage for filtering
    let mut coverages = Vec::new();
    for node_id in 0..graph.len() {
        let node = graph.get_node(node_id);
        coverages.push(*node.data());
    }
    
    let avg_coverage = coverages.iter().map(|&c| c as f64).sum::<f64>() / coverages.len() as f64;
    let min_coverage_threshold = (avg_coverage * 0.1).max(1.0) as u16;
    
    info!("Average coverage: {:.2}, min threshold: {}", avg_coverage, min_coverage_threshold);
    
    for node_id in 0..graph.len() {
        let node = graph.get_node(node_id);
        let coverage = *node.data();
        
        // Skip very low coverage nodes (likely errors)
        if coverage < min_coverage_threshold {
            continue;
        }
        
        // Count in-edges and out-edges separately
        let in_degree = node.l_edges().len();
        let out_degree = node.r_edges().len();
        
        debug!("Node {:?}: coverage={}, in={}, out={}", 
               node.sequence(), coverage, in_degree, out_degree);
        
        // Start candidates: no incoming edges (or very few compared to outgoing)
        if in_degree == 0 && out_degree > 0 {
            start_candidates.push(node_id);
            info!("Start candidate found: {:?} (coverage={})", node.sequence(), coverage);
        }
        
        // End candidates: no outgoing edges (or very few compared to incoming)
        if out_degree == 0 && in_degree > 0 {
            end_candidates.push(node_id);
            info!("End candidate found: {:?} (coverage={})", node.sequence(), coverage);
        }
    }
    
    (start_candidates, end_candidates)
}

/// Score a path based on multiple metrics
fn score_path(
    graph: &DiGraph<String, f64>,
    path: &[NodeIndex],
    total_weight: f64,
) -> f64 {
    if path.is_empty() {
        return 0.0;
    }
    
    // Calculate path length
    let path_length = path.iter()
        .map(|&idx| graph[idx].len())
        .sum::<usize>() as f64;
    
    // Calculate mean coverage (inverse of weight)
    let mean_coverage = 1.0 / (total_weight / path.len() as f64);
    
    // Simple scoring: prioritize length and coverage equally
    // Normalize length by expected amplicon size (5000bp)
    let normalized_length = (path_length / 5000.0).min(1.0);
    let normalized_coverage = (mean_coverage / 100.0).min(1.0);
    
    let score = 0.6 * normalized_length + 0.4 * normalized_coverage;
    
    debug!("Path score: {:.3} (len={:.0}, cov={:.1})", 
           score, path_length, mean_coverage);
    
    score
}

/// Find best path among multiple endpoint pairs
fn find_best_endpoint_pair<K: Kmer + Send + Sync + Debug + 'static>(
    graph: &DebruijnGraph<K, u16>,
    start_candidates: Vec<usize>,
    end_candidates: Vec<usize>,
) -> Result<PathFindingResult, String> {
    info!("Evaluating {} start x {} end candidate pairs", 
         start_candidates.len(), end_candidates.len());
    
    // Convert to petgraph once
    let (pg, _) = convert_to_petgraph(graph);
    
    // Limit number of pairs to evaluate
    const MAX_PAIRS: usize = 100;
    let mut evaluated_pairs = 0;
    
    let mut best_result: Option<(PathFindingResult, f64)> = None;
    
    for &start_id in &start_candidates {
        for &end_id in &end_candidates {
            if evaluated_pairs >= MAX_PAIRS {
                warn!("Reached maximum pairs limit ({})", MAX_PAIRS);
                break;
            }
            evaluated_pairs += 1;
            
            let start_seq = graph.get_node(start_id).sequence().to_string();
            let end_seq = graph.get_node(end_id).sequence().to_string();
            
            debug!("Evaluating path: {} -> {}", start_seq, end_seq);
            
            // Find nodes in petgraph
            let start_nodes: Vec<NodeIndex> = pg.node_indices()
                .filter(|&idx| pg[idx].contains(&start_seq))
                .collect();
            let end_nodes: Vec<NodeIndex> = pg.node_indices()
                .filter(|&idx| pg[idx].contains(&end_seq))
                .collect();
            
            if start_nodes.is_empty() || end_nodes.is_empty() {
                continue;
            }
            
            // Try to find path
            if let Some((path, total_weight)) = find_shortest_path(&pg, &start_nodes, &end_nodes) {
                let score = score_path(&pg, &path, total_weight);
                
                info!("Found path with score {:.3}", score);
                
                if best_result.is_none() || score > best_result.as_ref().unwrap().1 {
                    let sequences = extract_path_sequences(&pg, &path);
                    let mean_coverage = 1.0 / (total_weight / path.len() as f64);
                    let assembled_sequence = concatenate_path_sequences(&sequences, K::k());
                    
                    best_result = Some((PathFindingResult {
                        path: sequences,
                        total_weight,
                        mean_coverage,
                        assembled_sequence,
                    }, score));
                }
            }
        }
    }
    
    match best_result {
        Some((result, score)) => {
            info!("Selected best path with score {:.3}, length {}", 
                  score, result.assembled_sequence.len());
            Ok(result)
        }
        None => Err("No valid paths found between any endpoint pairs".to_string())
    }
}

/// Main entry point for automatic path-based assembly
pub fn assemble_with_auto_path_finding<K: Kmer + Send + Sync + Debug + 'static>(
    preliminary_graph: &DebruijnGraph<K, u16>,
) -> Result<PathFindingResult, String> {
    let _ = env_logger::try_init();
    info!("Starting automatic path finding assembly");
    
    // Find endpoint candidates
    let (start_candidates, end_candidates) = find_endpoint_candidates(preliminary_graph);
    
    match (start_candidates.len(), end_candidates.len()) {
        (0, _) => Err("No start candidates found - possibly circular or highly branched".to_string()),
        (_, 0) => Err("No end candidates found - possibly circular or highly branched".to_string()),
        (1, 1) => {
            // Ideal case - single path
            info!("Found exactly one start and one end candidate");
            let start_seq = preliminary_graph.get_node(start_candidates[0]).sequence().to_string();
            let end_seq = preliminary_graph.get_node(end_candidates[0]).sequence().to_string();
            assemble_with_path_finding(preliminary_graph, &start_seq, &end_seq)
        },
        _ => {
            // Multiple candidates - need to evaluate
            info!("Found {} start and {} end candidates", 
                  start_candidates.len(), end_candidates.len());
            find_best_endpoint_pair(preliminary_graph, start_candidates, end_candidates)
        }
    }
}
