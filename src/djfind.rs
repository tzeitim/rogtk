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
    }
}

impl AssemblyMethod {
    /// Parse a method string and optional anchors into an AssemblyMethod
    pub fn from_str(method: &str, start_anchor: Option<String>, end_anchor: Option<String>) 
        -> Result<Self, String> {
        match method {
            "compression" => {
                // For compression method, anchors should not be provided
                if start_anchor.is_some() || end_anchor.is_some() {
                    return Err("Anchor sequences should not be provided for compression method".to_string());
                }
                Ok(AssemblyMethod::Compression)
            },
            "shortest_path" => {
                // For shortest path, both anchors must be provided
                match (start_anchor, end_anchor) {
                    (Some(start), Some(end)) => Ok(AssemblyMethod::ShortestPath { 
                        start_anchor: start, 
                        end_anchor: end 
                    }),
                    (None, None) => {
                        Err("Both start_anchor and end_anchor are required for shortest_path method".to_string())
                    },
                    (None, Some(_)) => {
                        Err("start_anchor is required for shortest_path method".to_string())
                    },
                    (Some(_), None) => {
                        Err("end_anchor is required for shortest_path method".to_string())
                    }
                }
            },
            _ => Err("Invalid assembly method. Must be 'compression' or 'shortest_path'".to_string())
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
        if sequence.contains(start_seq) {
            start_nodes.push(node_idx);
        }
        if sequence.contains(end_seq) {
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
                                trace!("  Checking neighbor: {} with distance {}", graph[neighbor], neighbor_dist);
                                if neighbor_dist < best_dist {
                                    best_dist = neighbor_dist;
                                    best_prev = Some(neighbor);
                                    trace!("    New best neighbor: {} with dist {}", graph[neighbor], neighbor_dist);
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

