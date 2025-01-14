use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::algo::dijkstra;
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
                let weight = 2.0 / (from_coverage + to_coverage);
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
    
    (start_nodes, end_nodes)
}

// Find shortest path between anchors
pub fn find_shortest_path(
    graph: &DiGraph<String, f64>,
    start_nodes: &[NodeIndex],
    end_nodes: &[NodeIndex],
) -> Option<(Vec<NodeIndex>, f64)> {
    let _ = env_logger::try_init();
    info!("Starting shortest path search with {} start nodes and {} end nodes", 
          start_nodes.len(), end_nodes.len());
    
    let mut best_path = None;
    let mut min_total_weight = f64::INFINITY;
    
    for (i, &start) in start_nodes.iter().enumerate() {
        debug!("Analyzing start node {}/{}: {:?}", i + 1, start_nodes.len(), graph[start]);
        
        // Run Dijkstra's algorithm from this start node
        let distances = dijkstra(graph, start, None, |e| *e.weight());


        debug!("Completed Dijkstra's algorithm from start node {}, found {} reachable nodes", 
               i + 1, distances.len());
        
        // Check all possible end nodes
        for (j, &end) in end_nodes.iter().enumerate() {

            debug!("Distances computed for start node {} -> end node {}:", i + 1, j + 1);
            for (node, dist) in &distances {
                debug!("  Node {:?} -> Distance: {}", graph[*node], dist);
            }

            if let Some(weight) = distances.get(&end) {
                debug!("Found path to end node {}/{} with weight {}", 
                       j + 1, end_nodes.len(), weight);
                
                if *weight < min_total_weight {
                    debug!("New best path found with weight {} (previous best: {})", 
                           weight, min_total_weight);
                    
                    // Reconstruct path
                    let mut path = Vec::new();
                    let mut current = end;
                    path.push(current);
                    
                    while current != start {
                        let mut min_prev = None;
                        let mut min_weight = f64::INFINITY;
                        
                        for neighbor in graph.neighbors_directed(current, petgraph::Direction::Incoming) {
                            if let Some(dist) = distances.get(&neighbor) {
                                let edge_weight = graph.edge_weight(graph.find_edge(neighbor, current).unwrap()).unwrap();
                                if dist + edge_weight < min_weight {
                                    min_weight = dist + edge_weight;
                                    min_prev = Some(neighbor);
                                }
                            }
                        }
                        
                        if let Some(prev) = min_prev {
                            current = prev;
                            path.push(current);
                            trace!("Added node to path: {:?}", graph[current]);
                        } else {
                            warn!("Path reconstruction failed: couldn't find previous node");
                            break;
                        }
                    }
                    
                    path.reverse();
                    min_total_weight = *weight;
                    best_path = Some((path.clone(), *weight));
                    
                    info!("Updated best path: length={}, total_weight={}", 
                          path.len(), weight);
                }
            } else {
                debug!("No path found to end node {}/{}", j + 1, end_nodes.len());
            }
        }
    }
    
    match &best_path {
        Some((path, weight)) => {
            info!("Found optimal path: length={}, total_weight={}", 
                  path.len(), weight);
        }
        None => {
            warn!("No valid path found between any start and end nodes");
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

