use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::algo::dijkstra;
use std::collections::HashMap;
use debruijn::*;
use debruijn::graph::DebruijnGraph;
use log::*;
use std::fmt::Debug;
use serde::Deserialize;

//#[allow(dead_code)]
//pub struct PathFindingResult {
//    pub path: Vec<String>,
//    pub total_weight: f64,
//    pub mean_coverage: f64,
//    pub assembled_sequence: String,
//}


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

// Custom struct to hold node data with Drop trait
struct NodeData {
    sequence: String,
    coverage: f64,
}

impl Drop for NodeData {
    fn drop(&mut self) {
        // Custom cleanup logic if needed
        self.sequence.clear();
    }
}

fn convert_to_petgraph<K: Kmer, D: std::fmt::Debug + Copy>(
    graph: &DebruijnGraph<K, D>,
) -> (DiGraph<NodeData, f64>, HashMap<NodeIndex, usize>) {
    // Pre-allocate collections with capacity
    let node_count = graph.len();
    let mut pg = DiGraph::with_capacity(node_count, node_count * 2);  // Estimate edge count
    let mut node_map = HashMap::with_capacity(node_count);
    
    // First pass: create nodes
    // Store indices in a pre-allocated Vec to avoid repeated allocations
    let mut node_indices = Vec::with_capacity(node_count);
    
    for node_id in 0..node_count {
        let node = graph.get_node(node_id);
        // Create NodeData struct that owns the sequence string
        let node_data = NodeData {
            sequence: node.sequence().to_string(),  // One-time allocation
            coverage: match node.data() {
                d => format!("{:?}", d).parse::<f64>().unwrap_or(1.0)
            },
        };
        
        let idx = pg.add_node(node_data);
        node_map.insert(idx, node_id);
        node_indices.push(idx);
    }
    
    // Second pass: create edges
    // Use references to avoid string allocations in edge processing
    for (idx_pos, &from_idx) in node_indices.iter().enumerate() {
        let from_id = node_map[&from_idx];
        let from_node = graph.get_node(from_id);
        let from_coverage = match from_node.data() {
            d => format!("{:?}", d).parse::<f64>().unwrap_or(1.0)
        };
        
        // Process right edges in a single pass
        for (to_id, _, _) in from_node.r_edges() {
            let to_node = graph.get_node(to_id);
            let to_coverage = match to_node.data() {
                d => format!("{:?}", d).parse::<f64>().unwrap_or(1.0)
            };
            
            // Use direct index lookup instead of find
            if let Some(&to_idx) = node_indices.get(to_id) {
                let mean_coverage = (from_coverage + to_coverage) / 2.0;
                let weight = -mean_coverage.ln();
                pg.add_edge(from_idx, to_idx, weight);
            }
        }
    }
    
    // Clear temporary data structure
    node_indices.clear();
    
    (pg, node_map)
}

// Helper trait for finding nodes with specific sequences
trait GraphSearchExt {
    fn find_nodes_containing(&self, sequence: &str) -> Vec<NodeIndex>;
}

impl GraphSearchExt for DiGraph<NodeData, f64> {
    fn find_nodes_containing(&self, sequence: &str) -> Vec<NodeIndex> {
        self.node_indices()
            .filter(|&idx| self[idx].sequence.contains(sequence))
            .collect()
    }
}

// Find nodes containing specific sequences
fn find_anchor_nodes(
    graph: &DiGraph<NodeData, f64>,
    start_seq: &str,
    end_seq: &str,
) -> (Vec<NodeIndex>, Vec<NodeIndex>) {
    let start_nodes = graph.find_nodes_containing(start_seq);
    let end_nodes = graph.find_nodes_containing(end_seq);
    (start_nodes, end_nodes)
}

// Helper function to estimate path capacity
fn estimate_path_capacity(
    distances: &HashMap<NodeIndex, f64>,
    graph: &DiGraph<NodeData, f64>
) -> usize {
    // Estimate based on average distance and graph density
    let avg_distance = distances.values().sum::<f64>() / distances.len() as f64;
    let graph_density = graph.edge_count() as f64 / (graph.node_count() * graph.node_count()) as f64;
    
    // Heuristic formula: higher density or distance suggests longer paths
    ((avg_distance * graph_density * 2.0) as usize).max(10).min(graph.node_count())
}

// Helper struct for path reconstruction to avoid repeated allocations
struct PathReconstructor<'a> {
    graph: &'a DiGraph<NodeData, f64>,
    distances: &'a HashMap<NodeIndex, f64>,
    path: Vec<NodeIndex>,
}

impl<'a> PathReconstructor<'a> {
    fn new(graph: &'a DiGraph<NodeData, f64>, distances: &'a HashMap<NodeIndex, f64>) -> Self {
        let capacity = estimate_path_capacity(distances, graph);
        Self {
            graph,
            distances,
            path: Vec::with_capacity(capacity),
        }
    }

    fn reconstruct(&mut self, start: NodeIndex, end: NodeIndex) -> Option<Vec<NodeIndex>> {
        self.path.clear();  // Reuse existing allocation
        self.path.push(end);
        
        let mut current = end;
        
        while current != start {
            let next = self.find_next_node(current)?;
            current = next;
            self.path.push(current);
        }
        
        self.path.reverse();
        Some(self.path.clone())
    }

    fn find_next_node(&self, current: NodeIndex) -> Option<NodeIndex> {
        let mut min_prev = None;
        let mut min_weight = f64::INFINITY;
        
        for neighbor in self.graph.neighbors_directed(current, petgraph::Direction::Incoming) {
            if let Some(dist) = self.distances.get(&neighbor) {
                let edge_weight = self.graph
                    .edge_weight(self.graph.find_edge(neighbor, current)?)?;
                if dist + edge_weight < min_weight {
                    min_weight = dist + edge_weight;
                    min_prev = Some(neighbor);
                }
            }
        }
        
        min_prev
    }
}

pub fn find_shortest_path(
    graph: &DiGraph<NodeData, f64>,
    start_nodes: &[NodeIndex],
    end_nodes: &[NodeIndex],
) -> Option<(Vec<NodeIndex>, f64)> {
    let _ = env_logger::try_init();
    info!("Starting shortest path search with {} start nodes and {} end nodes", 
          start_nodes.len(), end_nodes.len());
    
    let mut best_path = None;
    let mut min_total_weight = f64::INFINITY;
    
    for (i, &start) in start_nodes.iter().enumerate() {
        debug!("Analyzing start node {}/{}: {:?}", i + 1, start_nodes.len(), graph[start].sequence);
        
        {   // Block for distances and path reconstruction scope
            let distances = dijkstra(graph, start, None, |e| *e.weight());
            let mut reconstructor = PathReconstructor::new(graph, &distances);
            
            debug!("Completed Dijkstra's algorithm from start node {}, found {} reachable nodes", 
                   i + 1, distances.len());
            
            for (j, &end) in end_nodes.iter().enumerate() {
                if let Some(weight) = distances.get(&end) {
                    debug!("Found path to end node {}/{} with weight {}", 
                           j + 1, end_nodes.len(), weight);
                    
                    if *weight < min_total_weight {
                        debug!("New best path found with weight {} (previous best: {})", 
                               weight, min_total_weight);
                        
                        if let Some(path) = reconstructor.reconstruct(start, end) {
                            min_total_weight = *weight;
                            best_path = Some((path, *weight));
                            
                            info!("Updated best path: length={}, total_weight={}", 
                                  best_path.as_ref().unwrap().0.len(), weight);
                        } else {
                            warn!("Path reconstruction failed");
                        }
                    }
                } else {
                    debug!("No path found to end node {}/{}", j + 1, end_nodes.len());
                }
            }
        } // distances and reconstructor are dropped here
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
pub fn extract_path_sequences(
    graph: &DiGraph<NodeData, f64>,
    path: &[NodeIndex],
) -> Vec<String> {
    path.iter()
        .map(|&idx| graph[idx].sequence.clone())
        .collect()
}
//#####################
#[derive(Debug)]
pub enum AssemblyHistory {
    Full {
        path: Vec<String>,
        total_weight: f64,
        mean_coverage: f64,
    },
    Minimal {
        node_count: usize,
        total_weight: f64,
        mean_coverage: f64,
    }
}

impl Drop for AssemblyHistory {
    fn drop(&mut self) {
        if let AssemblyHistory::Full { path, .. } = self {
            path.clear(); // Explicitly clear the vector before dropping
        }
    }
}

#[derive(Debug)]
pub struct PathFindingResult {
    pub assembled_sequence: String,
    pub history: AssemblyHistory,
}

impl PathFindingResult {
    // Constructor for full history
    pub fn with_full_history(
        path: Vec<String>,
        total_weight: f64,
        mean_coverage: f64,
        assembled_sequence: String,
    ) -> Self {
        Self {
            assembled_sequence,
            history: AssemblyHistory::Full {
                path,
                total_weight,
                mean_coverage,
            },
        }
    }

    // Constructor for minimal history
    pub fn with_minimal_history(
        node_count: usize,
        total_weight: f64,
        mean_coverage: f64,
        assembled_sequence: String,
    ) -> Self {
        Self {
            assembled_sequence,
            history: AssemblyHistory::Minimal {
                node_count,
                total_weight,
                mean_coverage,
            },
        }
    }

    // Helper method to get statistics regardless of history type
    pub fn get_stats(&self) -> (f64, f64) {
        match &self.history {
            AssemblyHistory::Full { total_weight, mean_coverage, .. } => {
                (*total_weight, *mean_coverage)
            },
            AssemblyHistory::Minimal { total_weight, mean_coverage, .. } => {
                (*total_weight, *mean_coverage)
            }
        }
    }

    // Helper method to get path length
    pub fn path_length(&self) -> usize {
        match &self.history {
            AssemblyHistory::Full { path, .. } => path.len(),
            AssemblyHistory::Minimal { node_count, .. } => *node_count,
        }
    }

    // Optional: Method to convert full history to minimal
    pub fn minimize_history(&mut self) {
        if let AssemblyHistory::Full {ref path, total_weight, mean_coverage } = std::mem::replace(&mut self.history, AssemblyHistory::Minimal {
            node_count: 0,
            total_weight: 0.0,
            mean_coverage: 0.0,
        }) {
            self.history = AssemblyHistory::Minimal {
                node_count: path.len(),
                total_weight,
                mean_coverage,
            };
        }
    }
}

pub fn assemble_with_path_finding<K: Kmer + Send + Sync + Debug + 'static>(
    preliminary_graph: &DebruijnGraph<K, u16>,
    start_anchor: &str,
    end_anchor: &str,
        keep_full_history: Option<bool>,  // Make it optional

) -> Result<PathFindingResult, String> {
    let _ = env_logger::try_init();
    info!("Converting de Bruijn graph to weighted digraph for path finding");
    
    let keep_history = keep_full_history.unwrap_or(true);

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

        let result = if keep_full_history.expect("expected a boolean") {
            PathFindingResult::with_full_history(
                sequences,
                total_weight,
                mean_coverage,
                assembled_sequence,
            )
        } else {
            PathFindingResult::with_minimal_history(
                sequences.len(),
                total_weight,
                mean_coverage,
                assembled_sequence,
            )
        };

        Ok(result)
    } else {
        Err("No valid path found between anchors".to_string())
    }
}
