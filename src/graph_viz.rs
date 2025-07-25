use std::fs::File;
use std::io::Write;
use std::path::Path;

use debruijn::*;
use debruijn::graph::DebruijnGraph;
use anyhow::Result;
use polars::prelude::*;
use log::info;

/// Convert graph into a single DataFrame where each row represents a node
/// and its associated edges
fn graph_to_dataframe<K: Kmer, D: std::fmt::Debug>(
    graph: &DebruijnGraph<K, D>
) -> Result<DataFrame> {
    let mut node_ids = Vec::new();
    let mut sequences = Vec::new();
    let mut node_data = Vec::new();
    let mut node_types = Vec::new();
    let mut outgoing_edges = Vec::new();
    let mut outgoing_directions = Vec::new();
    
    // Single iteration to collect all node and edge data
    for node_id in 0..graph.len() {
        let node = graph.get_node(node_id);
        let left_edges = node.l_edges();
        let right_edges = node.r_edges();
        
        info!(" in df {}", node.sequence().to_string());
        // Basic node information
        node_ids.push(node_id as i64);
        sequences.push(node.sequence().to_string());
        node_data.push(format!("{:?}", node.data()));
        
        // Determine node type based on edge connectivity
        let node_type = if left_edges.is_empty() && right_edges.is_empty() {
            "isolated"
        } else if left_edges.is_empty() || right_edges.is_empty() {
            "terminal"
        } else {
            "internal"
        };
        node_types.push(node_type);
        
        // Collect edge information as strings
        let mut edge_targets = Vec::new();
        let mut edge_dirs = Vec::new();
        
        for (to_id, incoming_dir, _) in node.r_edges() {
            edge_targets.push(to_id.to_string());
            edge_dirs.push(format!("{:?}", incoming_dir));
        }
        
        // Store edge information as comma-separated strings
        outgoing_edges.push(edge_targets.join(","));
        outgoing_directions.push(edge_dirs.join(","));
    }
    
    // Create DataFrame with Series
    let df = DataFrame::new(vec![
        Series::new("node_id".into(), node_ids),
        Series::new("sequence".into(), sequences),
        Series::new("node_type".into(), node_types),
        Series::new("coverage".into(), node_data),
        Series::new("outgoing_nodes".into(), outgoing_edges),
        Series::new("outgoing_directions".into(), outgoing_directions),
    ])?;
    
    Ok(df)
}

/// Export graph to both DOT and CSV formats
pub fn export_graph<K: Kmer, D: std::fmt::Debug>(
    graph: &DebruijnGraph<K, D>,
    path: &str,
    title: &str
) -> Result<()> {
    // Convert graph to DataFrame
    let mut graph_df = graph_to_dataframe(graph)?;
    
    // Calculate graph statistics directly
    let total_nodes = graph_df.height();
    let terminal_count = graph_df.column("node_type")?
        .str()?
        .into_iter()
        .filter(|x| *x == Some("terminal"))
        .count();
    let isolated_count = graph_df.column("node_type")?
        .str()?
        .into_iter()
        .filter(|x| *x == Some("isolated"))
        .count();
    
    // Export CSV file
    let csv_path = Path::new(path).with_extension("csv");
    CsvWriter::new(File::create(csv_path)?)
        .finish(&mut graph_df)?;
    
    // Export DOT visualization
    export_dot_from_dataframe(
        &graph_df,
        path,
        title,
        total_nodes,
        terminal_count,
        isolated_count
    )?;
    
    Ok(())
}


/// Generate DOT file from DataFrame
fn export_dot_from_dataframe(
    graph_df: &DataFrame,
    path: &str,
    title: &str,
    total_nodes: usize,
    terminal_nodes: usize,
    isolated_nodes: usize,
) -> Result<()> {
    let mut file = File::create(path)?;
    
    // Write DOT header with stats
    writeln!(file, "digraph {} {{", title)?;
    writeln!(file, "    label=\"{} de Bruijn Graph\\n", title)?;
    writeln!(file, "Nodes: {}  Terminal: {}  Isolated: {}\"",
             total_nodes, terminal_nodes, isolated_nodes)?;
    writeln!(file, "    labelloc=\"t\"")?;
    writeln!(file, "    node [shape=box]")?;
    
    // Write nodes and edges
    let node_ids = graph_df.column("node_id")?.i64()?;
    let sequences = graph_df.column("sequence")?.str()?;
    let node_types = graph_df.column("node_type")?.str()?;
    let cov = graph_df.column("coverage")?.str()?;
    let outgoing_nodes = graph_df.column("outgoing_nodes")?.str()?;
    let outgoing_directions = graph_df.column("outgoing_directions")?.str()?;
    
    for idx in 0..graph_df.height() {
        let node_id = node_ids.get(idx).unwrap();
        let sequence = sequences.get(idx).unwrap_or("");
        let node_type = node_types.get(idx).unwrap_or("");
        let node_data = cov.get(idx).unwrap_or("");
        let edges = outgoing_nodes.get(idx).unwrap_or("");
        let directions = outgoing_directions.get(idx).unwrap_or("");
        
        let color = match node_type {
            "isolated" | "terminal" => "#ff110030",
            _ => "#4895fa30"
        };
        
        writeln!(file,
            "    n{} [label=\"ID: {}\\nSeq: {}\\ncov: {}\", style=filled, fillcolor=\"{}\"]",
            node_id, node_id, sequence, node_data, color
        )?;
        
        // Write edges if present
        if !edges.is_empty() {
            let targets: Vec<&str> = edges.split(',').collect();
            let edge_dirs: Vec<&str> = directions.split(',').collect();
            
            for (target, direction) in targets.iter().zip(edge_dirs.iter()) {
                writeln!(file, "    n{} -> n{} [label=\"{}\"]",
                    node_id, target, direction)?;
            }
        }
    }
    
    writeln!(file, "}}")?;
    Ok(())
}
