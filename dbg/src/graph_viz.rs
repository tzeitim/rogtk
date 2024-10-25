use std::fs::File;
use std::io::Write;
use std::path::Path;

use debruijn::*;
use debruijn::graph::DebruijnGraph;
use anyhow::Result;

pub fn export_graph<K: Kmer, D: std::fmt::Debug>(
    graph: &DebruijnGraph<K, D>,
    path: &str,
    title: &str
) -> Result<()> {
    let mut file = File::create(Path::new(path))?;

    // Calculate graph stats for title
    let mut terminal_count = 0;
    let mut isolated_count = 0;
    for node_id in 0..graph.len() {
        let node = graph.get_node(node_id);
        let left_edges = node.l_edges();
        let right_edges = node.r_edges();
        
        if left_edges.is_empty() || right_edges.is_empty() {
            terminal_count += 1;
        }
        if left_edges.is_empty() && right_edges.is_empty() {
            isolated_count += 1;
        }
    }

    // Write DOT header with stats
    writeln!(file, "digraph {} {{", title)?;
    writeln!(file, "    label=\"{} de Bruijn Graph\\n", title)?;
    writeln!(file, "Nodes: {}  Terminal: {}  Isolated: {}\"", 
             graph.len(), terminal_count, isolated_count)?;
    writeln!(file, "    labelloc=\"t\"")?;
    writeln!(file, "    node [shape=box]")?;

    // Write nodes with standard Graphviz colors
    for node_id in 0..graph.len() {
        let node = graph.get_node(node_id);
        let seq = node.sequence().to_string();
        let data = node.data();
        
        let left_edges = node.l_edges();
        let right_edges = node.r_edges();
        
        let color = if left_edges.is_empty() && right_edges.is_empty() {
            "#ff110030"  // Isolated nodes
        } else if left_edges.is_empty() || right_edges.is_empty() {
            "#ff110030"  // Terminal nodes
        } else {
            "#4895fa30"  // Internal nodes
        };
        
        writeln!(file, "    n{} [label=\"ID: {}\\nSeq: {}\\nData: {:?}\" style=filled fillcolor=\"{}\"]", 
            node_id, node_id, seq, data, color)?;
    }

    // Write edges
    for from_id in 0..graph.len() {
        let node = graph.get_node(from_id);
        
        // Add right edges
        for (to_id, incoming_dir, _) in node.r_edges() {
            writeln!(file, "    n{} -> n{} [label=\"{:?}\"]",
                from_id, to_id, incoming_dir)?;
        }
    }

    writeln!(file, "}}")?;
    Ok(())
}
