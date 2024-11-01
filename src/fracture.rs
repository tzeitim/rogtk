use pyo3::prelude::*;
use std::path::Path;

use bio::io::fasta;
use boomphf::hashmap::BoomHashMap2;
use debruijn::*;
use debruijn::dna_string::*;
use debruijn::filter::*;
use debruijn::compression::*;
use debruijn::kmer::*;
use debruijn::graph::BaseGraph;
use log::*;
use env_logger;
use anyhow::Result;
use std::fmt::Debug;

use crate::graph_viz::export_graph;

fn estimate_k(sequences: &[String]) -> usize {
    // Safety checks
    if sequences.is_empty() {
        return 31; // Default k value if no sequences
    }
    
    // Calculate mean read length safely
    let mut total_length: usize = 0;
    let mut count = 0;
    
    for seq in sequences {
        if !seq.is_empty() {
            total_length = total_length.saturating_add(seq.len());
            count += 1;
        }
    }
    
    if count == 0 {
        return 31; // Default k value if all sequences are empty
    }
    
    // Calculate mean length safely using floating point
    let mean_length = (total_length as f64) / (count as f64);
    
    // Estimate k as mean_length/3, clamped to reasonable values
    let k = (mean_length / 3.0).round() as usize;
    
    // Ensure k is odd and within valid range (min 11, max 63)
    let k = if k % 2 == 0 { k - 1 } else { k };
    k.clamp(11, 63)
}

#[allow(dead_code)]
trait DbgInterface {
    fn node_count(&self) -> usize;
    fn terminal_count(&self) -> usize;
    fn isolated_count(&self) -> usize;
    fn compress_and_get_contigs(&self, min_size: usize) -> Vec<String>;
}

// Implement the trait for our DbgResult
impl<K: Kmer + Send + Sync + Debug + 'static> DbgInterface for DbgResult<K> {
    fn node_count(&self) -> usize { self.node_count }
    fn terminal_count(&self) -> usize { self.terminal_count }
    fn isolated_count(&self) -> usize { self.isolated_count }
    
    fn compress_and_get_contigs(&self, min_size: usize) -> Vec<String> {
        let spec = SimpleCompress::new(|d1: u16, d2: &u16| {
            debug!("Merging counts: d1={}, d2={}", d1, *d2);
            d1.saturating_add(*d2)
        });

        let compressed_graph = compress_kmers_with_hash(
            true, 
            &spec,
            &self.valid_kmers
        ).finish();

        debug!("Compressed graph has {} nodes", compressed_graph.len());

        let mut contigs = Vec::new();
        for (i, node) in compressed_graph.iter_nodes().enumerate() {
            let seq = node.sequence();
            if seq.len() >= min_size {
                let contig = seq.to_string();
                debug!("Found contig {}: length={}, sequence={}", i, contig.len(), contig);
                contigs.push(contig);
            }
        }
        contigs
    }
}

// Simplified struct without unused fields
struct DbgResult<K: Kmer> {
    node_count: usize,
    terminal_count: usize,
    isolated_count: usize,
    valid_kmers: BoomHashMap2<K, Exts, u16>
}

fn analyze_dbg<K: Kmer + Send + Sync + Debug + 'static>(
    seq_tuples: &[(DnaString, Exts, ())],
    min_coverage: usize
) -> Option<DbgResult<K>> {
    info!("Using kmer size k={}", K::k());
    let (valid_kmers, _) = filter_kmers::<K, _, _, _, _>(
        seq_tuples,
        &Box::new(CountFilter::new(min_coverage)),
        true,
        true,
        4
    );

    // Count important graph features
    let mut terminal_count = 0;
    let mut isolated_count = 0;
    let mut coverage_dist = std::collections::HashMap::new();

    for (kmer, exts, count) in valid_kmers.iter() {
        debug!("\nKmer: {}", kmer.to_string());
        debug!("Coverage: {}", count);
        debug!("Extensions: {:?}", exts);
        
        let left_exts = exts.get(Dir::Left);
        let right_exts = exts.get(Dir::Right);
        
        debug!("Left extensions: {:?}", left_exts);
        debug!("Right extensions: {:?}", right_exts);

        // Track terminal and isolated nodes
        if left_exts.is_empty() || right_exts.is_empty() {
            terminal_count += 1;
        }
        if left_exts.is_empty() && right_exts.is_empty() {
            isolated_count += 1;
        }

        // Track coverage distribution
        *coverage_dist.entry(*count).or_insert(0) += 1;
    }

    debug!("\nCoverage distribution:");
    for (cov, count) in coverage_dist.iter() {
        debug!("Coverage {}: {} kmers", cov, count);
    }

    Some(DbgResult {
        node_count: valid_kmers.len(),
        terminal_count,
        isolated_count,
        valid_kmers,
    })
}

/// Read sequences from a FASTA file
fn read_fasta_sequences(fasta_path: &Path) -> Result<Vec<DnaString>> {
    let reader = fasta::Reader::from_file(fasta_path)?;
    let sequences: Vec<DnaString> = reader
        .records()
        .filter_map(|record| {
            let record = match record {
                Ok(rec) => rec,
                Err(e) => {
                    warn!("Error reading FASTA record: {}", e);
                    return None;
                }
            };
            
            let seq = record.seq();
            // Convert to uppercase and validate
            let seq_str = String::from_utf8_lossy(seq).to_uppercase();
            if seq_str.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                Some(DnaString::from_dna_string(&seq_str))
            } else {
                warn!("Skipping sequence with invalid characters: {}", record.id());
                None
            }
        })
        .collect();

    Ok(sequences)
}
pub fn assemble_sequences(
    sequences: Vec<String>, 
    k: usize, 
    min_coverage: usize, 
    export_graphs: Option<bool>,
    only_largest: Option<bool>,
    min_length: Option<usize>,
    auto_k: Option<bool>,
    prefix: Option<String>,
) -> Result<Vec<String>> {
    info!("Starting assembly of {} sequences with k={}, min_coverage={}", sequences.len(), k, min_coverage);
    
    let k = if auto_k.unwrap_or(false) {
        let estimated_k = estimate_k(&sequences);
        info!("Automatically estimated k = {} based on mean read length", estimated_k);
        estimated_k
    } else {
        k
    };


    if k > 64 {
        error!("K-mer size {} not supported (maximum is 64)", k);
        return Ok(Vec::new());
    }

    // Convert sequences to DnaString
    let sequences: Vec<DnaString> = sequences.into_iter()
        .filter_map(|seq| {
            // Convert to uppercase and validate
            let seq = seq.to_uppercase();
            if seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                Some(DnaString::from_dna_string(&seq))
            } else {
                warn!("Skipping sequence with invalid characters");
                None
            }
        })
        .collect();

    debug!("Processed {} valid sequences", sequences.len());
    
    if sequences.is_empty() {
        warn!("No valid sequences to process!");
        return Ok(Vec::new());
    }

    let seq_tuples: Vec<(DnaString, Exts, ())> = sequences.into_iter()
        .map(|s| (s, Exts::empty(), ()))
        .collect();

    info!("Building De Bruijn graph with k={}", k);
    
    let prefix = prefix.unwrap_or_else(|| "assembly".to_string());
    
    let contigs = match k {
        k if k <= 4 => assemble_with_k::<Kmer4>(&seq_tuples, min_coverage, &prefix, export_graphs),
        k if k <= 8 => assemble_with_k::<Kmer8>(&seq_tuples, min_coverage, &prefix, export_graphs),
        k if k <= 16 => assemble_with_k::<Kmer16>(&seq_tuples, min_coverage, &prefix, export_graphs),
        k if k <= 32 => assemble_with_k::<Kmer32>(&seq_tuples, min_coverage, &prefix, export_graphs), 
        k if k <= 64 => assemble_with_k::<Kmer64>(&seq_tuples, min_coverage, &prefix, export_graphs),
        _ => {
            error!("K-mer size {} not supported. Please use k <= 64", k);
            Ok(Vec::new())
        }
    }?;

    // Filter by minimum length if specified
    let min_length = min_length.unwrap_or(0);
    let filtered_contigs: Vec<String> = contigs.into_iter()
        .filter(|contig| contig.len() >= min_length)
        .collect();

    if filtered_contigs.is_empty() {
        warn!("No contigs found with length >= {}", min_length);
        return Ok(Vec::new());
    }

    // Return only the largest contig if requested
    if only_largest.unwrap_or(false) {
        if let Some(largest) = filtered_contigs.into_iter()
            .max_by_key(|contig| contig.len()) {
            Ok(vec![largest])
        } else {
            Ok(Vec::new())
        }
    } else {
        Ok(filtered_contigs)
    }
}

pub fn assemble_fasta(fasta_path: &Path, k: usize, min_coverage: usize, export_graphs: Option<bool>) -> Result<Vec<String>> {
    info!("Starting assembly of {} with k={}, min_coverage={}", 
          fasta_path.display(), k, min_coverage);
    
    if k > 64 {
        error!("K-mer size {} not supported (maximum is 64)", k);
        return Ok(Vec::new());
    }

    let sequences = read_fasta_sequences(fasta_path)?;
    debug!("Loaded {} valid sequences from FASTA", sequences.len());
    
    if sequences.is_empty() {
        warn!("No valid sequences to process!");
        return Ok(Vec::new());
    }

    let seq_tuples: Vec<(DnaString, Exts, ())> = sequences.into_iter()
        .map(|s| (s, Exts::empty(), ()))
        .collect();

    info!("Building De Bruijn graph with k={}", k);
    
    // Generate output prefix from input filename
    let prefix = fasta_path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("assembly");
    
    match k {
        k if k <= 4 => assemble_with_k::<Kmer4>(&seq_tuples, min_coverage, prefix, export_graphs),
        k if k <= 8 => assemble_with_k::<Kmer8>(&seq_tuples, min_coverage, prefix, export_graphs),
        k if k <= 16 => assemble_with_k::<Kmer16>(&seq_tuples, min_coverage, prefix, export_graphs),
        k if k <= 32 => assemble_with_k::<Kmer32>(&seq_tuples, min_coverage, prefix, export_graphs), 
        k if k <= 64 => assemble_with_k::<Kmer64>(&seq_tuples, min_coverage, prefix, export_graphs),
        _ => {
            error!("K-mer size {} not supported. Please use k <= 64", k);
            Ok(Vec::new())
        }
    }
}

fn assemble_with_k<K: Kmer + Send + Sync + Debug + 'static>(
    seq_tuples: &[(DnaString, Exts, ())],
    min_coverage: usize,
    prefix: &str,
    export_graphs: Option<bool>,
) -> Result<Vec<String>> {
    let should_export = export_graphs.unwrap_or(true);
    
    // Get initial kmers and stats
    let preliminary_stats = analyze_dbg::<K>(seq_tuples, min_coverage)
        .ok_or_else(|| anyhow::anyhow!("Failed to create graph"))?;

    info!("Initial k-mer statistics:");
    info!("  Total kmers: {}", preliminary_stats.node_count);
    info!("  Terminal kmers: {}", preliminary_stats.terminal_count);
    info!("  Isolated kmers: {}", preliminary_stats.isolated_count);

    // Build initial uncompressed graph directly from kmers
    let mut preliminary_graph: BaseGraph<K, u16> = BaseGraph::new(true);
    for (kmer, exts, count) in preliminary_stats.valid_kmers.iter() {
        let kmer_bytes: Vec<u8> = (0..K::k()).map(|i| kmer.get(i)).collect();
        preliminary_graph.add(&kmer_bytes, *exts, *count);
    }
    let preliminary_graph = preliminary_graph.finish();

    // Export uncompressed preliminary graph if enabled
    if should_export {
        let prelim_path = format!("{}_preliminary.dot", prefix);
        export_graph(&preliminary_graph, &prelim_path, "Preliminary ")?;
        info!("Exported preliminary graph to {}", prelim_path);
    }

    // Build compressed graph
    let compressed_graph = {
        let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1.saturating_add(*d2));
        compress_graph(true, &spec, preliminary_graph, None)
    };

    // Export compressed graph if enabled
    if should_export {
        let comp_path = format!("{}_compressed.dot", prefix);
        export_graph(&compressed_graph, &comp_path, "Compressed")?;
        info!("Exported compressed graph to {}", comp_path);
    }

    // Extract contigs from compressed graph
    let mut contigs = Vec::new();
    for node in compressed_graph.iter_nodes() {
        let seq = node.sequence();
        if seq.len() >= K::k() {
            contigs.push(seq.to_string());
        }
    }

    info!("Assembly complete. Found {} contigs", contigs.len());
    Ok(contigs)
}

#[pyfunction]
#[pyo3(signature = (fasta_path, k, min_coverage, min_length=200, export_graphs=None))]
pub fn fracture_fasta(
    fasta_path: String, 
    k: usize, 
    min_coverage: usize,
    min_length: Option<usize>,
    export_graphs: Option<bool>
) -> PyResult<String> {
    // Initialize logger
    //env_logger::init();
    let _ = env_logger::try_init();
    
    // Use provided min_length or default to 200
    let min_length = min_length.unwrap_or(200);
    
    // Perform assembly
    let contigs = assemble_fasta(Path::new(&fasta_path), k, min_coverage, export_graphs)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    
    // Find the largest contig
    let largest_contig = contigs.into_iter()
        .filter(|contig| contig.len() >= min_length)
        .max_by_key(|contig| contig.len())
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("No contigs found with length >= {}", min_length)
        ))?
        .to_string();
    
    Ok(largest_contig)
}

#[pyfunction]
#[pyo3(signature = (sequences, k, min_coverage, min_length=200, export_graphs=None, only_largest=None, auto_k=None, prefix=None))]
/// Assemble sequences using a de Bruijn graph approach
/// 
/// Args:
///     sequences (List[str]): List of DNA sequences to assemble
///     k (int): K-mer size for graph construction (used if auto_k=False)
///     min_coverage (int): Minimum k-mer coverage threshold
///     min_length (int, optional): Minimum contig length to return. Defaults to 200.
///     export_graphs (bool, optional): Whether to export graph visualization files. Defaults to None.
///     only_largest (bool, optional): Return only the largest contig. Defaults to None.
///     auto_k (bool, optional): Automatically estimate optimal k-mer size. Defaults to True.
///     prefix (str, optional): Prefix for output files. Defaults to "assembly".
/// 
/// Returns:
///     str: Assembled contigs separated by newlines, or single contig if only_largest=True

pub fn fracture_sequences(
    sequences: Vec<String>, 
    k: usize, 
    min_coverage: usize,
    min_length: Option<usize>,
    export_graphs: Option<bool>,
    only_largest: Option<bool>,
    auto_k: Option<bool>,
    prefix: Option<String>,
) -> PyResult<String> {
    // Try to initialize logger, ignore if already initialized
    let _ = env_logger::try_init();
    
    // Perform assembly
    let contigs = assemble_sequences(
        sequences, 
        k, 
        min_coverage, 
        export_graphs,
        only_largest,
        min_length,
        auto_k,
        prefix
    ).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    if contigs.is_empty() {
        return Ok(String::new());
    }

    if only_largest.unwrap_or(false) {
        Ok(contigs[0].clone())
    } else {
        Ok(contigs.join("\n"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_fasta() -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        let reads = vec![
            ">read1",
            "ATGCATGCATGCTAGCTGATCGATCGTAGCTAGCTAGCTGATCGATCGTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGT",
            ">read2", 
            "GTAGCTAGCTAGCTGATCGATCGTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGT",
            ">read3",
            "GTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGTACGTAGCTACGTACGTACGTAG",
            ">read4",
            "TACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCT"
        ];

        for line in reads {
            writeln!(file, "{}", line).unwrap();
        }
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_fasta_assembly() {
        let _ = env_logger::try_init();
        let test_file = create_test_fasta();
        
        // Use k=20 for testing instead of k=47
        let contigs = assemble_fasta(test_file.path(), 20, 1).unwrap();
        
        assert!(!contigs.is_empty(), "Should produce at least one contig");
        
        if let Some(contig) = contigs.first() {
            println!("First contig: {}", contig);
            assert!(contig.len() > 150, "Should produce a contig longer than 150bp");
        }
    }
}
