
use bio::io::fasta;
use boomphf::hashmap::BoomHashMap2;
use debruijn::*;
use debruijn::dna_string::*;
use debruijn::filter::*;
use debruijn::compression::*;
use debruijn::kmer::*;
use log::*;
use env_logger;
use std::path::Path;
use anyhow::Result;
use std::fmt::Debug;

// Simplified trait with only needed methods
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

pub fn assemble_fasta(fasta_path: &Path, k: usize, min_coverage: usize) -> Result<Vec<String>> {
    info!("Starting assembly of {} with k={}, min_coverage={}", 
          fasta_path.display(), k, min_coverage);
    
    if k > 64 {
        error!("K-mer size {} not supported (maximum is 64)", k);
        return Ok(Vec::new());
    }

    // 1. Read and validate sequences from FASTA
    let sequences = read_fasta_sequences(fasta_path)?;
    debug!("Loaded {} valid sequences from FASTA", sequences.len());
    
    // Print each sequence for debugging
    for (i, seq) in sequences.iter().enumerate() {
        debug!("Sequence {}: length={}, content={:?}", i, seq.len(), seq);
    }
    
    if sequences.is_empty() {
        warn!("No valid sequences to process!");
        return Ok(Vec::new());
    }

    // 2. Create sequence tuples with empty extensions
    let seq_tuples: Vec<(DnaString, Exts, ())> = sequences.into_iter()
        .map(|s| (s, Exts::empty(), ()))
        .collect();

    // 3. Filter kmers and build initial DBG
    info!("Building De Bruijn graph with k={}", k);
    
    let dbg_stats: Option<Box<dyn DbgInterface>> = match k {
        k if k <= 4 => analyze_dbg::<Kmer4>(&seq_tuples, min_coverage)
            .map(|stats| Box::new(stats) as Box<dyn DbgInterface>),
        k if k <= 8 => analyze_dbg::<Kmer8>(&seq_tuples, min_coverage)
            .map(|stats| Box::new(stats) as Box<dyn DbgInterface>),
        k if k <= 16 => analyze_dbg::<Kmer16>(&seq_tuples, min_coverage)
            .map(|stats| Box::new(stats) as Box<dyn DbgInterface>),
        k if k <= 32 => analyze_dbg::<Kmer32>(&seq_tuples, min_coverage)
            .map(|stats| Box::new(stats) as Box<dyn DbgInterface>),
        k if k <= 64 => analyze_dbg::<Kmer64>(&seq_tuples, min_coverage)
            .map(|stats| Box::new(stats) as Box<dyn DbgInterface>),
        _ => {
            error!("K-mer size {} not supported. Please use k <= 64", k);
            None
        }
    };

    let dbg_stats = match dbg_stats {
        Some(stats) => stats,
        None => {
            error!("Failed to create DBG with k={}", k);
            return Ok(Vec::new());
        }
    };

    info!("DBG Statistics:");
    info!("Total nodes: {}", dbg_stats.node_count());
    info!("Terminal nodes: {}", dbg_stats.terminal_count());
    info!("Isolated nodes: {}", dbg_stats.isolated_count());

    // Get contigs using the trait method
    let contigs = dbg_stats.compress_and_get_contigs(k);

    info!("Assembly complete. Found {} contigs", contigs.len());
    Ok(contigs)
}

fn main() -> Result<()> {
    // Initialize logger
    env_logger::init();

    let fasta_path = std::env::args().nth(1).expect("Please provide a FASTA file path");
    let k: usize = std::env::args().nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20); // Default to k=20 instead of 47
    let min_coverage: usize = std::env::args().nth(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);

    let contigs = assemble_fasta(Path::new(&fasta_path), k, min_coverage)?;
    
    println!("\nAssembled {} contigs:", contigs.len());
    for (i, contig) in contigs.iter().enumerate() {
        println!("Contig {}: {} (length={})", i+1, contig, contig.len());
    }

    Ok(())
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
