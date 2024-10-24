use bio::io::fasta;
use debruijn::*;
use debruijn::dna_string::*;
use debruijn::filter::*;
use debruijn::compression::*;
use debruijn::kmer::*;
use log::*;
use env_logger;
use std::path::Path;
use anyhow::Result;

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
    
    // Print each sequence and its k-mers for debugging
    for (i, seq) in sequences.iter().enumerate() {
        debug!("Sequence {}: length={}, content={:?}", i, seq.len(), seq);
        
        if seq.len() >= k {
            debug!("First k-mer from sequence {}: {:?}", i, seq.slice(0, k));
        }
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
    
    let (valid_kmers, all_kmers) = filter_kmers::<Kmer64, _, _, _, _>(
        &seq_tuples,
        &Box::new(CountFilter::new(min_coverage)),
        false,  // stranded
        true,   // report all kmers  
        4       // memory size
    );

    debug!("Found {} total k-mers before filtering", all_kmers.len());
    debug!("Found {} valid k-mers after filtering", valid_kmers.len());

    if valid_kmers.is_empty() {
        warn!("No valid k-mers found after filtering! Try reducing min_coverage");
        return Ok(Vec::new());
    }

    // 4. Compress the graph
    info!("Compressing graph...");
    let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1.saturating_add(*d2));
    let compressed_graph = compress_kmers_with_hash(
        false,
        &spec, 
        &valid_kmers
    ).finish();

    debug!("Compressed graph has {} nodes", compressed_graph.len());

    // 5. Extract contigs
    let mut contigs = Vec::new();
    for (i, node) in compressed_graph.iter_nodes().enumerate() {
        let seq = node.sequence();
        if seq.len() >= k {
            let contig = seq.to_string();
            debug!("Found contig {}: length={}, sequence={}", i, contig.len(), contig);
            contigs.push(contig);
        }
    }

    info!("Assembly complete. Found {} contigs", contigs.len());
    Ok(contigs)
}

fn main() -> Result<()> {
    // Initialize logger
    env_logger::init();

    let fasta_path = std::env::args().nth(1).expect("Please provide a FASTA file path");
    let k: usize = std::env::args().nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(47);
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
        let contigs = assemble_fasta(test_file.path(), 47, 1).unwrap();
        
        assert!(!contigs.is_empty(), "Should produce at least one contig");
        
        if let Some(contig) = contigs.first() {
            println!("First contig: {}", contig);
            assert!(contig.len() > 150, "Should produce a contig longer than 150bp");
        }
    }
}
