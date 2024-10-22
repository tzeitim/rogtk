use debruijn::*;
use debruijn::dna_string::*;
use debruijn::filter::*;
use debruijn::compression::*;
use debruijn::kmer::*;
use log::*;
use env_logger;

pub fn assemble_reads(reads: Vec<&str>, k: usize, min_coverage: usize) -> Vec<String> {
    info!("Starting assembly with k={}, min_coverage={}", k, min_coverage);
    
    if k > 64 {
        error!("K-mer size {} not supported (maximum is 64)", k);
        return Vec::new();
    }

    // 1. Convert reads to DnaString format and validate
    let sequences: Vec<DnaString> = reads.iter()
        .filter_map(|r| {
            // Convert to uppercase first
            let r_upper = r.to_uppercase();
            if r_upper.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                Some(DnaString::from_dna_string(&r_upper))
            } else {
                warn!("Skipping read with invalid characters: {}", r);
                None
            }
        })
        .collect();
    
    debug!("Converted {} valid reads to DnaString format", sequences.len());
    
    // Print each sequence and its k-mers for debugging
    for (i, seq) in sequences.iter().enumerate() {
        debug!("Sequence {}: length={}, content={:?}", i, seq.len(), seq);
        
        // Debug: Show first k-mer from each sequence
        if seq.len() >= k {
            debug!("First k-mer from sequence {}: {:?}", i, seq.slice(0, k));
        }
    }
    
    if sequences.is_empty() {
        warn!("No valid sequences to process!");
        return Vec::new();
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
        return Vec::new();
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
    contigs
}

fn main() {
    // Initialize logger
    env_logger::init();

    // Create overlapping 120bp reads that should assemble into a longer contig
    // Each read overlaps with the next by 90bp to ensure good coverage
    let reads = vec![
        "ATGCATGCATGCTAGCTGATCGATCGTAGCTAGCTAGCTGATCGATCGTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGT",
        "GTAGCTAGCTAGCTGATCGATCGTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGT",
        "GTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGTACGTAGCTACGTACGTACGTAG",
        "TACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCT"
    ];

    info!("Starting assembly of {} reads", reads.len());
    
    // Use k=47 for longer read assembly
    let contigs = assemble_reads(reads, 47, 1);
    
    println!("\nAssembled {} contigs:", contigs.len());
    for (i, contig) in contigs.iter().enumerate() {
        println!("Contig {}: {} (length={})", i+1, contig, contig.len());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_long_read_assembly() {
        // Initialize logger for test
        let _ = env_logger::try_init();
        
        let reads = vec![
            "ATGCATGCATGCTAGCTGATCGATCGTAGCTAGCTAGCTGATCGATCGTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGT",
            "GTAGCTAGCTAGCTGATCGATCGTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGT",
            "GTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGTACGTAGCTACGTACGTACGTAG",
            "TACGTACGTACGTAGCTAGCTGATCGTAGCTACGTAGCTAGCTAGCTGATCGTACGTACGTAGCTGATCGATCGTAGCTACGTACGTACGTAGCTACGTACGTACGTAGCTAGCTGATCGTAGCT"
        ];

        let contigs = assemble_reads(reads, 47, 1);
        assert!(!contigs.is_empty(), "Should produce at least one contig");
        
        if let Some(contig) = contigs.first() {
            println!("First contig: {}", contig);
            assert!(contig.len() > 150, "Should produce a contig longer than 150bp");
        }
    }
}
