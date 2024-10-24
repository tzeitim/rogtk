use dbg_uncompressed::*;
use debruijn::*;
use log::*;
use env_logger;

fn main() {
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(log::LevelFilter::Debug)
        .init();

    // Example DNA reads with overlaps
    let reads = vec![
        "ACGTACGT",
        "CGTACGTA", 
        "GTACGTAC",
        "TACGTACG"
    ];

    info!("Building uncompressed De Bruijn graph from {} reads", reads.len());
    for (i, read) in reads.iter().enumerate() {
        info!("Read {}: {} (length={})", i, read, read.len());
    }
    
    // Build graph with k=4 and minimum coverage of 1
    let k = 4;
    let min_coverage = 1;
    
    match build_uncompressed_dbg(&reads, k, min_coverage) {
        Some(dbg) => {
            println!("\nUncompressed De Bruijn Graph Statistics:");
            println!("Number of k-mers: {}", dbg.len());
            
            if dbg.len() > 0 {
                println!("\nDetailed k-mer information:");
                for i in 0..dbg.len() {
                    let kmer = dbg.get_kmer(i);
                    let exts = dbg.get_extensions(i);
                    let count = dbg.get_count(i);
                    
                    println!("\nK-mer node {}:", i);
                    println!("  Sequence: {}", kmer.to_string());
                    println!("  Coverage: {}", count);
                    
                    print!("  Left extensions: ");
                    for base in 0..4 {
                        if exts.has_ext(Dir::Left, base) {
                            print!("{} ", bits_to_base(base));
                        }
                    }
                    println!();
                    
                    print!("  Right extensions: ");
                    for base in 0..4 {
                        if exts.has_ext(Dir::Right, base) {
                            print!("{} ", bits_to_base(base));
                        }
                    }
                    println!();
                }
            }
        },
        None => {
            error!("Failed to build De Bruijn graph!");
        }
    }
}
