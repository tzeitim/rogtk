
use std::collections::HashMap;

#[derive(Debug)]
pub struct ComplexityScore {
    pub shannon_entropy: f64,
    pub linguistic_complexity: f64,
    pub homopolymer_fraction: f64,
    pub dinucleotide_entropy: f64,
    pub combined_score: f64,
}


impl ComplexityScore {
    pub fn is_low_complexity(&self, threshold: f64) -> bool {
        self.combined_score < threshold
    }
}

pub fn calculate_umi_complexity(umi: &str) -> ComplexityScore {
    let shannon = shannon_entropy(umi);
    let linguistic = linguistic_complexity(umi);
    let homopolymer = homopolymer_fraction(umi);
    let dinuc = dinucleotide_entropy(umi);
    
    // Combine metrics (weights can be adjusted)
    let combined = 0.3 * shannon 
                 + 0.3 * linguistic 
                 + 0.2 * (1.0 - homopolymer)  // Invert so higher is better
                 + 0.2 * dinuc;
    
    ComplexityScore {
        shannon_entropy: shannon,
        linguistic_complexity: linguistic,
        homopolymer_fraction: homopolymer,
        dinucleotide_entropy: dinuc,
        combined_score: combined,
    }
}

fn shannon_entropy(seq: &str) -> f64 {
    let mut counts = [0u32; 4];
    let mut total = 0u32;
    
    for byte in seq.bytes() {
        total += 1;
        match byte {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {} // Handle invalid characters if needed
        }
    }
    
    if total == 0 {
        return 0.0;
    }
    
    let mut entropy = 0.0;
    for &count in &counts {
        if count > 0 {
            let p = count as f64 / total as f64;
            entropy -= p * p.log2();
        }
    }
    
    entropy
}


// Linguistic complexity: ratio of unique k-mers to total possible k-mers
fn linguistic_complexity(seq: &str) -> f64 {
    if seq.len() < 3 {
        return 0.0;
    }
    
    let k = 3; // tri-nucleotides
    let mut kmers = HashMap::new();
    
    for window in seq.as_bytes().windows(k) {
        *kmers.entry(window).or_insert(0) += 1;
    }
    
    let unique_kmers = kmers.len() as f64;
    let max_possible = (seq.len() - k + 1).min(4_usize.pow(k as u32)) as f64;
    
    unique_kmers / max_possible
}

// Fraction of sequence that's part of homopolymer runs (3+ identical bases)
fn homopolymer_fraction(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    
    let bytes = seq.as_bytes();
    let mut in_homopolymer = 0;
    let mut i = 0;
    
    while i < bytes.len() {
        let current = bytes[i];
        let mut run_length = 1;
        
        while i + run_length < bytes.len() && bytes[i + run_length] == current {
            run_length += 1;
        }
        
        if run_length >= 3 {
            in_homopolymer += run_length;
        }
        
        i += run_length;
    }
    
    in_homopolymer as f64 / seq.len() as f64
}

// Entropy of dinucleotide frequencies
fn dinucleotide_entropy(seq: &str) -> f64 {
    if seq.len() < 2 {
        return 0.0;
    }
    
    let mut counts = HashMap::new();
    let bytes = seq.as_bytes();
    
    for window in bytes.windows(2) {
        *counts.entry((window[0], window[1])).or_insert(0) += 1;
    }
    
    let total = (seq.len() - 1) as f64;
    let mut entropy = 0.0;
    
    for &count in counts.values() {
        let p = count as f64 / total;
        entropy -= p * p.log2();
    }
    
    // Normalize to [0, 1] range (max possible is log2(16) = 4)
    entropy / 4.0
}

// Additional useful metric: longest homopolymer run
fn longest_homopolymer_run(seq: &str) -> usize {
    if seq.is_empty() {
        return 0;
    }
    
    let bytes = seq.as_bytes();
    let mut max_run = 1;
    let mut current_run = 1;
    
    for i in 1..bytes.len() {
        if bytes[i] == bytes[i-1] {
            current_run += 1;
            max_run = max_run.max(current_run);
        } else {
            current_run = 1;
        }
    }
    
    max_run
}

// Dust score - another common low-complexity measure
fn dust_score(seq: &str, window_size: usize) -> f64 {
    if seq.len() < window_size {
        return 0.0;
    }
    
    let mut total_score = 0.0;
    let bytes = seq.as_bytes();
    
    for i in 0..=seq.len() - window_size {
        let window = &bytes[i..i + window_size];
        let mut triplet_counts = HashMap::new();
        
        for j in 0..=window_size - 3 {
            let triplet = &window[j..j + 3];
            *triplet_counts.entry(triplet).or_insert(0) += 1;
        }
        
        let mut window_score = 0.0;
        for &count in triplet_counts.values() {
            if count > 1 {
                window_score += (count * (count - 1)) as f64 / 2.0;
            }
        }
        
        total_score += window_score;
    }
    
    // Normalize
    total_score / (seq.len() - window_size + 1) as f64
}
