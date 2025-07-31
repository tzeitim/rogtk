# HTSlib Parallel Optimization Instructions

*Critical Header Resolution Fix and Additional Performance Optimizations*

## 🔍 **The Problem: Chromosome Name Resolution**

### **Current Broken Implementation**
In our optimized parallel architecture, look at this code in `process_htslib_records_to_batch()` (around line 1551):

```rust
// 2. Chromosome - use tid to resolve reference name
let chrom = if hts_record.tid() >= 0 {
    // We need the header to resolve chromosome names, but we can't access it here
    // For now, use the tid directly as string - this is a limitation we need to fix
    Some(format!("chr_{}", hts_record.tid()))  // ❌ WRONG!
} else {
    None
};
```

### **What's Actually Happening**
1. **HTSlib stores chromosomes as integers**: `tid` (target ID) like `0, 1, 2, 3...`
2. **We're creating fake names**: `"chr_0", "chr_1", "chr_2"` instead of real names
3. **Real names should be**: `"chr1", "chr2", "chrX", "chrMT"` etc.

### **Why This Is Really Bad**
```python
# What we're outputting now (WRONG):
record_1: chrom = "chr_0"    # Should be "chr1" 
record_2: chrom = "chr_1"    # Should be "chr2"
record_3: chrom = "chr_23"   # Should be "chrX"

# This makes the output data UNUSABLE for downstream analysis!
```

## 🏗️ **The Architecture Problem**

### **Why We Can't Access the Header**
Our parallel architecture looks like this:

```rust
Main Thread (has header) → [Queue] → Worker Threads (NO header access)
     ↓                                        ↓
Has HeaderView with                    Only gets Vec<hts_bam::Record>
chromosome mapping                     Can't resolve tid → name
```

**The Issue**: Worker threads process `hts_bam::Record` objects but don't have access to the `HeaderView` needed to convert `tid` integers to actual chromosome names.

### **Current vs Correct Workflow**

#### ❌ **Current Broken Flow**:
```rust
// Worker thread receives:
hts_record.tid() = 0
// But has no way to know that tid=0 means "chr1"
// So it creates: format!("chr_{}", 0) = "chr_0"  // WRONG!
```

#### ✅ **Correct Flow Should Be**:
```rust
// Worker thread should:
hts_record.tid() = 0
header.tid2name(0) = "chr1"  // But worker can't access header!
// Output: "chr1"  // CORRECT!
```

## 💡 **The Solution: Pre-Built Chromosome Lookup**

### **Step 1: Build Lookup Table in Main Thread**
```rust
// Add this helper function to bam.rs:
#[cfg(feature = "htslib")]
fn build_chromosome_lookup(header: &hts_bam::HeaderView) -> Arc<Vec<String>> {
    let mut lookup = Vec::new();
    for tid in 0..header.target_count() {
        let chrom_name = String::from_utf8_lossy(header.tid2name(tid as u32)).into_owned();
        lookup.push(chrom_name);
    }
    Arc::new(lookup)
}

// Example result:
// lookup[0] = "chr1"
// lookup[1] = "chr2" 
// lookup[23] = "chrX"
// lookup[24] = "chrY"
```

### **Step 2: Pass Lookup to Worker Threads**
```rust
// Modify process_htslib_records_to_batch to accept lookup:
#[cfg(feature = "htslib")]
fn process_htslib_records_to_batch(
    hts_records: &[hts_bam::Record],
    chrom_lookup: &Arc<Vec<String>>,  // ← ADD THIS
    include_sequence: bool,
    include_quality: bool,
    batch_id: usize,
) -> Result<RecordBatch, Box<dyn std::error::Error + Send + Sync>> {
```

### **Step 3: Use Lookup in Worker**
```rust
// Replace the broken chromosome resolution code with:
let chrom = if hts_record.tid() >= 0 {
    let tid = hts_record.tid() as usize;
    if tid < chrom_lookup.len() {
        Some(chrom_lookup[tid].clone())  // ✅ CORRECT!
    } else {
        None  // Invalid tid
    }
} else {
    None  // Unmapped read
};
```

## 🛠️ **Complete Implementation Steps**

### **1. Add Helper Function**
Add this function after the existing helper functions in `bam.rs`:

```rust
// Helper function to build chromosome lookup table
#[cfg(feature = "htslib")]
fn build_chromosome_lookup(header: &hts_bam::HeaderView) -> Arc<Vec<String>> {
    let mut lookup = Vec::new();
    for tid in 0..header.target_count() {
        let chrom_name = String::from_utf8_lossy(header.tid2name(tid as u32)).into_owned();
        lookup.push(chrom_name);
    }
    Arc::new(lookup)
}
```

### **2. Modify Function Signature**
Update `process_htslib_records_to_batch()` signature:

```rust
// OLD:
fn process_htslib_records_to_batch(
    hts_records: &[hts_bam::Record],
    include_sequence: bool,
    include_quality: bool, 
    batch_id: usize,
) -> Result<RecordBatch, Box<dyn std::error::Error + Send + Sync>>

// NEW:
fn process_htslib_records_to_batch(
    hts_records: &[hts_bam::Record],
    chrom_lookup: &Arc<Vec<String>>,  // ← Add this parameter
    include_sequence: bool,
    include_quality: bool,
    batch_id: usize,
) -> Result<RecordBatch, Box<dyn std::error::Error + Send + Sync>>
```

### **3. Build Lookup in Main Function**
In `bam_to_arrow_ipc_htslib_parallel()`, after creating the header:

```rust
// Find this existing code:
let header_view = hts_reader.header().clone();

// Add this line after it:
let chrom_lookup = build_chromosome_lookup(&header_view);
```

### **4. Pass Lookup to Workers**
Update the worker spawn code to pass the lookup:

```rust
// Find the existing worker spawn code and modify the process call:
let batch_result = process_htslib_records_to_batch(
    &hts_records,
    &chrom_lookup,  // ← Add this parameter
    include_sequence,
    include_quality,
    batch_id,
);
```

### **5. Update Chromosome Resolution Code**
In `process_htslib_records_to_batch()`, replace the broken chromosome code:

```rust
// REPLACE THIS (around line 1551):
let chrom = if hts_record.tid() >= 0 {
    // We need the header to resolve chromosome names, but we can't access it here
    // For now, use the tid directly as string - this is a limitation we need to fix
    Some(format!("chr_{}", hts_record.tid()))
} else {
    None
};

// WITH THIS:
let chrom = if hts_record.tid() >= 0 {
    let tid = hts_record.tid() as usize;
    if tid < chrom_lookup.len() {
        Some(chrom_lookup[tid].clone())
    } else {
        None  // Invalid tid
    }
} else {
    None  // Unmapped read
};
```

## 📊 **Expected Performance Impact**

### **Current Performance Cost**
```rust
// Every record does this (EXPENSIVE):
Some(format!("chr_{}", hts_record.tid()))  // Allocates new String every time
```

**For 1M records**: 1,000,000 × `format!()` calls = millions of string allocations

### **Optimized Performance**
```rust
// Every record does this (CHEAP):
Some(chrom_lookup[tid].clone())  // Just clones existing String
```

**For 1M records**: 1,000,000 × simple clone operations = much faster

### **Expected Results**
```python
# Before fix:
139,958 records/sec with WRONG chromosome names ("chr_0", "chr_1")

# After fix: 
~160,000 records/sec with CORRECT chromosome names ("chr1", "chr2", "chrX")
```

**Expected Improvement**: 10-20% performance gain + correct data output

## 🚀 **Additional High-Impact Optimizations**

Once the header resolution is fixed, these are the next highest-impact optimizations:

### **1. String Interning for Repeated Values** (Potential 1.15x gain)
```rust
use string_cache::DefaultAtom;
use std::collections::HashMap;

struct StringCache {
    chroms: HashMap<i32, DefaultAtom>,
    quality_patterns: HashMap<Vec<u8>, DefaultAtom>,
}

// Reuse common strings like "chr1", "chr2", repeated quality patterns
```

### **2. SIMD-Accelerated Base Decoding** (Potential 1.2-1.3x gain)
```rust
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

fn decode_bases_simd(sequence_bytes: &[u8]) -> String {
    // Process 16 bytes (32 bases) at once using AVX2
    // Replace the current byte-by-byte loop
}
```

### **3. Zero-Copy Quality Score Processing** (Potential 1.1x gain)
```rust
// Replace the current char-by-char conversion with:
fn quality_to_string_fast(qual: &[u8]) -> String {
    let mut result = String::with_capacity(qual.len());
    unsafe {
        let bytes = result.as_mut_vec();
        bytes.extend_from_slice(qual);
        for byte in bytes {
            *byte += b'!';
        }
    }
    result
}
```

### **4. Adaptive Batch Sizing** (Potential 1.1x gain)
```rust
struct AdaptiveBatcher {
    target_memory_per_batch: usize,
    avg_record_size: f64,
    current_batch_size: usize,
}

// Automatically adjust batch sizes based on record sizes
```

## 🎯 **Implementation Priority**

### **Phase 1: Critical Fix (Immediate - 1 day)**
1. **Header resolution fix** - Correctness + 15% performance gain

### **Phase 2: High-Impact Optimizations (1-2 weeks)**
2. **String interning** - 15% additional gain
3. **SIMD base decoding** - 20-30% for sequence-heavy workloads
4. **Zero-copy quality scores** - 10-15% gain

### **Phase 3: Advanced Optimizations (1 month)**
5. **Adaptive batch sizing** - 10-15% gain
6. **Memory pooling** - Reduced GC pressure
7. **Lock-free queues** - 5-10% synchronization improvement

## 🏆 **Potential Combined Impact**

Implementing all optimizations could achieve:

```python
# Current: 139,958 records/sec
# Phase 1: ~160,000 records/sec (+14% - header fix)
# Phase 2: ~200,000 records/sec (+25% - string optimizations)  
# Phase 3: ~240,000 records/sec (+20% - advanced optimizations)

# Total potential: 240,000 records/sec (1.7x additional improvement)
# Overall system: 6x improvement over original 40k baseline
```

## 🔧 **Testing the Fix**

After implementing the header resolution fix:

1. **Correctness Test**:
```python
# Check that output contains real chromosome names
import polars as pl
df = pl.read_ipc("output.arrow")
print(df['chrom'].unique())  # Should show ["chr1", "chr2", "chrX"] not ["chr_0", "chr_1"]
```

2. **Performance Test**:
```python
# Run simple benchmark to confirm ~15% improvement
python simple_benchmark.py --bam-path test.bam
```

## 🚨 **Why This Is Critical Priority**

1. **Data Correctness**: Currently outputting WRONG chromosome names!
2. **Easy Fix**: ~30 lines of code changes
3. **Big Performance Gain**: 10-20% improvement  
4. **Zero Risk**: Simple, well-understood optimization
5. **Foundation**: Enables other string optimizations

The header resolution fix should be implemented immediately as it addresses both correctness and performance issues with minimal risk and effort.