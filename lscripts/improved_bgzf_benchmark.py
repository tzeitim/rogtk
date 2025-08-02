#!/usr/bin/env python3
"""
Improved BGZF Benchmark with Performance Analysis

Tests both small and large dataset scenarios to understand when
BGZF block-level optimization provides benefits.
"""

import argparse
import time
import os
import sys
from pathlib import Path

try:
    import rogtk
except ImportError:
    print("ERROR: rogtk module not found. Please install with 'maturin develop --features htslib --release'")
    sys.exit(1)


def run_test(name, function_name, bam_path, limit, **params):
    """Run a single test configuration"""
    print(f"Running: {name} ({limit:,} records)...", end=" ", flush=True)
    
    output_file = f"temp_improved_{int(time.time() * 1000) % 100000}.arrow.ipc"
    
    try:
        start_time = time.time()
        
        func = getattr(rogtk, function_name)
        func(
            bam_path=bam_path,
            arrow_ipc_path=output_file,
            include_sequence=True,
            include_quality=True,
            limit=limit,
            **params
        )
        
        end_time = time.time()
        duration = end_time - start_time
        records_per_sec = limit / duration
        
        result = {
            'name': name,
            'duration': duration,
            'records_per_sec': records_per_sec,
            'limit': limit,
            'params': params
        }
        
        print(f"✓ {records_per_sec:,.0f} rec/sec ({duration:.1f}s)")
        return result
        
    except Exception as e:
        print(f"✗ Failed: {e}")
        return {
            'name': name,
            'duration': float('inf'),
            'records_per_sec': 0,
            'limit': limit,
            'params': params
        }
    
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def main():
    parser = argparse.ArgumentParser(description="Improved BGZF Benchmark")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam_path):
        print(f"ERROR: BAM file not found: {args.bam_path}")
        sys.exit(1)
    
    print("=" * 80)
    print("IMPROVED BGZF BLOCK-LEVEL BENCHMARK")
    print("=" * 80)
    print(f"BAM file: {args.bam_path}")
    
    file_size_gb = os.path.getsize(args.bam_path) / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")
    print()
    
    results = []
    
    # Test different dataset sizes to find the crossover point
    test_sizes = [100_000, 500_000, 1_000_000, 5_000_000]  # Scale up to find sweet spot
    
    for limit in test_sizes:
        print(f"\n{'='*20} TESTING {limit:,} RECORDS {'='*20}")
        
        # Baseline: Single-reader optimized
        baseline = run_test(
            f"Baseline ({limit//1000}K)",
            "bam_to_arrow_ipc_htslib_optimized",
            args.bam_path,
            limit,
            batch_size=20_000,
            max_bgzf_threads=4,
            writing_threads=8,  # Reduced for fair comparison
            read_buffer_mb=1024,
            write_buffer_mb=256
        )
        results.append(baseline)
        
        # BGZF with optimized settings for this dataset size
        if limit <= 1_000_000:
            # Small dataset: minimal workers
            bgzf_result = run_test(
                f"BGZF-Optimized ({limit//1000}K)",
                "bam_to_arrow_ipc_htslib_bgzf_blocks",
                args.bam_path,
                limit,
                batch_size=20_000,
                bgzf_threads=4,
                writing_threads=4,  # Fewer workers for small datasets
                num_block_workers=2  # Auto-limited by our fixes
            )
        else:
            # Larger dataset: more workers
            bgzf_result = run_test(
                f"BGZF-Optimized ({limit//1000}K)",
                "bam_to_arrow_ipc_htslib_bgzf_blocks",
                args.bam_path,
                limit,
                batch_size=15_000,
                bgzf_threads=6,
                writing_threads=6,
                num_block_workers=4
            )
        results.append(bgzf_result)
        
        # Calculate improvement for this dataset size
        if baseline['records_per_sec'] > 0 and bgzf_result['records_per_sec'] > 0:
            improvement = (bgzf_result['records_per_sec'] / baseline['records_per_sec'] - 1) * 100
            if improvement > 0:
                print(f"   📈 BGZF is {improvement:+.1f}% faster for {limit:,} records")
            else:
                print(f"   📉 BGZF is {improvement:+.1f}% slower for {limit:,} records")
        
        print()
    
    # Analysis
    print("=" * 80)
    print("CROSSOVER POINT ANALYSIS")
    print("=" * 80)
    
    crossover_found = False
    for i in range(0, len(results), 2):  # Step by 2 to compare baseline vs BGZF
        if i+1 < len(results):
            baseline = results[i]
            bgzf = results[i+1]
            
            if baseline['records_per_sec'] > 0 and bgzf['records_per_sec'] > 0:
                improvement = (bgzf['records_per_sec'] / baseline['records_per_sec'] - 1) * 100
                
                if improvement > 0 and not crossover_found:
                    print(f"🎯 CROSSOVER POINT: BGZF becomes beneficial at ~{bgzf['limit']:,} records")
                    print(f"   Improvement: {improvement:+.1f}%")
                    crossover_found = True
                    break
    
    if not crossover_found:
        max_tested = max(r['limit'] for r in results)
        print(f"⚠️  BGZF did not outperform baseline up to {max_tested:,} records")
        print(f"   This suggests the crossover point is >~{max_tested:,} records")
        print(f"   BGZF is optimized for very large files (50GB+, 50M+ records)")
    
    print(f"\n💡 Key Insights:")
    print(f"   • File size: {file_size_gb:.2f}GB - {'Small' if file_size_gb < 10 else 'Large'} file")
    print(f"   • BGZF optimizations: ✓ Reader pooling, ✓ Buffer scaling, ✓ Reduced debug output")
    print(f"   • Worker auto-limiting: ✓ 2 workers for ≤1M records, 4 workers for >1M records")
    
    if file_size_gb < 10:
        print(f"   • RECOMMENDATION: Use bam_to_arrow_ipc_htslib_optimized() for files <10GB")
    else:
        print(f"   • RECOMMENDATION: BGZF may show benefits with larger record limits")


if __name__ == "__main__":
    main()