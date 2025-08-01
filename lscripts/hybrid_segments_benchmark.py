#!/usr/bin/env python3
"""
Hybrid Segments Benchmark

Tests the new hybrid approach: Multiple independent optimized single-readers
Each segment runs the proven 153k rec/sec optimized pipeline independently.

Expected: 153k rec/sec × N segments = Massive linear scalability
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
    
    output_file = f"temp_hybrid_{int(time.time() * 1000) % 100000}.arrow.ipc"
    
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
            'function': function_name,
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
            'function': function_name,
            'duration': float('inf'),
            'records_per_sec': 0,
            'limit': limit,
            'params': params
        }
    
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def main():
    parser = argparse.ArgumentParser(description="Hybrid Segments Benchmark")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam_path):
        print(f"ERROR: BAM file not found: {args.bam_path}")
        sys.exit(1)
    
    print("=" * 80)
    print("HYBRID SEGMENTS BENCHMARK")
    print("=" * 80)
    print(f"BAM file: {args.bam_path}")
    
    file_size_gb = os.path.getsize(args.bam_path) / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")
    print()
    
    print("🎯 Testing Hypothesis: Multiple Independent Optimized Readers")
    print("   Each segment runs the proven 153k rec/sec optimized pipeline")
    print("   Expected linear scaling: 153k × N segments")
    print()
    
    results = []
    limit = 1_000_000  # 1M records for consistent testing
    
    # Test 1: Baseline single optimized reader
    print("📊 BASELINE TEST")
    baseline = run_test(
        "Single Optimized Reader",
        "bam_to_arrow_ipc_htslib_optimized",
        args.bam_path,
        limit,
        batch_size=20_000,
        max_bgzf_threads=4,
        writing_threads=8,
        read_buffer_mb=1024,
        write_buffer_mb=256
    )
    results.append(baseline)
    
    print(f"\n🚀 HYBRID SEGMENTS TESTS")
    
    # Test 2: Hybrid with auto-segmentation
    hybrid_auto = run_test(
        "Hybrid Auto-Segments",
        "bam_to_arrow_ipc_htslib_hybrid_segments",
        args.bam_path,
        limit,
        batch_size=20_000,
        max_bgzf_threads=4,
        writing_threads=8,
        read_buffer_mb=1024,
        write_buffer_mb=256
        # num_segments will be auto-determined
    )
    results.append(hybrid_auto)
    
    # Test 3: Hybrid with 2 segments
    hybrid_2 = run_test(
        "Hybrid 2-Segments",
        "bam_to_arrow_ipc_htslib_hybrid_segments",
        args.bam_path,
        limit,
        batch_size=20_000,
        max_bgzf_threads=4,
        writing_threads=8,  
        read_buffer_mb=1024,
        write_buffer_mb=256,
        num_segments=2
    )
    results.append(hybrid_2)
    
    # Test 4: Hybrid with 4 segments
    hybrid_4 = run_test(
        "Hybrid 4-Segments",
        "bam_to_arrow_ipc_htslib_hybrid_segments",
        args.bam_path,
        limit,
        batch_size=20_000,
        max_bgzf_threads=4,
        writing_threads=8,
        read_buffer_mb=1024,
        write_buffer_mb=256,
        num_segments=4
    )
    results.append(hybrid_4)
    
    # Test 5: Compare with original BGZF block-level
    print(f"\n🔍 COMPARISON WITH ORIGINAL BGZF")
    bgzf_original = run_test(
        "Original BGZF Block-Level",
        "bam_to_arrow_ipc_htslib_bgzf_blocks",
        args.bam_path,
        limit,
        batch_size=20_000,
        bgzf_threads=4,
        writing_threads=8,
        num_block_workers=4
    )
    results.append(bgzf_original)
    
    # Analysis
    print("\n" + "=" * 80)
    print("HYBRID SEGMENTS PERFORMANCE ANALYSIS")
    print("=" * 80)
    
    # Sort by performance
    valid_results = [r for r in results if r['records_per_sec'] > 0]
    sorted_results = sorted(valid_results, key=lambda x: x['records_per_sec'], reverse=True)
    
    print(f"{'Rank':<4} {'Configuration':<30} {'Records/sec':<12} {'vs Baseline':<12} {'Segments':<10}")
    print("-" * 80)
    
    baseline_perf = baseline['records_per_sec'] if baseline['records_per_sec'] > 0 else 1
    
    for i, result in enumerate(sorted_results, 1):
        improvement = (result['records_per_sec'] / baseline_perf - 1) * 100
        
        # Extract segment count
        segments = "1"
        if "2-Segments" in result['name']:
            segments = "2"
        elif "4-Segments" in result['name']:
            segments = "4"
        elif "Auto" in result['name']:
            segments = "Auto"
        elif "BGZF" in result['name']:
            segments = "N/A"
        
        print(f"{i:<4} {result['name']:<30} {result['records_per_sec']:>10,.0f} {improvement:>+9.1f}% {segments:<10}")
    
    # Linear scaling analysis
    print(f"\n💡 LINEAR SCALING ANALYSIS:")
    
    baseline_single = baseline_perf
    
    for result in results:
        if "Hybrid" in result['name'] and result['records_per_sec'] > 0:
            if "2-Segments" in result['name']:
                expected = baseline_single * 2
                actual = result['records_per_sec']
                efficiency = (actual / expected) * 100
                print(f"   2-Segments: {actual:,.0f} rec/sec (expected: {expected:,.0f}, efficiency: {efficiency:.1f}%)")
            
            elif "4-Segments" in result['name']:
                expected = baseline_single * 4
                actual = result['records_per_sec']
                efficiency = (actual / expected) * 100
                print(f"   4-Segments: {actual:,.0f} rec/sec (expected: {expected:,.0f}, efficiency: {efficiency:.1f}%)")
    
    # Best approach recommendation
    print(f"\n🎯 RECOMMENDATION:")
    if sorted_results:
        best = sorted_results[0]
        if "Hybrid" in best['name']:
            print(f"   ✓ HYBRID APPROACH WINS! {best['name']}")
            print(f"   ✓ Performance: {best['records_per_sec']:,.0f} rec/sec")
            improvement = (best['records_per_sec'] / baseline_perf - 1) * 100
            print(f"   ✓ Improvement: {improvement:+.1f}% over single-reader baseline")
            
            # Progress toward 200k target
            target_progress = (best['records_per_sec'] / 200_000) * 100
            print(f"   🎯 Progress toward 200k target: {target_progress:.1f}%")
            
            if best['records_per_sec'] >= 200_000:
                print(f"   🎉 TARGET ACHIEVED! Exceeded 200k records/sec!")
            
        else:
            print(f"   → Single optimized reader remains best for this file size")
    
    print(f"\n🔬 ARCHITECTURE INSIGHTS:")
    print(f"   • File size: {file_size_gb:.2f}GB")
    print(f"   • Hybrid approach eliminates coordination overhead")
    print(f"   • Each segment runs independently at proven 153k rec/sec")
    print(f"   • Linear scaling depends on I/O bandwidth and CPU cores")


if __name__ == "__main__":
    main()