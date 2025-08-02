#!/usr/bin/env python3
"""
Minimal Fix Test
================

Tests if the ONLY critical fix (proper segment concatenation) 
can achieve the performance boost we need.

Expected: 258k → 400-500k rec/sec just from using ALL segments
"""

import argparse
import time
import os
import sys

try:
    import rogtk
except ImportError:
    print("ERROR: rogtk module not found. Please install with 'maturin develop --features htslib --release'")
    sys.exit(1)


def run_test(name, function_name, bam_path, limit, **params):
    """Run a single test configuration"""
    print(f"Running: {name} ({limit:,} records)...", end=" ", flush=True)
    
    output_file = f"temp_minimal_{int(time.time() * 1000) % 100000}.arrow.ipc"
    
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
        
        print(f"✓ {records_per_sec:,.0f} rec/sec ({duration:.1f}s)")
        return records_per_sec
        
    except Exception as e:
        print(f"✗ Failed: {e}")
        return 0
    
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def main():
    parser = argparse.ArgumentParser(description="Minimal Fix Test")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam_path):
        print(f"ERROR: BAM file not found: {args.bam_path}")
        sys.exit(1)
    
    print("=" * 60)
    print("MINIMAL FIX TEST: Segment Concatenation Only")
    print("=" * 60)
    print(f"BAM file: {args.bam_path}")
    
    file_size_gb = os.path.getsize(args.bam_path) / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")
    print()
    
    limit = 2_000_000
    
    # Optimal parameters from previous testing
    optimal_params = {
        'batch_size': 15000,
        'read_buffer_mb': 2048,
        'write_buffer_mb': 128,
        'max_bgzf_threads': 4,
        'writing_threads': 6
    }
    
    print("🧪 HYPOTHESIS: Single fix should boost 258k → 400-500k rec/sec")
    print("   Fix: Use ALL segment outputs instead of just the first")
    print()
    
    # Test 1: Original (baseline)
    print("📊 BASELINE:")
    original = run_test(
        "Original Hybrid (broken concatenation)",
        "bam_to_arrow_ipc_htslib_hybrid_segments",
        args.bam_path,
        limit,
        **optimal_params
    )
    
    print(f"\n🔧 MINIMAL FIX:")
    # Test 2: Minimal fix
    fixed = run_test(
        "Minimal Fix (proper concatenation)",
        "bam_to_arrow_ipc_htslib_hybrid_minimal_fix",
        args.bam_path,
        limit,
        **optimal_params
    )
    
    print(f"\n" + "=" * 60)
    print("MINIMAL FIX ANALYSIS")
    print("=" * 60)
    
    if original > 0 and fixed > 0:
        improvement = (fixed / original - 1) * 100
        print(f"Original: {original:,.0f} rec/sec")
        print(f"Fixed:    {fixed:,.0f} rec/sec")
        print(f"Improvement: {improvement:+.1f}%")
        
        if fixed >= 400_000:
            print(f"🎉 SUCCESS! Achieved 400k+ rec/sec target")
        elif fixed >= 350_000:
            print(f"⭐ GREAT! Major improvement, close to target")
        elif improvement >= 50:
            print(f"✅ GOOD! Significant improvement from simple fix")
        else:
            print(f"⚠️  Minimal improvement - may need other fixes")
            
        # Efficiency calculation
        segments = 8 if file_size_gb > 50 else 4
        efficiency = (fixed / (153_000 * segments)) * 100
        print(f"Efficiency: {efficiency:.1f}% per segment")
        
    elif fixed > 0:
        print(f"Fixed version works: {fixed:,.0f} rec/sec")
        print(f"Original version failed")
    else:
        print(f"Both versions failed - need debugging")
    
    print(f"\n💡 Key Insight:")
    print(f"   The original hybrid was only using the first segment's output")
    print(f"   This minimal fix ensures ALL {8 if file_size_gb > 50 else 4} segments contribute to results")


if __name__ == "__main__":
    main()