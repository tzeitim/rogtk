#!/usr/bin/env python3
"""
Advanced Optimization Benchmark
===============================

Tests the new zero-copy hybrid optimized approach against the 40% efficiency target.
Expected performance leap: 309k rec/sec → 500-600k rec/sec (20.3% → 40%+ efficiency)

Key optimizations tested:
1. Zero-copy memory pools
2. SIMD-optimized batch processing  
3. Proper multi-segment concatenation
4. Lock-free result aggregation
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
    """Run a single test configuration with detailed performance tracking"""
    print(f"Running: {name} ({limit:,} records)...", end=" ", flush=True)
    
    output_file = f"temp_advanced_{int(time.time() * 1000) % 100000}.arrow.ipc"
    
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
            'params': params,
            'error': str(e)
        }
    
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def main():
    parser = argparse.ArgumentParser(description="Advanced Optimization Benchmark")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    parser.add_argument("--limit", type=int, default=2_000_000, help="Number of records to test (default: 2M)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam_path):
        print(f"ERROR: BAM file not found: {args.bam_path}")
        sys.exit(1)
    
    print("=" * 80)
    print("ADVANCED OPTIMIZATION BENCHMARK")
    print("=" * 80)
    print(f"BAM file: {args.bam_path}")
    
    file_size_gb = os.path.getsize(args.bam_path) / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")
    print(f"Test limit: {args.limit:,} records")
    print()
    
    print("🎯 TESTING HYPOTHESIS: Advanced Optimizations → 40%+ Efficiency")
    print("   Target: 500-600k rec/sec (up from 309k baseline)")
    print("   Key optimizations: Zero-copy pools, SIMD batching, proper concatenation")
    print()
    
    results = []
    
    # Optimal parameters from previous optimization
    optimal_params = {
        'batch_size': 15000,
        'read_buffer_mb': 2048,
        'write_buffer_mb': 128,
        'max_bgzf_threads': 4,
        'writing_threads': 6
    }
    
    print("📊 BASELINE COMPARISON")
    print("-" * 40)
    
    # Test 1: Original hybrid segments (baseline from parameter optimization)
    baseline = run_test(
        "Original Hybrid Segments",
        "bam_to_arrow_ipc_htslib_hybrid_segments",
        args.bam_path,
        args.limit,
        **optimal_params
    )
    results.append(baseline)
    
    print(f"\n🚀 ADVANCED OPTIMIZATIONS")
    print("-" * 40)
    
    # Test 2: New zero-copy hybrid optimized approach
    optimized = run_test(
        "Zero-Copy Hybrid Optimized",
        "bam_to_arrow_ipc_htslib_hybrid_optimized",
        args.bam_path,
        args.limit,
        **optimal_params
    )
    results.append(optimized)
    
    # Test 3: Auto-segmentation with optimizations
    auto_optimized = run_test(
        "Auto-Segment Optimized",
        "bam_to_arrow_ipc_htslib_hybrid_optimized",
        args.bam_path,
        args.limit,
        **{**optimal_params, 'num_segments': None}  # Auto-determine segments
    )
    results.append(auto_optimized)
    
    # Test 4: High-segment optimized (push the limits)
    high_segment = run_test(
        "High-Segment Optimized",
        "bam_to_arrow_ipc_htslib_hybrid_optimized", 
        args.bam_path,
        args.limit,
        **{**optimal_params, 'num_segments': 8}  # Force 8 segments
    )
    results.append(high_segment)
    
    # Analysis
    print("\n" + "=" * 80)
    print("ADVANCED OPTIMIZATION PERFORMANCE ANALYSIS")
    print("=" * 80)
    
    # Sort by performance
    valid_results = [r for r in results if r['records_per_sec'] > 0]
    sorted_results = sorted(valid_results, key=lambda x: x['records_per_sec'], reverse=True)
    
    print(f"{'Rank':<4} {'Configuration':<30} {'Records/sec':<12} {'vs Baseline':<12} {'Efficiency':<12}")
    print("-" * 80)
    
    baseline_perf = baseline['records_per_sec'] if baseline['records_per_sec'] > 0 else 1
    theoretical_max = 153_000  # Per segment theoretical max
    
    for i, result in enumerate(sorted_results, 1):
        improvement = (result['records_per_sec'] / baseline_perf - 1) * 100
        
        # Calculate efficiency assuming segments
        segments = result['params'].get('num_segments')
        if segments is None:
            # Auto-determined segments based on file size
            if file_size_gb < 1.0:
                segments = 1
            elif file_size_gb < 10.0:
                segments = 2
            elif file_size_gb < 50.0:
                segments = 4
            else:
                segments = 8
        elif segments == 0:  # Original might not have this param
            segments = 4  # Reasonable default
        
        efficiency = (result['records_per_sec'] / (theoretical_max * segments)) * 100
        
        print(f"{i:<4} {result['name']:<30} {result['records_per_sec']:>10,.0f} {improvement:>+9.1f}% {efficiency:>9.1f}%")
    
    print(f"\n💡 OPTIMIZATION IMPACT ANALYSIS:")
    
    # Initialize default values
    best_efficiency = 0
    best = None
    
    if len(valid_results) >= 1:
        best = sorted_results[0]
        baseline_result = baseline
        
        print(f"   📈 Best Configuration: {best['name']}")
        print(f"   📊 Performance: {best['records_per_sec']:,.0f} rec/sec")
        
        if baseline_result['records_per_sec'] > 0:
            improvement = (best['records_per_sec'] / baseline_result['records_per_sec'] - 1) * 100
            print(f"   🚀 Improvement: {improvement:+.1f}% over baseline")
        
        # Efficiency calculation for best result
        best_segments = best['params'].get('num_segments')
        if best_segments is None:
            if file_size_gb < 1.0:
                best_segments = 1
            elif file_size_gb < 10.0:
                best_segments = 2
            elif file_size_gb < 50.0:
                best_segments = 4
            else:
                best_segments = 8
        elif best_segments == 0:
            best_segments = 4
            
        best_efficiency = (best['records_per_sec'] / (theoretical_max * best_segments)) * 100
        print(f"   ⚡ Efficiency: {best_efficiency:.1f}% per segment")
        
        # Target achievement
        target_achievement = (best['records_per_sec'] / 500_000) * 100
        print(f"   🎯 Target Achievement: {target_achievement:.1f}% of 500k rec/sec goal")
        
        if best['records_per_sec'] >= 500_000:
            print(f"   🎉 BREAKTHROUGH! Exceeded 500k records/sec target!")
            if best_efficiency >= 40.0:
                print(f"   🏆 EFFICIENCY TARGET ACHIEVED! Reached {best_efficiency:.1f}% (40%+ goal)")
        elif best['records_per_sec'] >= 400_000:
            print(f"   ⭐ MAJOR PROGRESS! Close to 500k target")
    else:
        print(f"   ⚠️  All optimization tests failed - need to fix implementation issues")
        
    print(f"\n🔬 OPTIMIZATION BREAKDOWN:")
    print(f"   • Zero-copy memory pools: Eliminate record cloning overhead")
    print(f"   • SIMD-optimized batching: Vectorized field extraction")
    print(f"   • Proper concatenation: Use all segment outputs")
    print(f"   • Lock-free aggregation: Reduce thread synchronization")
    print(f"   • File size: {file_size_gb:.2f}GB optimal for parallel processing")
    
    # Save detailed results
    results_file = f"advanced_optimization_results_{int(time.time())}.json"
    
    import json
    detailed_results = {
        'test_info': {
            'bam_file': args.bam_path,
            'test_limit': args.limit,
            'file_size_gb': file_size_gb,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'optimization_target': '40%+ efficiency (500-600k rec/sec)'
        },
        'baseline_performance': {
            'records_per_sec': baseline_perf,
            'efficiency_percent': (baseline_perf / (theoretical_max * 4)) * 100  # Assume 4 segments
        },
        'results': results,
        'best_result': sorted_results[0] if sorted_results else None,
        'target_achievement': {
            'efficiency_target': 40.0,
            'throughput_target': 500_000,
            'achieved_efficiency': best_efficiency if valid_results else 0,
            'achieved_throughput': sorted_results[0]['records_per_sec'] if sorted_results else 0
        }
    }
    
    with open(results_file, 'w') as f:
        json.dump(detailed_results, f, indent=2)
    
    print(f"\n💾 Detailed results saved to: {results_file}")
    
    if sorted_results and sorted_results[0]['records_per_sec'] >= 500_000:
        print(f"\n🎊 MISSION ACCOMPLISHED!")
        print(f"   Advanced optimizations successfully achieved 500k+ rec/sec target!")
    else:
        print(f"\n🔄 NEXT STEPS:")
        print(f"   Consider additional optimizations:")
        print(f"   • Custom Arrow builders for faster serialization")
        print(f"   • Memory-mapped I/O for large files")
        print(f"   • CPU-specific SIMD intrinsics")


if __name__ == "__main__":
    main()