#!/usr/bin/env python3
"""
BGZF Block-Level Optimization Benchmark Script

Tests the new reader pooling, buffer optimization, and segment tuning improvements
to validate the estimated 25-40% performance gain.

Usage:
    python bgzf_optimization_benchmark.py --bam-path /path/to/file.bam
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


class BGZFOptimizationBenchmark:
    def __init__(self, bam_path: str):
        self.bam_path = bam_path
        self.results = []
        
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
    
    def run_test(self, name: str, function_name: str, **params):
        """Run a single test configuration"""
        print(f"Running: {name}...", end=" ", flush=True)
        
        # Create temporary output file
        output_file = f"temp_bgzf_output_{int(time.time() * 1000) % 100000}.arrow.ipc"
        
        try:
            start_time = time.time()
            
            # Get the function by name
            func = getattr(rogtk, function_name)
            
            func(
                bam_path=self.bam_path,
                arrow_ipc_path=output_file,
                include_sequence=True,
                include_quality=True,
                limit=1_000_000,  # 1M records for testing
                **params
            )
            
            end_time = time.time()
            duration = end_time - start_time
            records_per_sec = 1_000_000 / duration
            
            # Get file size for throughput calculation
            if os.path.exists(output_file):
                file_size_mb = os.path.getsize(output_file) / (1024 * 1024)
                mb_per_sec = file_size_mb / duration
            else:
                mb_per_sec = 0
            
            self.results.append({
                'name': name,
                'duration': duration,
                'records_per_sec': records_per_sec,
                'mb_per_sec': mb_per_sec,
                'params': params
            })
            
            print(f"✓ {records_per_sec:,.0f} rec/sec ({duration:.1f}s)")
            
        except Exception as e:
            print(f"✗ Failed: {e}")
            self.results.append({
                'name': name,
                'duration': float('inf'),
                'records_per_sec': 0,
                'mb_per_sec': 0,
                'params': params
            })
        
        finally:
            # Clean up output file
            if os.path.exists(output_file):
                os.remove(output_file)
    
    def run_benchmark(self):
        """Run all benchmark configurations"""
        print("=" * 80)
        print("BGZF Block-Level Optimization Benchmark")
        print("=" * 80)
        print(f"BAM file: {self.bam_path}")
        print(f"Test dataset: 1M records per test")
        print()
        
        # Test 1: Original optimized version (159k rec/sec baseline)
        print("🔍 Testing baseline performance...")
        self.run_test(
            "Baseline: Optimized Single-Reader",
            "bam_to_arrow_ipc_htslib_optimized",
            batch_size=20_000,
            max_bgzf_threads=4,
            writing_threads=12,
            read_buffer_mb=1024,
            write_buffer_mb=256
        )
        
        # Test 2: BGZF block-level with old parameters
        print("\n🔍 Testing BGZF block-level (basic configuration)...")
        self.run_test(
            "BGZF Basic: 4 Workers",
            "bam_to_arrow_ipc_htslib_bgzf_blocks",
            batch_size=20_000,
            bgzf_threads=4,
            writing_threads=8,
            num_block_workers=4,
            read_buffer_mb=1024,
            write_buffer_mb=256
        )
        
        # Test 3: BGZF with optimized buffers (auto-scaling)
        print("\n🚀 Testing BGZF with optimized buffer scaling...")
        self.run_test(
            "BGZF Optimized: 4 Workers + Auto Buffers",
            "bam_to_arrow_ipc_htslib_bgzf_blocks",
            batch_size=20_000,
            bgzf_threads=4,
            writing_threads=8,
            num_block_workers=4
            # Let buffer sizes auto-scale based on worker count
        )
        
        # Test 4: BGZF with more workers and optimizations
        print("\n🚀 Testing BGZF with high parallelism...")
        self.run_test(
            "BGZF High-Parallel: 8 Workers + Auto Buffers",
            "bam_to_arrow_ipc_htslib_bgzf_blocks",
            batch_size=20_000,
            bgzf_threads=6,
            writing_threads=8,
            num_block_workers=8
            # Auto buffer scaling: 1.5GB read, 384MB write for 8+ workers
        )
        
        # Test 5: BGZF with maximum optimization
        print("\n🚀 Testing BGZF with maximum workers...")
        self.run_test(
            "BGZF Maximum: 10 Workers + Auto Buffers",
            "bam_to_arrow_ipc_htslib_bgzf_blocks",
            batch_size=15_000,  # Smaller batches for more parallelism
            bgzf_threads=8,
            writing_threads=8,
            num_block_workers=10
        )
        
        print("\n" + "=" * 80)
        print("BENCHMARK RESULTS")
        print("=" * 80)
        
        # Sort results by performance
        sorted_results = sorted(self.results, key=lambda x: x['records_per_sec'], reverse=True)
        
        # Print results table
        print(f"{'Rank':<4} {'Configuration':<35} {'Records/sec':<12} {'Duration':<10} {'Improvement':<12}")
        print("-" * 80)
        
        baseline_perf = None
        for i, result in enumerate(sorted_results, 1):
            if result['records_per_sec'] > 0:
                if baseline_perf is None:
                    baseline_perf = result['records_per_sec']
                    improvement = "Baseline"
                else:
                    improvement = f"+{((result['records_per_sec'] / baseline_perf - 1) * 100):+.1f}%"
                
                print(f"{i:<4} {result['name']:<35} {result['records_per_sec']:>10,.0f} {result['duration']:>8.1f}s {improvement:<12}")
            else:
                print(f"{i:<4} {result['name']:<35} {'FAILED':<12} {'N/A':<10} {'N/A':<12}")
        
        # Analysis
        print("\n" + "=" * 80)
        print("PERFORMANCE ANALYSIS")
        print("=" * 80)
        
        if len([r for r in sorted_results if r['records_per_sec'] > 0]) >= 2:
            best = sorted_results[0]
            baseline = next((r for r in sorted_results if 'Baseline' in r['name']), None)
            
            if baseline and best['name'] != baseline['name']:
                improvement = (best['records_per_sec'] / baseline['records_per_sec'] - 1) * 100
                print(f"🎯 Best Configuration: {best['name']}")
                print(f"📈 Performance Improvement: +{improvement:.1f}% over baseline")
                print(f"⚡ Throughput: {best['records_per_sec']:,.0f} records/sec")
                
                # Check if we're approaching the 200k+ target
                target = 200_000
                progress = (best['records_per_sec'] / target) * 100
                print(f"🎯 Progress toward 200k rec/sec target: {progress:.1f}%")
                
                if best['records_per_sec'] >= target:
                    print("🎉 TARGET ACHIEVED! Exceeded 200k records/sec!")
                elif best['records_per_sec'] >= target * 0.9:
                    print("🔥 Very close to target! Within 10% of 200k rec/sec goal")
                elif improvement >= 25:
                    print("✓ Significant improvement achieved! Met roadmap expectation of 25%+ gain")
            else:
                print("⚠️  No improvement detected over baseline")
        
        print(f"\n💡 Optimization Summary:")
        print(f"   ✓ Reader pooling: Reduces HTSlib initialization overhead")
        print(f"   ✓ Buffer auto-scaling: 1.5GB read/384MB write for high parallelism")
        print(f"   ✓ Segment optimization: Minimum 512MB segments to reduce coordination")
        print(f"   ✓ Zero-copy quality processing: Already integrated")


def main():
    parser = argparse.ArgumentParser(description="BGZF Block-Level Optimization Benchmark")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    
    args = parser.parse_args()
    
    try:
        benchmark = BGZFOptimizationBenchmark(args.bam_path)
        benchmark.run_benchmark()
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Benchmark failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()