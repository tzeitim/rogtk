#!/usr/bin/env python3
"""
Simple HTSlib Parallel Benchmark Script

This script runs a focused set of parameter tests and displays results
in a clean, easy-to-read table format.

Usage:
    python simple_benchmark.py --bam-path /path/to/file.bam
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


class SimpleBenchmark:
    def __init__(self, bam_path: str):
        self.bam_path = bam_path
        self.results = []
        
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
    
    def run_test(self, name: str, **params):
        """Run a single test configuration"""
        print(f"Running: {name}...", end=" ", flush=True)
        
        # Create temporary output file
        output_file = f"temp_output_{int(time.time() * 1000) % 100000}.arrow.ipc"
        
        try:
            start_time = time.time()
            
            rogtk.bam_to_arrow_ipc_htslib_parallel(
                bam_path=self.bam_path,
                arrow_ipc_path=output_file,
                include_sequence=True,
                include_quality=True,
                limit=1_000_000,  # 1M records for quick testing
                **params
            )
            
            duration = time.time() - start_time
            throughput = 1_000_000 / duration
            
            # Clean up
            if os.path.exists(output_file):
                os.unlink(output_file)
            
            result = {
                'name': name,
                'throughput': throughput,
                'duration': duration,
                'status': 'SUCCESS',
                **params
            }
            
            print(f"{throughput:,.0f} rec/sec ({duration:.1f}s)")
            
        except Exception as e:
            # Clean up
            if os.path.exists(output_file):
                os.unlink(output_file)
            
            result = {
                'name': name,
                'throughput': 0,
                'duration': 0,
                'status': f'FAILED: {str(e)[:50]}',
                **params
            }
            
            print(f"FAILED: {str(e)[:50]}")
        
        self.results.append(result)
        return result
    
    def run_focused_benchmark(self):
        """Run a focused set of tests to find optimal configuration"""
        print("=" * 80)
        print("HTSlib Parallel Performance Benchmark")
        print(f"BAM file: {self.bam_path}")
        print(f"Test size: 1,000,000 records per test")
        print("=" * 80)
        print()
        
        # Test configurations
        tests = [
            # Baseline tests
            ("Baseline (small)", {
                "batch_size": 50_000,
                "bgzf_threads": 2,
                "writing_threads": 1,
                "read_buffer_mb": 64,
                "write_buffer_mb": 32
            }),
            
            # Current best known config
            ("Current Best", {
                "batch_size": 50_000,
                "bgzf_threads": 4,
                "writing_threads": 8,
                "read_buffer_mb": 256,
                "write_buffer_mb": 64
            }),
            
            # Batch size variations
            ("Small Batches", {
                "batch_size": 25_000,
                "bgzf_threads": 4,
                "writing_threads": 8,
                "read_buffer_mb": 256,
                "write_buffer_mb": 64
            }),
            
            ("Large Batches", {
                "batch_size": 100_000,
                "bgzf_threads": 4,
                "writing_threads": 8,
                "read_buffer_mb": 256,
                "write_buffer_mb": 64
            }),
            
            # Thread scaling tests
            ("More Processing", {
                "batch_size": 50_000,
                "bgzf_threads": 4,
                "writing_threads": 12,
                "read_buffer_mb": 256,
                "write_buffer_mb": 64
            }),
            
            ("Max Processing", {
                "batch_size": 50_000,
                "bgzf_threads": 4,
                "writing_threads": 16,
                "read_buffer_mb": 256,
                "write_buffer_mb": 64
            }),
            
            # BGZF thread tests
            ("More BGZF", {
                "batch_size": 50_000,
                "bgzf_threads": 6,
                "writing_threads": 8,
                "read_buffer_mb": 256,
                "write_buffer_mb": 64
            }),
            
            ("Max BGZF", {
                "batch_size": 50_000,
                "bgzf_threads": 8,
                "writing_threads": 8,
                "read_buffer_mb": 256,
                "write_buffer_mb": 64
            }),
            
            # Buffer size tests
            ("Huge Buffers", {
                "batch_size": 50_000,
                "bgzf_threads": 4,
                "writing_threads": 8,
                "read_buffer_mb": 512,
                "write_buffer_mb": 128
            }),
            
            ("Small Buffers", {
                "batch_size": 50_000,
                "bgzf_threads": 4,
                "writing_threads": 8,
                "read_buffer_mb": 128,
                "write_buffer_mb": 32
            }),
            
            # Extreme configurations
            ("High Parallel", {
                "batch_size": 25_000,
                "bgzf_threads": 6,
                "writing_threads": 12,
                "read_buffer_mb": 512,
                "write_buffer_mb": 128
            }),
            
            ("Conservative", {
                "batch_size": 100_000,
                "bgzf_threads": 2,
                "writing_threads": 4,
                "read_buffer_mb": 128,
                "write_buffer_mb": 64
            }),
        ]
        
        # Run all tests
        for test_name, params in tests:
            self.run_test(test_name, **params)
            time.sleep(1)  # Brief pause between tests
        
        print()
        self.display_results()
    
    def display_results(self):
        """Display results in a clean table format"""
        print("=" * 120)
        print("BENCHMARK RESULTS SUMMARY")
        print("=" * 120)
        
        # Sort by throughput
        successful_results = [r for r in self.results if r['status'] == 'SUCCESS']
        failed_results = [r for r in self.results if r['status'] != 'SUCCESS']
        
        successful_results.sort(key=lambda x: x['throughput'], reverse=True)
        
        if successful_results:
            print("\n🏆 SUCCESSFUL CONFIGURATIONS (sorted by performance):")
            print()
            print("| Rank | Configuration    | Throughput (rec/sec) | Batch Size | BGZF | Writing | Read Buf | Write Buf | Duration |")
            print("|------|------------------|---------------------|------------|------|---------|----------|-----------|----------|")
            
            for i, result in enumerate(successful_results, 1):
                print(f"| {i:2d}   | {result['name']:<16} | {result['throughput']:>15,.0f} | "
                      f"{result['batch_size']:>7,} | {result['bgzf_threads']:>4} | "
                      f"{result['writing_threads']:>7} | {result['read_buffer_mb']:>6}MB | "
                      f"{result['write_buffer_mb']:>7}MB | {result['duration']:>6.1f}s |")
            
            # Performance analysis
            best = successful_results[0]
            worst = successful_results[-1]
            improvement = best['throughput'] / worst['throughput']
            
            print(f"\n📊 PERFORMANCE ANALYSIS:")
            print(f"   Best:        {best['name']} - {best['throughput']:,.0f} records/sec")
            print(f"   Worst:       {worst['name']} - {worst['throughput']:,.0f} records/sec")
            print(f"   Improvement: {improvement:.1f}x faster")
            print(f"   Range:       {worst['throughput']:,.0f} - {best['throughput']:,.0f} records/sec")
            
            # Optimal configuration
            print(f"\n🎯 OPTIMAL CONFIGURATION:")
            print(f"   Configuration: {best['name']}")
            print(f"   Expected throughput: {best['throughput']:,.0f} records/sec")
            print(f"   Time for 2M records: ~{2_000_000/best['throughput']:.1f} seconds")
            print()
            print("   Python code:")
            print("   rogtk.bam_to_arrow_ipc_htslib_parallel(")
            print(f"       batch_size={best['batch_size']:,},")
            print(f"       bgzf_threads={best['bgzf_threads']},")
            print(f"       writing_threads={best['writing_threads']},")
            print(f"       read_buffer_mb={best['read_buffer_mb']},")
            print(f"       write_buffer_mb={best['write_buffer_mb']},")
            print("   )")
            
            # Parameter insights
            print(f"\n💡 KEY INSIGHTS:")
            
            # Analyze batch size impact
            batch_sizes = {}
            for r in successful_results:
                bs = r['batch_size']
                if bs not in batch_sizes:
                    batch_sizes[bs] = []
                batch_sizes[bs].append(r['throughput'])
            
            best_batch = max(batch_sizes.items(), key=lambda x: sum(x[1])/len(x[1]))
            print(f"   Optimal batch size: {best_batch[0]:,} (avg {sum(best_batch[1])/len(best_batch[1]):,.0f} rec/sec)")
            
            # Analyze thread impact
            writing_threads = {}
            for r in successful_results:
                wt = r['writing_threads']
                if wt not in writing_threads:
                    writing_threads[wt] = []
                writing_threads[wt].append(r['throughput'])
            
            best_writing = max(writing_threads.items(), key=lambda x: sum(x[1])/len(x[1]))
            print(f"   Optimal writing threads: {best_writing[0]} (avg {sum(best_writing[1])/len(best_writing[1]):,.0f} rec/sec)")
            
            # Analyze buffer impact
            read_buffers = {}
            for r in successful_results:
                rb = r['read_buffer_mb']
                if rb not in read_buffers:
                    read_buffers[rb] = []
                read_buffers[rb].append(r['throughput'])
            
            best_read_buf = max(read_buffers.items(), key=lambda x: sum(x[1])/len(x[1]))
            print(f"   Optimal read buffer: {best_read_buf[0]}MB (avg {sum(best_read_buf[1])/len(best_read_buf[1]):,.0f} rec/sec)")
        
        if failed_results:
            print(f"\n❌ FAILED CONFIGURATIONS:")
            for result in failed_results:
                print(f"   {result['name']}: {result['status']}")
        
        print("\n" + "=" * 120)


def main():
    parser = argparse.ArgumentParser(description="Simple HTSlib parallel benchmark")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    
    args = parser.parse_args()
    
    try:
        benchmark = SimpleBenchmark(args.bam_path)
        benchmark.run_focused_benchmark()
        
    except KeyboardInterrupt:
        print("\n\nBenchmark interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()