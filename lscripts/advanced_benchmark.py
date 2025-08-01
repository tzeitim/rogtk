#!/usr/bin/env python3
"""
Advanced HTSlib Parallel Benchmark Script

Based on findings that buffer sizes and small batches are key,
this script explores the parameter space around optimal configurations
to find the absolute maximum performance.

Key Findings to Explore:
- Huge buffers (512MB+) provide massive gains
- Small batches (25k) outperform large ones
- 8-12 writing threads optimal
- BGZF threads plateau at 4

Usage:
    python advanced_benchmark.py --bam-path /path/to/file.bam
"""

import argparse
import time
import os
import sys
from pathlib import Path
import itertools

try:
    import rogtk
except ImportError:
    print("ERROR: rogtk module not found. Please install with 'maturin develop --features htslib --release'")
    sys.exit(1)


class AdvancedBenchmark:
    def __init__(self, bam_path: str, test_size: int = 1_000_000):
        self.bam_path = bam_path
        self.test_size = test_size
        self.results = []
        
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
    
    def run_test(self, name: str, **params):
        """Run a single test configuration"""
        print(f"Testing: {name:<25}", end=" ", flush=True)
        
        # Create temporary output file
        output_file = f"temp_output_{int(time.time() * 1000) % 100000}.arrow.ipc"
        
        try:
            start_time = time.time()
            
            rogtk.bam_to_arrow_ipc_htslib_parallel(
                bam_path=self.bam_path,
                arrow_ipc_path=output_file,
                include_sequence=True,
                include_quality=True,
                limit=self.test_size,
                **params
            )
            
            duration = time.time() - start_time
            throughput = self.test_size / duration
            
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
            
            print(f"→ {throughput:>8,.0f} rec/sec ({duration:>5.1f}s)")
            
        except Exception as e:
            # Clean up
            if os.path.exists(output_file):
                os.unlink(output_file)
            
            result = {
                'name': name,
                'throughput': 0,
                'duration': 0,
                'status': f'FAILED: {str(e)[:30]}',
                **params
            }
            
            print(f"→ FAILED: {str(e)[:40]}")
        
        self.results.append(result)
        return result
    
    def run_buffer_exploration(self):
        """Explore extreme buffer sizes - the key optimization"""
        print("\n" + "=" * 80)
        print("🚀 BUFFER SIZE EXPLORATION (Key Optimization)")
        print("=" * 80)
        
        # Based on findings: huge buffers are the key!
        buffer_tests = [
            # Baseline from previous best
            ("Current Best", {
                "batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8,
                "read_buffer_mb": 512, "write_buffer_mb": 128
            }),
            
            # Even larger read buffers
            ("Massive Read 1GB", {
                "batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8,
                "read_buffer_mb": 1024, "write_buffer_mb": 128
            }),
            
            ("Extreme Read 2GB", {
                "batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8,
                "read_buffer_mb": 2048, "write_buffer_mb": 128
            }),
            
            # Larger write buffers
            ("Huge Write 256MB", {
                "batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8,
                "read_buffer_mb": 512, "write_buffer_mb": 256
            }),
            
            ("Massive Write 512MB", {
                "batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8,
                "read_buffer_mb": 512, "write_buffer_mb": 512
            }),
            
            # Both buffers extreme
            ("Both Massive", {
                "batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8,
                "read_buffer_mb": 1024, "write_buffer_mb": 256
            }),
            
            ("Both Extreme", {
                "batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8,
                "read_buffer_mb": 2048, "write_buffer_mb": 512
            }),
        ]
        
        for name, params in buffer_tests:
            self.run_test(name, **params)
            time.sleep(1)
    
    def run_batch_size_optimization(self):
        """Explore smaller batch sizes - second key finding"""
        print("\n" + "=" * 80)
        print("⚡ BATCH SIZE OPTIMIZATION (Parallelism vs Memory)")
        print("=" * 80)
        
        # Use best buffer config found so far
        best_buffer_config = {
            "bgzf_threads": 4, "writing_threads": 10,
            "read_buffer_mb": 1024, "write_buffer_mb": 256
        }
        
        batch_tests = [
            # Even smaller batches for more parallelism
            ("Tiny Batches 10k", {"batch_size": 10_000, **best_buffer_config}),
            ("Small Batches 15k", {"batch_size": 15_000, **best_buffer_config}),
            ("Optimal Batches 20k", {"batch_size": 20_000, **best_buffer_config}),
            ("Current Best 25k", {"batch_size": 25_000, **best_buffer_config}),
            ("Medium Batches 30k", {"batch_size": 30_000, **best_buffer_config}),
            ("Larger Batches 40k", {"batch_size": 40_000, **best_buffer_config}),
        ]
        
        for name, params in batch_tests:
            self.run_test(name, **params)
            time.sleep(1)
    
    def run_threading_fine_tuning(self):
        """Fine-tune threading around optimal range"""
        print("\n" + "=" * 80)
        print("🔧 THREADING FINE-TUNING (CPU Utilization)")
        print("=" * 80)
        
        # Use best config found so far
        best_config = {
            "batch_size": 20_000,  # Will adjust based on batch results
            "read_buffer_mb": 1024, "write_buffer_mb": 256
        }
        
        threading_tests = [
            # BGZF threading (I/O bound)
            ("BGZF 3 threads", {"bgzf_threads": 3, "writing_threads": 10, **best_config}),
            ("BGZF 4 threads", {"bgzf_threads": 4, "writing_threads": 10, **best_config}),
            ("BGZF 5 threads", {"bgzf_threads": 5, "writing_threads": 10, **best_config}),
            ("BGZF 6 threads", {"bgzf_threads": 6, "writing_threads": 10, **best_config}),
            
            # Writing threading (CPU bound)
            ("Writing 8 threads", {"bgzf_threads": 4, "writing_threads": 8, **best_config}),
            ("Writing 10 threads", {"bgzf_threads": 4, "writing_threads": 10, **best_config}),
            ("Writing 12 threads", {"bgzf_threads": 4, "writing_threads": 12, **best_config}),
            ("Writing 14 threads", {"bgzf_threads": 4, "writing_threads": 14, **best_config}),
            ("Writing 16 threads", {"bgzf_threads": 4, "writing_threads": 16, **best_config}),
        ]
        
        for name, params in threading_tests:
            self.run_test(name, **params)
            time.sleep(1)
    
    def run_extreme_configurations(self):
        """Test extreme configurations that might push limits"""
        print("\n" + "=" * 80)
        print("🚀 EXTREME CONFIGURATIONS (Pushing the Limits)")
        print("=" * 80)
        
        extreme_tests = [
            # Maximum reasonable configuration
            ("Maximum Reasonable", {
                "batch_size": 15_000, "bgzf_threads": 4, "writing_threads": 12,
                "read_buffer_mb": 2048, "write_buffer_mb": 512
            }),
            
            # Memory-optimized extreme
            ("Memory Extreme", {
                "batch_size": 10_000, "bgzf_threads": 6, "writing_threads": 16,
                "read_buffer_mb": 4096, "write_buffer_mb": 1024
            }),
            
            # Parallelism extreme
            ("Parallelism Extreme", {
                "batch_size": 5_000, "bgzf_threads": 4, "writing_threads": 20,
                "read_buffer_mb": 1024, "write_buffer_mb": 256
            }),
            
            # Balanced extreme
            ("Balanced Extreme", {
                "batch_size": 12_000, "bgzf_threads": 5, "writing_threads": 14,
                "read_buffer_mb": 1536, "write_buffer_mb": 384
            }),
            
            # Conservative extreme (less memory pressure)
            ("Conservative Extreme", {
                "batch_size": 20_000, "bgzf_threads": 3, "writing_threads": 8,
                "read_buffer_mb": 1024, "write_buffer_mb": 128
            }),
        ]
        
        for name, params in extreme_tests:
            self.run_test(name, **params)
            time.sleep(2)  # Longer pause for extreme tests
    
    def run_validation_tests(self):
        """Run validation tests on the absolute best configurations"""
        print("\n" + "=" * 80)
        print("✅ VALIDATION TESTS (Best Configurations)")
        print("=" * 80)
        
        # Get top 3 configurations from all previous tests
        successful = [r for r in self.results if r['status'] == 'SUCCESS']
        top_configs = sorted(successful, key=lambda x: x['throughput'], reverse=True)[:3]
        
        if not top_configs:
            print("No successful configurations to validate!")
            return
        
        print(f"Running validation tests on top {len(top_configs)} configurations...")
        print("Using larger test size for validation...")
        
        # Run validation with larger test size
        original_test_size = self.test_size
        self.test_size = 2_000_000  # 2M records for validation
        
        for i, config in enumerate(top_configs, 1):
            params = {k: v for k, v in config.items() 
                     if k not in ['name', 'throughput', 'duration', 'status']}
            self.run_test(f"Validation #{i}", **params)
            time.sleep(2)
        
        # Restore original test size
        self.test_size = original_test_size
    
    def analyze_and_display_results(self):
        """Comprehensive analysis and display of all results"""
        print("\n" + "=" * 120)
        print("📊 COMPREHENSIVE RESULTS ANALYSIS")
        print("=" * 120)
        
        successful = [r for r in self.results if r['status'] == 'SUCCESS']
        failed = [r for r in self.results if r['status'] != 'SUCCESS']
        
        if not successful:
            print("❌ No successful configurations found!")
            return
        
        # Sort by performance
        successful.sort(key=lambda x: x['throughput'], reverse=True)
        
        # Display top performers
        print("\n🏆 TOP 10 CONFIGURATIONS:")
        print(f"{'Rank':<4} {'Configuration':<25} {'Throughput (rec/sec)':<20} {'Batch':<8} {'BGZF':<5} {'Write':<6} {'Read Buf':<9} {'Write Buf':<10} {'Time':<6}")
        print("-" * 120)
        
        for i, result in enumerate(successful[:10], 1):
            print(f"{i:<4} {result['name']:<25} {result['throughput']:>15,.0f} "
                  f"{result['batch_size']:>7,} {result['bgzf_threads']:>4} "
                  f"{result['writing_threads']:>5} {result['read_buffer_mb']:>7}MB "
                  f"{result['write_buffer_mb']:>8}MB {result['duration']:>5.1f}s")
        
        # Performance statistics
        best = successful[0]
        worst = successful[-1]
        median = successful[len(successful)//2]
        
        print(f"\n📈 PERFORMANCE STATISTICS:")
        print(f"   Best:         {best['throughput']:>8,.0f} rec/sec ({best['name']})")
        print(f"   Worst:        {worst['throughput']:>8,.0f} rec/sec ({worst['name']})")
        print(f"   Median:       {median['throughput']:>8,.0f} rec/sec")
        print(f"   Improvement:  {best['throughput']/worst['throughput']:>8.1f}x")
        print(f"   Std Dev:      {self.calculate_std_dev([r['throughput'] for r in successful]):>8,.0f} rec/sec")
        
        # Parameter impact analysis
        self.analyze_parameter_impacts(successful)
        
        # Ultimate recommendation
        print(f"\n🎯 ULTIMATE OPTIMAL CONFIGURATION:")
        print(f"   Configuration: {best['name']}")
        print(f"   Peak Performance: {best['throughput']:,.0f} records/sec")
        print(f"   Estimated time for 10M records: {10_000_000/best['throughput']:.1f} seconds")
        print()
        print("   🐍 Python Code (Copy & Paste):")
        print("   rogtk.bam_to_arrow_ipc_htslib_parallel(")
        print("       bam_path=bam_path,")
        print("       arrow_ipc_path=output_path,")
        print(f"       batch_size={best['batch_size']:,},")
        print(f"       bgzf_threads={best['bgzf_threads']},")
        print(f"       writing_threads={best['writing_threads']},")
        print(f"       read_buffer_mb={best['read_buffer_mb']},")
        print(f"       write_buffer_mb={best['write_buffer_mb']},")
        print("   )")
        
        if failed:
            print(f"\n❌ FAILED CONFIGURATIONS ({len(failed)}):")
            for result in failed:
                print(f"   {result['name']}: {result['status']}")
        
        print("\n" + "=" * 120)
    
    def analyze_parameter_impacts(self, results):
        """Analyze the impact of each parameter"""
        print(f"\n🔍 PARAMETER IMPACT ANALYSIS:")
        
        params = ['batch_size', 'bgzf_threads', 'writing_threads', 'read_buffer_mb', 'write_buffer_mb']
        
        for param in params:
            param_groups = {}
            for r in results:
                value = r[param]
                if value not in param_groups:
                    param_groups[value] = []
                param_groups[value].append(r['throughput'])
            
            # Calculate averages and find best
            param_averages = {v: sum(throughputs)/len(throughputs) 
                            for v, throughputs in param_groups.items()}
            
            best_value = max(param_averages.items(), key=lambda x: x[1])
            worst_value = min(param_averages.items(), key=lambda x: x[1])
            
            print(f"   {param.replace('_', ' ').title():<20}: "
                  f"Best={best_value[0]} ({best_value[1]:,.0f} rec/sec), "
                  f"Worst={worst_value[0]} ({worst_value[1]:,.0f} rec/sec), "
                  f"Impact={best_value[1]/worst_value[1]:.1f}x")
    
    def calculate_std_dev(self, values):
        """Calculate standard deviation of throughput values"""
        if len(values) < 2:
            return 0
        mean = sum(values) / len(values)
        variance = sum((x - mean) ** 2 for x in values) / len(values)
        return variance ** 0.5
    
    def run_complete_benchmark(self):
        """Run the complete advanced benchmark suite"""
        print("🚀 ADVANCED HTSLIB PARALLEL BENCHMARK")
        print("=" * 80)
        print(f"BAM file: {self.bam_path}")
        print(f"Test size: {self.test_size:,} records per configuration")
        print("Exploring parameter space around key findings...")
        print()
        
        start_time = time.time()
        
        # Run all benchmark phases
        self.run_buffer_exploration()
        self.run_batch_size_optimization()
        self.run_threading_fine_tuning()
        self.run_extreme_configurations()
        self.run_validation_tests()
        
        total_time = time.time() - start_time
        
        print(f"\n⏱️  Total benchmark time: {total_time:.1f} seconds")
        print(f"📊 Total configurations tested: {len(self.results)}")
        
        # Final analysis
        self.analyze_and_display_results()


def main():
    parser = argparse.ArgumentParser(description="Advanced HTSlib parallel benchmark")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    parser.add_argument("--test-size", type=int, default=1_000_000, 
                       help="Number of records per test (default: 1M)")
    parser.add_argument("--quick", action="store_true", 
                       help="Run only buffer and batch optimization (skip extreme tests)")
    
    args = parser.parse_args()
    
    try:
        benchmark = AdvancedBenchmark(args.bam_path, args.test_size)
        
        if args.quick:
            print("🏃 QUICK MODE: Running buffer and batch optimization only")
            benchmark.run_buffer_exploration()
            benchmark.run_batch_size_optimization()
        else:
            benchmark.run_complete_benchmark()
        
        if not args.quick:
            benchmark.analyze_and_display_results()
            
    except KeyboardInterrupt:
        print("\n\nBenchmark interrupted by user")
        if benchmark.results:
            print("Displaying partial results...")
            benchmark.analyze_and_display_results()
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()