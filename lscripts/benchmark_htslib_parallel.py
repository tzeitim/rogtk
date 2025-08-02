#!/usr/bin/env python3
"""
HTSlib Parallel Parameter Exploration Script

This script systematically tests different parameter combinations for the 
bam_to_arrow_ipc_htslib_parallel function to find optimal configurations
and document performance characteristics.

Usage:
    python benchmark_htslib_parallel.py --bam-path /path/to/file.bam --output-dir ./benchmarks

Results are saved to CSV and markdown files for analysis.
"""

import argparse
import csv
import itertools
import os
import time
from pathlib import Path
import sys
from typing import Dict, List, Tuple, Any
import json

try:
    import rogtk
except ImportError:
    print("ERROR: rogtk module not found. Please install with 'maturin develop --features htslib --release'")
    sys.exit(1)


class BenchmarkConfig:
    """Configuration for benchmark parameters to test"""
    
    # Parameter ranges to test
    BATCH_SIZES = [10_000, 25_000, 50_000, 100_000, 200_000, 500_000]
    BGZF_THREADS = [1, 2, 4, 6, 8]
    WRITING_THREADS = [1, 2, 4, 6, 8, 12, 16]
    READ_BUFFER_MB = [32, 64, 128, 256, 512]
    WRITE_BUFFER_MB = [16, 32, 64, 128, 256]
    
    # Test limits
    QUICK_TEST_LIMIT = 500_000   # For quick parameter sweep
    FULL_TEST_LIMIT = 2_000_000  # For detailed analysis of top configs
    
    @classmethod
    def get_quick_test_combinations(cls) -> List[Dict[str, Any]]:
        """Generate a reduced set of combinations for quick testing"""
        # Test key combinations first
        combinations = []
        
        # Core configurations to test
        base_configs = [
            # Baseline tests
            {"batch_size": 50_000, "bgzf_threads": 2, "writing_threads": 1, "read_buffer_mb": 64, "write_buffer_mb": 32},
            {"batch_size": 50_000, "bgzf_threads": 4, "writing_threads": 1, "read_buffer_mb": 128, "write_buffer_mb": 64},
            
            # Threading scaling tests
            {"batch_size": 50_000, "bgzf_threads": 4, "writing_threads": 4, "read_buffer_mb": 256, "write_buffer_mb": 64},
            {"batch_size": 50_000, "bgzf_threads": 4, "writing_threads": 8, "read_buffer_mb": 256, "write_buffer_mb": 64},
            {"batch_size": 50_000, "bgzf_threads": 4, "writing_threads": 12, "read_buffer_mb": 256, "write_buffer_mb": 64},
            
            # Batch size scaling tests  
            {"batch_size": 25_000, "bgzf_threads": 4, "writing_threads": 8, "read_buffer_mb": 256, "write_buffer_mb": 64},
            {"batch_size": 100_000, "bgzf_threads": 4, "writing_threads": 8, "read_buffer_mb": 256, "write_buffer_mb": 64},
            {"batch_size": 200_000, "bgzf_threads": 4, "writing_threads": 8, "read_buffer_mb": 256, "write_buffer_mb": 64},
            
            # Buffer size tests
            {"batch_size": 50_000, "bgzf_threads": 4, "writing_threads": 8, "read_buffer_mb": 128, "write_buffer_mb": 32},
            {"batch_size": 50_000, "bgzf_threads": 4, "writing_threads": 8, "read_buffer_mb": 512, "write_buffer_mb": 128},
            
            # BGZF thread tests
            {"batch_size": 50_000, "bgzf_threads": 1, "writing_threads": 8, "read_buffer_mb": 256, "write_buffer_mb": 64},
            {"batch_size": 50_000, "bgzf_threads": 6, "writing_threads": 8, "read_buffer_mb": 256, "write_buffer_mb": 64},
            {"batch_size": 50_000, "bgzf_threads": 8, "writing_threads": 8, "read_buffer_mb": 256, "write_buffer_mb": 64},
        ]
        
        return base_configs
    
    @classmethod 
    def get_full_test_combinations(cls, top_configs: List[Dict[str, Any]], num_variations: int = 5) -> List[Dict[str, Any]]:
        """Generate comprehensive test combinations around the best performing configs"""
        combinations = []
        
        for config in top_configs[:3]:  # Test variations of top 3 configs
            base_batch = config["batch_size"]
            base_bgzf = config["bgzf_threads"] 
            base_writing = config["writing_threads"]
            base_read_buf = config["read_buffer_mb"]
            base_write_buf = config["write_buffer_mb"]
            
            # Add the base config
            combinations.append(config.copy())
            
            # Vary each parameter around the optimal
            variations = [
                # Batch size variations
                {**config, "batch_size": max(10_000, base_batch // 2)},
                {**config, "batch_size": min(500_000, base_batch * 2)},
                
                # Thread variations  
                {**config, "bgzf_threads": max(1, base_bgzf - 1)},
                {**config, "bgzf_threads": min(8, base_bgzf + 1)},
                {**config, "writing_threads": max(1, base_writing - 2)},
                {**config, "writing_threads": min(16, base_writing + 2)},
                
                # Buffer variations
                {**config, "read_buffer_mb": base_read_buf // 2},
                {**config, "read_buffer_mb": min(512, base_read_buf * 2)},
                {**config, "write_buffer_mb": base_write_buf // 2}, 
                {**config, "write_buffer_mb": min(256, base_write_buf * 2)},
            ]
            
            combinations.extend(variations)
        
        # Remove duplicates
        seen = set()
        unique_combinations = []
        for combo in combinations:
            key = tuple(sorted(combo.items()))
            if key not in seen:
                seen.add(key)
                unique_combinations.append(combo)
        
        return unique_combinations


class BenchmarkRunner:
    """Runs benchmarks and collects results"""
    
    def __init__(self, bam_path: str, output_dir: str):
        self.bam_path = bam_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = []
        
        # Verify BAM file exists
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
    
    def run_single_benchmark(self, config: Dict[str, Any], limit: int) -> Dict[str, Any]:
        """Run a single benchmark configuration"""
        print(f"Testing: batch={config['batch_size']:,}, bgzf={config['bgzf_threads']}, "
              f"writing={config['writing_threads']}, read_buf={config['read_buffer_mb']}MB, "
              f"write_buf={config['write_buffer_mb']}MB")
        
        # Generate unique output filename
        timestamp = int(time.time() * 1000) % 100000
        output_file = self.output_dir / f"test_output_{timestamp}.arrow.ipc"
        
        try:
            start_time = time.time()
            
            # Run the benchmark
            rogtk.bam_to_arrow_ipc_htslib_parallel(
                bam_path=self.bam_path,
                arrow_ipc_path=str(output_file),
                batch_size=config["batch_size"],
                include_sequence=True,
                include_quality=True,
                bgzf_threads=config["bgzf_threads"],
                writing_threads=config["writing_threads"], 
                read_buffer_mb=config["read_buffer_mb"],
                write_buffer_mb=config["write_buffer_mb"],
                limit=limit
            )
            
            end_time = time.time()
            duration = end_time - start_time
            throughput = limit / duration
            
            # Get output file size
            file_size_mb = output_file.stat().st_size / (1024 * 1024) if output_file.exists() else 0
            
            result = {
                **config,
                "limit": limit,
                "duration_seconds": round(duration, 3),
                "throughput_records_per_sec": round(throughput, 0),
                "output_size_mb": round(file_size_mb, 2),
                "status": "success"
            }
            
            # Clean up output file
            if output_file.exists():
                output_file.unlink()
            
            print(f"  → {throughput:,.0f} records/sec in {duration:.1f}s")
            return result
            
        except Exception as e:
            print(f"  → FAILED: {str(e)}")
            result = {
                **config,
                "limit": limit,
                "duration_seconds": None,
                "throughput_records_per_sec": None,
                "output_size_mb": None,
                "status": f"failed: {str(e)}"
            }
            
            # Clean up output file if it exists
            if output_file.exists():
                output_file.unlink()
                
            return result
    
    def run_quick_benchmark(self) -> List[Dict[str, Any]]:
        """Run quick parameter sweep"""
        print("=" * 80)
        print("QUICK PARAMETER SWEEP")
        print("=" * 80)
        
        configs = BenchmarkConfig.get_quick_test_combinations()
        results = []
        
        for i, config in enumerate(configs, 1):
            print(f"\n[{i}/{len(configs)}] Quick test:")
            result = self.run_single_benchmark(config, BenchmarkConfig.QUICK_TEST_LIMIT)
            results.append(result)
            self.results.append(result)
            
            # Small delay to prevent overheating
            time.sleep(1)
        
        return results
    
    def run_detailed_benchmark(self, top_configs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Run detailed benchmark on top performing configurations"""
        print("\n" + "=" * 80)
        print("DETAILED BENCHMARKING OF TOP CONFIGURATIONS")
        print("=" * 80)
        
        configs = BenchmarkConfig.get_full_test_combinations(top_configs)
        results = []
        
        for i, config in enumerate(configs, 1):
            print(f"\n[{i}/{len(configs)}] Detailed test:")
            result = self.run_single_benchmark(config, BenchmarkConfig.FULL_TEST_LIMIT)
            results.append(result)
            self.results.append(result)
            
            # Longer delay for detailed tests
            time.sleep(2)
        
        return results
    
    def save_results(self):
        """Save results to CSV and JSON files"""
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        
        # Save to CSV
        csv_file = self.output_dir / f"benchmark_results_{timestamp}.csv"
        if self.results:
            with open(csv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=self.results[0].keys())
                writer.writeheader()
                writer.writerows(self.results)
            print(f"\nResults saved to: {csv_file}")
        
        # Save to JSON for detailed analysis
        json_file = self.output_dir / f"benchmark_results_{timestamp}.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"Results saved to: {json_file}")
        
        return csv_file, json_file


class ResultsAnalyzer:
    """Analyzes benchmark results and generates reports"""
    
    def __init__(self, results: List[Dict[str, Any]]):
        self.results = [r for r in results if r["status"] == "success"]
    
    def get_top_configs(self, n: int = 5) -> List[Dict[str, Any]]:
        """Get top N performing configurations"""
        sorted_results = sorted(self.results, 
                              key=lambda x: x["throughput_records_per_sec"] or 0, 
                              reverse=True)
        return sorted_results[:n]
    
    def analyze_parameter_impact(self, parameter: str) -> Dict[str, float]:
        """Analyze the impact of a specific parameter on performance"""
        param_groups = {}
        for result in self.results:
            param_value = result[parameter]
            if param_value not in param_groups:
                param_groups[param_value] = []
            param_groups[param_value].append(result["throughput_records_per_sec"])
        
        # Calculate average throughput for each parameter value
        param_averages = {}
        for param_value, throughputs in param_groups.items():
            valid_throughputs = [t for t in throughputs if t is not None]
            if valid_throughputs:
                param_averages[param_value] = sum(valid_throughputs) / len(valid_throughputs)
        
        return param_averages
    
    def generate_markdown_report(self, output_file: str):
        """Generate a comprehensive markdown report"""
        with open(output_file, 'w') as f:
            f.write("# HTSlib Parallel Performance Benchmark Results\n\n")
            f.write(f"**Generated on:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Total configurations tested:** {len(self.results)}\n\n")
            
            # Top performing configurations
            f.write("## 🏆 Top Performing Configurations\n\n")
            top_configs = self.get_top_configs(10)
            
            f.write("| Rank | Throughput (rec/sec) | Batch Size | BGZF Threads | Writing Threads | Read Buffer (MB) | Write Buffer (MB) | Duration (s) |\n")
            f.write("|------|---------------------|------------|--------------|-----------------|------------------|-------------------|---------------|\n")
            
            for i, config in enumerate(top_configs, 1):
                f.write(f"| {i} | {config['throughput_records_per_sec']:,.0f} | "
                       f"{config['batch_size']:,} | {config['bgzf_threads']} | "
                       f"{config['writing_threads']} | {config['read_buffer_mb']} | "
                       f"{config['write_buffer_mb']} | {config['duration_seconds']} |\n")
            
            # Parameter impact analysis
            f.write("\n## 📊 Parameter Impact Analysis\n\n")
            
            parameters = ["batch_size", "bgzf_threads", "writing_threads", "read_buffer_mb", "write_buffer_mb"]
            for param in parameters:
                f.write(f"### {param.replace('_', ' ').title()}\n\n")
                param_impact = self.analyze_parameter_impact(param)
                
                # Sort by performance
                sorted_impact = sorted(param_impact.items(), key=lambda x: x[1], reverse=True)
                
                f.write("| Value | Avg Throughput (rec/sec) | Relative Performance |\n")
                f.write("|-------|-------------------------|----------------------|\n")
                
                if sorted_impact:
                    best_throughput = sorted_impact[0][1]
                    for value, throughput in sorted_impact:
                        relative_perf = (throughput / best_throughput) * 100
                        f.write(f"| {value} | {throughput:,.0f} | {relative_perf:.1f}% |\n")
                
                f.write("\n")
            
            # Recommendations
            f.write("## 🎯 Recommendations\n\n")
            if top_configs:
                best = top_configs[0]
                f.write("### Optimal Configuration\n\n")
                f.write("```python\n")
                f.write("rogtk.bam_to_arrow_ipc_htslib_parallel(\n")
                f.write("    bam_path=bam_path,\n")
                f.write("    arrow_ipc_path=output_path,\n")
                f.write(f"    batch_size={best['batch_size']:,},\n")
                f.write(f"    bgzf_threads={best['bgzf_threads']},\n")
                f.write(f"    writing_threads={best['writing_threads']},\n")
                f.write(f"    read_buffer_mb={best['read_buffer_mb']},\n")
                f.write(f"    write_buffer_mb={best['write_buffer_mb']},\n")
                f.write(")\n")
                f.write("```\n\n")
                f.write(f"**Expected Performance:** {best['throughput_records_per_sec']:,.0f} records/sec\n\n")
            
            # All results table
            f.write("## 📋 Complete Results\n\n")
            f.write("| Batch Size | BGZF Threads | Writing Threads | Read Buffer (MB) | Write Buffer (MB) | Throughput (rec/sec) | Duration (s) | Status |\n")
            f.write("|------------|--------------|-----------------|------------------|-------------------|---------------------|--------------|--------|\n")
            
            # Sort all results by throughput
            all_sorted = sorted(self.results, key=lambda x: x["throughput_records_per_sec"] or 0, reverse=True)
            for result in all_sorted:
                throughput = f"{result['throughput_records_per_sec']:,.0f}" if result['throughput_records_per_sec'] else "FAILED"
                duration = f"{result['duration_seconds']}" if result['duration_seconds'] else "N/A"
                
                f.write(f"| {result['batch_size']:,} | {result['bgzf_threads']} | "
                       f"{result['writing_threads']} | {result['read_buffer_mb']} | "
                       f"{result['write_buffer_mb']} | {throughput} | {duration} | {result['status']} |\n")


def main():
    parser = argparse.ArgumentParser(description="Benchmark HTSlib parallel performance")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    parser.add_argument("--output-dir", default="./benchmarks", help="Output directory for results")
    parser.add_argument("--quick-only", action="store_true", help="Run only quick tests")
    parser.add_argument("--detailed-only", action="store_true", help="Run only detailed tests (requires previous results)")
    
    args = parser.parse_args()
    
    print("HTSlib Parallel Performance Benchmark")
    print("=" * 80)
    print(f"BAM file: {args.bam_path}")
    print(f"Output directory: {args.output_dir}")
    print()
    
    runner = BenchmarkRunner(args.bam_path, args.output_dir)
    
    if not args.detailed_only:
        # Run quick benchmark
        quick_results = runner.run_quick_benchmark()
        
        # Analyze quick results to find top performers
        analyzer = ResultsAnalyzer(quick_results)
        top_configs = analyzer.get_top_configs(5)
        
        print(f"\n🏆 Top 5 configurations from quick test:")
        for i, config in enumerate(top_configs, 1):
            print(f"{i}. {config['throughput_records_per_sec']:,.0f} rec/sec - "
                  f"batch={config['batch_size']:,}, bgzf={config['bgzf_threads']}, "
                  f"writing={config['writing_threads']}, buffers={config['read_buffer_mb']}/{config['write_buffer_mb']}MB")
        
        if not args.quick_only:
            # Run detailed benchmark on top configurations
            detailed_results = runner.run_detailed_benchmark(top_configs)
    
    # Save all results
    csv_file, json_file = runner.save_results()
    
    # Generate analysis report
    final_analyzer = ResultsAnalyzer(runner.results)
    report_file = runner.output_dir / f"performance_report_{time.strftime('%Y%m%d_%H%M%S')}.md"
    final_analyzer.generate_markdown_report(str(report_file))
    
    print(f"\n📊 Performance report generated: {report_file}")
    print(f"\n🎯 Best configuration found:")
    
    top_config = final_analyzer.get_top_configs(1)[0]
    print(f"  Throughput: {top_config['throughput_records_per_sec']:,.0f} records/sec")
    print(f"  Batch size: {top_config['batch_size']:,}")
    print(f"  BGZF threads: {top_config['bgzf_threads']}")
    print(f"  Writing threads: {top_config['writing_threads']}")
    print(f"  Read buffer: {top_config['read_buffer_mb']}MB")
    print(f"  Write buffer: {top_config['write_buffer_mb']}MB")


if __name__ == "__main__":
    main()