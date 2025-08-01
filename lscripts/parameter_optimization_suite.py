#!/usr/bin/env python3
"""
Comprehensive Parameter Optimization Suite for Hybrid Segments

This script systematically tests different parameter combinations to optimize
the hybrid segments approach from current 25.8% efficiency to 40%+ efficiency.

CONTEXT FOR NEW CLAUDE SESSION:
- We achieved 395k rec/sec with hybrid segments (197% of 200k target)
- Current efficiency: 25.8% per segment (39.5k rec/sec per segment vs 153k theoretical)
- Goal: Optimize to 40%+ efficiency → 500-600k rec/sec total throughput
- Architecture: Multiple independent optimized readers processing BAM segments in parallel

Usage: python parameter_optimization_suite.py --bam-path /path/to/file.bam
Output: Detailed JSON report with optimal parameters for maximum performance
"""

import rogtk
import time
import json
import sys
import os
import argparse
from datetime import datetime
from pathlib import Path
import itertools

class ParameterOptimizer:
    def __init__(self, bam_path, output_dir="optimization_results"):
        self.bam_path = bam_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Test dataset size - large enough to show meaningful performance differences
        self.test_limit = 2_000_000  # 2M records (proven to show best scaling)
        self.test_segments = 10      # 10 segments (proven optimal configuration)
        
        # Results storage
        self.results = {
            "test_info": {
                "bam_file": bam_path,
                "test_limit": self.test_limit,
                "test_segments": self.test_segments,
                "timestamp": datetime.now().isoformat(),
                "context": "Optimizing hybrid segments from 25.8% to 40%+ efficiency per segment"
            },
            "baseline_performance": {
                "current_total_throughput": 395140,  # rec/sec
                "current_per_segment_throughput": 39514,  # rec/sec  
                "current_efficiency_percent": 25.8,
                "theoretical_max_per_segment": 153000,  # rec/sec
                "target_efficiency_percent": 40.0,
                "target_total_throughput": 612000  # 153k * 10 segments * 40%
            },
            "optimization_phases": {},
            "best_configurations": {},
            "recommendations": {}
        }
    
    def test_configuration(self, phase_name, config_name, **params):
        """Test a specific parameter configuration"""
        print(f"  Testing {config_name}...", end=" ", flush=True)
        
        output_file = self.output_dir / f"test_{phase_name}_{config_name.replace(' ', '_')}.arrow"
        
        try:
            start_time = time.time()
            
            rogtk.bam_to_arrow_ipc_htslib_hybrid_segments(
                bam_path=self.bam_path,
                arrow_ipc_path=str(output_file),
                include_sequence=True,
                include_quality=True,
                limit=self.test_limit,
                num_segments=self.test_segments,
                **params
            )
            
            duration = time.time() - start_time
            throughput = self.test_limit / duration
            per_segment_throughput = throughput / self.test_segments
            efficiency = (per_segment_throughput / 153000) * 100
            
            result = {
                "config_name": config_name,
                "parameters": params,
                "duration_seconds": round(duration, 2),
                "total_throughput_rec_per_sec": int(throughput),
                "per_segment_throughput_rec_per_sec": int(per_segment_throughput),
                "efficiency_percent": round(efficiency, 1),
                "improvement_vs_baseline": round(((throughput / 395140) - 1) * 100, 1),
                "status": "success"
            }
            
            print(f"✓ {throughput:,.0f} rec/sec ({efficiency:.1f}% efficiency)")
            
            # Cleanup
            if output_file.exists():
                output_file.unlink()
                
            return result
            
        except Exception as e:
            print(f"✗ Failed: {e}")
            
            # Cleanup on error
            if output_file.exists():
                output_file.unlink()
                
            return {
                "config_name": config_name,
                "parameters": params,
                "status": "failed",
                "error": str(e)
            }
    
    def phase_7a_buffer_optimization(self):
        """Phase 7A: Buffer Size Optimization"""
        print("\n🔧 PHASE 7A: Buffer Size Optimization")
        print("=" * 60)
        print("Testing different read/write buffer combinations for optimal I/O performance")
        print()
        
        # Buffer size combinations to test (in MB)
        read_buffers = [512, 1024, 1536, 2048]    # 512MB to 2GB
        write_buffers = [128, 256, 384, 512]      # 128MB to 512MB
        
        # Baseline configuration for comparison
        baseline_params = {
            "batch_size": 25000,
            "max_bgzf_threads": 6,
            "writing_threads": 10
        }
        
        results = []
        
        for read_mb, write_mb in itertools.product(read_buffers, write_buffers):
            config_name = f"Read{read_mb}MB_Write{write_mb}MB"
            
            test_params = baseline_params.copy()
            test_params.update({
                "read_buffer_mb": read_mb,
                "write_buffer_mb": write_mb
            })
            
            result = self.test_configuration("7A_buffer", config_name, **test_params)
            results.append(result)
        
        self.results["optimization_phases"]["phase_7a_buffer_optimization"] = {
            "description": "Testing read/write buffer size combinations",
            "test_matrix": {
                "read_buffers_mb": read_buffers,
                "write_buffers_mb": write_buffers,
                "total_combinations": len(results)
            },
            "results": results,
            "best_config": max([r for r in results if r["status"] == "success"], 
                             key=lambda x: x["total_throughput_rec_per_sec"], default=None)
        }
        
        print(f"\n📊 Phase 7A Summary: Tested {len(results)} buffer combinations")
        return results
    
    def phase_7b_thread_optimization(self):
        """Phase 7B: Thread Configuration Tuning"""
        print("\n🔧 PHASE 7B: Thread Configuration Tuning")
        print("=" * 60)
        print("Testing BGZF vs writing thread combinations for optimal CPU utilization")
        print()
        
        # Use best buffer configuration from Phase 7A
        best_buffer_config = self.results["optimization_phases"]["phase_7a_buffer_optimization"]["best_config"]
        if not best_buffer_config:
            print("⚠️  Using default buffers (Phase 7A had no successful results)")
            buffer_params = {"read_buffer_mb": 1024, "write_buffer_mb": 256}
        else:
            buffer_params = {k: v for k, v in best_buffer_config["parameters"].items() 
                           if k in ["read_buffer_mb", "write_buffer_mb"]}
        
        # Thread combinations to test
        bgzf_threads = [4, 6, 8, 10]       # BGZF decompression threads
        writing_threads = [6, 8, 10, 12]   # Writing/processing threads
        
        baseline_params = {
            "batch_size": 25000,
            **buffer_params
        }
        
        results = []
        
        for bgzf_t, write_t in itertools.product(bgzf_threads, writing_threads):
            # Skip configurations where total threads exceed reasonable limits
            if bgzf_t + write_t > 20:
                continue
                
            config_name = f"BGZF{bgzf_t}_Write{write_t}"
            
            test_params = baseline_params.copy()
            test_params.update({
                "max_bgzf_threads": bgzf_t,
                "writing_threads": write_t
            })
            
            result = self.test_configuration("7B_threads", config_name, **test_params)
            results.append(result)
        
        self.results["optimization_phases"]["phase_7b_thread_optimization"] = {
            "description": "Testing BGZF vs writing thread combinations",
            "inherited_from_7a": buffer_params,
            "test_matrix": {
                "bgzf_threads": bgzf_threads,
                "writing_threads": writing_threads,
                "total_combinations": len(results)
            },
            "results": results,
            "best_config": max([r for r in results if r["status"] == "success"], 
                             key=lambda x: x["total_throughput_rec_per_sec"], default=None)
        }
        
        print(f"\n📊 Phase 7B Summary: Tested {len(results)} thread combinations")
        return results
    
    def phase_7c_batch_optimization(self):
        """Phase 7C: Batch Size Optimization"""
        print("\n🔧 PHASE 7C: Batch Size Optimization")
        print("=" * 60)
        print("Testing different batch sizes for optimal memory usage and processing efficiency")
        print()
        
        # Use best configuration from Phase 7B
        best_thread_config = self.results["optimization_phases"]["phase_7b_thread_optimization"]["best_config"]
        if not best_thread_config:
            print("⚠️  Using default threads (Phase 7B had no successful results)")
            base_params = {
                "max_bgzf_threads": 6,
                "writing_threads": 10,
                "read_buffer_mb": 1024,
                "write_buffer_mb": 256
            }
        else:
            base_params = best_thread_config["parameters"].copy()
        
        # Batch sizes to test
        batch_sizes = [10000, 15000, 20000, 25000, 30000, 40000, 50000]
        
        results = []
        
        for batch_size in batch_sizes:
            config_name = f"Batch{batch_size}"
            
            test_params = base_params.copy()
            test_params.update({"batch_size": batch_size})
            
            result = self.test_configuration("7C_batch", config_name, **test_params)
            results.append(result)
        
        self.results["optimization_phases"]["phase_7c_batch_optimization"] = {
            "description": "Testing different batch sizes for optimal processing",
            "inherited_from_7b": {k: v for k, v in base_params.items() if k != "batch_size"},
            "test_matrix": {
                "batch_sizes": batch_sizes,
                "total_combinations": len(results)
            },
            "results": results,
            "best_config": max([r for r in results if r["status"] == "success"], 
                             key=lambda x: x["total_throughput_rec_per_sec"], default=None)
        }
        
        print(f"\n📊 Phase 7C Summary: Tested {len(results)} batch sizes")
        return results
    
    def analyze_results(self):
        """Analyze all results and provide recommendations"""
        print("\n🎯 FINAL ANALYSIS")
        print("=" * 60)
        
        # Find overall best configuration
        all_successful_results = []
        for phase_name, phase_data in self.results["optimization_phases"].items():
            successful_results = [r for r in phase_data["results"] if r["status"] == "success"]
            all_successful_results.extend(successful_results)
        
        if not all_successful_results:
            print("❌ No successful optimization results found")
            return
        
        # Best overall configuration
        best_overall = max(all_successful_results, key=lambda x: x["total_throughput_rec_per_sec"])
        
        # Calculate improvements
        baseline_throughput = 395140
        improvement = ((best_overall["total_throughput_rec_per_sec"] / baseline_throughput) - 1) * 100
        
        self.results["best_configurations"]["overall_best"] = best_overall
        self.results["best_configurations"]["improvement_analysis"] = {
            "baseline_throughput": baseline_throughput,
            "optimized_throughput": best_overall["total_throughput_rec_per_sec"],
            "improvement_percent": round(improvement, 1),
            "efficiency_improvement": f"{25.8:.1f}% → {best_overall['efficiency_percent']:.1f}%",
            "target_achievement": f"{(best_overall['efficiency_percent'] / 40.0) * 100:.1f}% of 40% target"
        }
        
        # Recommendations
        recommendations = []
        
        if best_overall["efficiency_percent"] >= 40:
            recommendations.append("🎉 TARGET ACHIEVED! Efficiency exceeded 40% target")
        elif best_overall["efficiency_percent"] >= 35:
            recommendations.append("🔥 Very close to 40% efficiency target")
        elif best_overall["efficiency_percent"] >= 30:
            recommendations.append("✅ Solid improvement, consider advanced optimizations next")
        else:
            recommendations.append("⚠️  Consider additional optimization approaches")
        
        if improvement >= 50:
            recommendations.append("🚀 Exceptional optimization gains achieved")
        elif improvement >= 25:
            recommendations.append("⚡ Strong optimization improvements")
        elif improvement >= 10:
            recommendations.append("📈 Meaningful performance gains")
        
        self.results["recommendations"]["next_steps"] = recommendations
        self.results["recommendations"]["production_config"] = best_overall["parameters"]
        
        # Print summary
        print(f"🏆 BEST CONFIGURATION FOUND:")
        print(f"   Throughput: {best_overall['total_throughput_rec_per_sec']:,} rec/sec")
        print(f"   Efficiency: {best_overall['efficiency_percent']:.1f}% per segment")
        print(f"   Improvement: +{improvement:.1f}% vs baseline")
        print(f"   Configuration: {best_overall['config_name']}")
        print()
        print("📋 OPTIMAL PARAMETERS:")
        for param, value in best_overall["parameters"].items():
            print(f"   {param}: {value}")
        print()
        print("💡 RECOMMENDATIONS:")
        for rec in recommendations:
            print(f"   {rec}")
    
    def save_results(self):
        """Save detailed results to JSON file"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = self.output_dir / f"parameter_optimization_results_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        print(f"\n💾 Results saved to: {output_file}")
        print(f"📁 Working directory: {self.output_dir.absolute()}")
        
        # Also save a summary file for quick reference
        summary_file = self.output_dir / f"optimization_summary_{timestamp}.txt"
        with open(summary_file, 'w') as f:
            f.write("HYBRID SEGMENTS PARAMETER OPTIMIZATION SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Test Date: {self.results['test_info']['timestamp']}\n")
            f.write(f"BAM File: {self.results['test_info']['bam_file']}\n")
            f.write(f"Test Dataset: {self.results['test_info']['test_limit']:,} records\n")
            f.write(f"Test Segments: {self.results['test_info']['test_segments']}\n\n")
            
            if "overall_best" in self.results["best_configurations"]:
                best = self.results["best_configurations"]["overall_best"]
                improvement = self.results["best_configurations"]["improvement_analysis"]
                
                f.write("BEST CONFIGURATION FOUND:\n")
                f.write("-" * 30 + "\n")
                f.write(f"Throughput: {best['total_throughput_rec_per_sec']:,} rec/sec\n")
                f.write(f"Per-segment: {best['per_segment_throughput_rec_per_sec']:,} rec/sec\n")
                f.write(f"Efficiency: {best['efficiency_percent']:.1f}% per segment\n")
                f.write(f"Improvement: +{improvement['improvement_percent']:.1f}% vs baseline\n\n")
                
                f.write("OPTIMAL PARAMETERS:\n")
                f.write("-" * 20 + "\n")
                for param, value in best["parameters"].items():
                    f.write(f"{param}: {value}\n")
                
                f.write(f"\nFor new Claude session: Use these parameters with hybrid segments\n")
                f.write(f"to achieve {best['total_throughput_rec_per_sec']:,} rec/sec performance.\n")
        
        return output_file
    
    def run_full_optimization(self):
        """Run complete parameter optimization suite"""
        print("🚀 HYBRID SEGMENTS PARAMETER OPTIMIZATION SUITE")
        print("=" * 70)
        print(f"BAM File: {self.bam_path}")
        print(f"Test Dataset: {self.test_limit:,} records with {self.test_segments} segments")
        print(f"Current Performance: 395k rec/sec (25.8% efficiency per segment)")
        print(f"Target: 40%+ efficiency → 500-600k rec/sec total throughput")
        print()
        
        # Run optimization phases
        try:
            self.phase_7a_buffer_optimization()
            self.phase_7b_thread_optimization()  
            self.phase_7c_batch_optimization()
            
            self.analyze_results()
            output_file = self.save_results()
            
            print(f"\n🎯 OPTIMIZATION COMPLETE!")
            print(f"📊 Full results: {output_file}")
            print(f"🔧 Ready for Phase 8 advanced optimizations")
            
        except Exception as e:
            print(f"\n❌ Optimization failed: {e}")
            # Still save partial results
            self.save_results()
            raise

def main():
    parser = argparse.ArgumentParser(description="Hybrid Segments Parameter Optimization Suite")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file for testing")
    parser.add_argument("--output-dir", default="optimization_results", 
                       help="Directory for output files (default: optimization_results)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam_path):
        print(f"❌ BAM file not found: {args.bam_path}")
        sys.exit(1)
    
    optimizer = ParameterOptimizer(args.bam_path, args.output_dir)
    optimizer.run_full_optimization()

if __name__ == "__main__":
    main()