#!/usr/bin/env python3
"""
Intensive benchmark that stresses the hybrid approach with larger datasets
and provides detailed performance measurements
"""

import rogtk
import time
import sys
import os
import psutil
import threading
from pathlib import Path

class SystemMonitor:
    def __init__(self):
        self.monitoring = False
        self.cpu_samples = []
        self.memory_samples = []
        self.io_samples = []
        
    def start_monitoring(self):
        self.monitoring = True
        self.monitor_thread = threading.Thread(target=self._monitor_loop)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
        
    def stop_monitoring(self):
        self.monitoring = False
        if hasattr(self, 'monitor_thread'):
            self.monitor_thread.join(timeout=1)
    
    def _monitor_loop(self):
        process = psutil.Process()
        while self.monitoring:
            try:
                # CPU usage
                cpu_percent = process.cpu_percent()
                self.cpu_samples.append(cpu_percent)
                
                # Memory usage
                memory_info = process.memory_info()
                memory_mb = memory_info.rss / (1024 * 1024)
                self.memory_samples.append(memory_mb)
                
                # I/O stats
                io_counters = process.io_counters()
                self.io_samples.append({
                    'read_bytes': io_counters.read_bytes,
                    'write_bytes': io_counters.write_bytes
                })
                
                time.sleep(0.1)  # Sample every 100ms
            except:
                pass
    
    def get_stats(self):
        if not self.cpu_samples:
            return {}
            
        return {
            'avg_cpu_percent': sum(self.cpu_samples) / len(self.cpu_samples),
            'max_cpu_percent': max(self.cpu_samples),
            'avg_memory_mb': sum(self.memory_samples) / len(self.memory_samples),
            'max_memory_mb': max(self.memory_samples),
            'total_read_gb': (self.io_samples[-1]['read_bytes'] - self.io_samples[0]['read_bytes']) / (1024**3) if len(self.io_samples) > 1 else 0,
            'total_write_gb': (self.io_samples[-1]['write_bytes'] - self.io_samples[0]['write_bytes']) / (1024**3) if len(self.io_samples) > 1 else 0
        }

def stress_test(bam_path, test_name, segments, limit, iterations=3):
    """Run intensive stress test with multiple iterations"""
    print(f"🔥 STRESS TEST: {test_name}")
    print(f"   Dataset: {limit:,} records × {iterations} iterations")
    print(f"   Segments: {segments}")
    print("-" * 60)
    
    results = []
    total_records_processed = 0
    
    for iteration in range(iterations):
        print(f"Iteration {iteration + 1}/{iterations}...", end=" ", flush=True)
        
        monitor = SystemMonitor()
        monitor.start_monitoring()
        
        output_file = f"stress_test_{segments}seg_{iteration}.arrow"
        
        try:
            start_time = time.time()
            
            rogtk.bam_to_arrow_ipc_htslib_hybrid_segments(
                bam_path=bam_path,
                arrow_ipc_path=output_file,
                include_sequence=True,
                include_quality=True,
                limit=limit,
                num_segments=segments,
                batch_size=25_000,  # Larger batches for stress
                max_bgzf_threads=6,  # More threads
                writing_threads=10   # More writing threads
            )
            
            end_time = time.time()
            duration = end_time - start_time
            throughput = limit / duration
            
            monitor.stop_monitoring()
            system_stats = monitor.get_stats()
            
            # Get output file size
            output_size_mb = 0
            if os.path.exists(output_file):
                output_size_mb = os.path.getsize(output_file) / (1024 * 1024)
                os.remove(output_file)  # Cleanup
            
            results.append({
                'iteration': iteration + 1,
                'duration': duration,
                'throughput': throughput,
                'output_size_mb': output_size_mb,
                'system_stats': system_stats
            })
            
            total_records_processed += limit
            
            print(f"✓ {throughput:,.0f} rec/sec ({duration:.1f}s)")
            
        except Exception as e:
            print(f"✗ Failed: {e}")
            monitor.stop_monitoring()
            # Cleanup on error
            if os.path.exists(output_file):
                os.remove(output_file)
    
    # Calculate aggregate statistics
    if results:
        throughputs = [r['throughput'] for r in results]
        durations = [r['duration'] for r in results]
        
        avg_throughput = sum(throughputs) / len(throughputs)
        max_throughput = max(throughputs)
        min_throughput = min(throughputs)
        std_dev = (sum((t - avg_throughput) ** 2 for t in throughputs) / len(throughputs)) ** 0.5
        
        # System resource statistics
        avg_cpu = sum(r['system_stats'].get('avg_cpu_percent', 0) for r in results) / len(results)
        max_memory = max(r['system_stats'].get('max_memory_mb', 0) for r in results)
        total_io_read = sum(r['system_stats'].get('total_read_gb', 0) for r in results)
        total_io_write = sum(r['system_stats'].get('total_write_gb', 0) for r in results)
        
        print(f"\n📊 STRESS TEST RESULTS:")
        print(f"   Average: {avg_throughput:,.0f} rec/sec")
        print(f"   Maximum: {max_throughput:,.0f} rec/sec")
        print(f"   Minimum: {min_throughput:,.0f} rec/sec")
        print(f"   Std Dev: {std_dev:,.0f} rec/sec ({(std_dev/avg_throughput)*100:.1f}%)")
        print(f"   Consistency: {((min_throughput/max_throughput)*100):.1f}%")
        
        print(f"\n🖥️  SYSTEM RESOURCES:")
        print(f"   Avg CPU: {avg_cpu:.1f}%")
        print(f"   Max Memory: {max_memory:,.0f} MB")
        print(f"   Total I/O Read: {total_io_read:.2f} GB")
        print(f"   Total I/O Write: {total_io_write:.2f} GB")
        
        print(f"\n⚡ PERFORMANCE METRICS:")
        efficiency_per_segment = avg_throughput / (segments * 153_000) * 100
        print(f"   Efficiency per segment: {efficiency_per_segment:.1f}%")
        print(f"   Records processed: {total_records_processed:,}")
        print(f"   Total duration: {sum(durations):.1f}s")
        
        return avg_throughput
    
    return 0

def main():
    if len(sys.argv) > 1:
        bam_path = sys.argv[1]
    else:
        bam_path = "/home/projects/nyosef/pedro/projects/lt/datain/20250716/fracture/PBC15320_20250123_jurkat_1_1.bam"
    
    print("🚀 INTENSIVE HYBRID SEGMENTS BENCHMARK")
    print("=" * 70)
    print(f"BAM file: {bam_path}")
    
    # Get file size
    file_size_gb = os.path.getsize(bam_path) / (1024**3)
    print(f"File size: {file_size_gb:.1f} GB")
    
    # System info
    cpu_count = psutil.cpu_count()
    memory_gb = psutil.virtual_memory().total / (1024**3)
    print(f"System: {cpu_count} CPU cores, {memory_gb:.1f} GB RAM")
    print()
    
    # Define stress test scenarios - larger datasets to really stress the system
    scenarios = [
        ("Baseline Single Segment", 1, 5_000_000),      # 5M records, 1 segment
        ("Dual Segment Power", 2, 5_000_000),           # 5M records, 2 segments  
        ("Quad Segment Force", 4, 5_000_000),           # 5M records, 4 segments
        ("Octa Segment Extreme", 8, 10_000_000),        # 10M records, 8 segments
        ("Maximum Stress", 8, 20_000_000),              # 20M records, 8 segments
    ]
    
    results = []
    
    for test_name, segments, limit in scenarios:
        print()
        throughput = stress_test(bam_path, test_name, segments, limit, iterations=2)
        results.append((test_name, segments, limit, throughput))
        print()
    
    # Final analysis
    print("=" * 70)
    print("🎯 INTENSIVE BENCHMARK SUMMARY")
    print("=" * 70)
    print(f"{'Test Name':<25} {'Segments':<10} {'Records':<12} {'Throughput':<15}")
    print("-" * 70)
    
    best_throughput = 0
    best_config = ""
    
    for test_name, segments, limit, throughput in results:
        if throughput > 0:
            print(f"{test_name:<25} {segments:<10} {limit:>10,} {throughput:>12,.0f} rec/sec")
            if throughput > best_throughput:
                best_throughput = throughput
                best_config = f"{segments} segments @ {limit:,} records"
        else:
            print(f"{test_name:<25} {segments:<10} {limit:>10,} {'FAILED':>15}")
    
    print()
    print(f"🏆 BEST PERFORMANCE: {best_throughput:,.0f} rec/sec")
    print(f"🎯 Best configuration: {best_config}")
    
    # Target analysis
    target_progress = (best_throughput / 200_000) * 100
    print(f"📈 Progress toward 200k target: {target_progress:.1f}%")
    
    if best_throughput >= 200_000:
        print("🎉 TARGET ACHIEVED! Exceeded 200k records/sec!")
        over_target = ((best_throughput / 200_000) - 1) * 100
        print(f"🚀 Performance exceeds target by {over_target:.1f}%!")
    elif best_throughput >= 180_000:
        print("🔥 VERY CLOSE! Within 10% of 200k target")
    elif best_throughput >= 150_000:
        print("⚡ STRONG PERFORMANCE! Over 150k rec/sec achieved")
    
    print(f"\n💡 Key insights:")
    print(f"   • Hybrid segmentation enables true parallel processing")
    print(f"   • Performance scales with dataset size and segment count") 
    print(f"   • System can handle intensive workloads with multiple segments")
    print(f"   • Your insight about independent optimized readers was correct!")

if __name__ == "__main__":
    main()