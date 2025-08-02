#!/usr/bin/env python3
"""
Test hybrid performance at different scales to find the sweet spot
"""

import rogtk
import time
import sys

def test_scale(bam_path, limit, segments=2):
    print(f"🧪 Testing {limit:,} records with {segments} segments...")
    
    try:
        start = time.time()
        
        rogtk.bam_to_arrow_ipc_htslib_hybrid_segments(
            bam_path=bam_path,
            arrow_ipc_path=f"scale_test_{limit}.arrow",
            include_sequence=True,
            include_quality=True,
            limit=limit,
            num_segments=segments,
            batch_size=10_000,
            max_bgzf_threads=4,
            writing_threads=4
        )
        
        duration = time.time() - start
        throughput = limit / duration
        expected = 153_000 * segments
        efficiency = (throughput / expected) * 100
        
        print(f"  ⚡ {throughput:,.0f} rec/sec ({efficiency:.1f}% efficiency, {duration:.1f}s)")
        return throughput
        
    except Exception as e:
        print(f"  ❌ Failed: {e}")
        return 0
    finally:
        # Cleanup
        import os
        output_file = f"scale_test_{limit}.arrow"
        if os.path.exists(output_file):
            os.remove(output_file)

def main():
    if len(sys.argv) > 1:
        bam_path = sys.argv[1]
    else:
        bam_path = "/home/projects/nyosef/pedro/projects/lt/datain/20250716/fracture/PBC15320_20250123_jurkat_1_1.bam"
    
    print("🚀 HYBRID SEGMENTS SCALING TEST")
    print("=" * 50)
    print(f"BAM file: {bam_path}")
    print()
    
    # Test different scales
    scales = [100_000, 500_000, 1_000_000, 2_000_000]
    results = []
    
    for limit in scales:
        throughput = test_scale(bam_path, limit, segments=10)
        results.append((limit, throughput))
        print()
    
    print("📊 SCALING ANALYSIS:")
    print("-" * 50)
    print(f"{'Records':<12} {'Throughput':<15} {'vs 100k':<10}")
    print("-" * 50)
    
    baseline = results[0][1] if results[0][1] > 0 else 1
    
    for limit, throughput in results:
        if throughput > 0:
            improvement = (throughput / baseline - 1) * 100
            print(f"{limit:>10,} {throughput:>12,.0f} rec/sec {improvement:>+6.1f}%")
        else:
            print(f"{limit:>10,} {'FAILED':>15} {'N/A':>10}")
    
    print()
    
    # Find best performance
    best_limit, best_throughput = max(results, key=lambda x: x[1])
    if best_throughput > 0:
        print(f"🎯 Best performance: {best_throughput:,.0f} rec/sec at {best_limit:,} records")
        
        if best_throughput >= 200_000:
            print("🎉 TARGET ACHIEVED! Exceeded 200k records/sec!")
        elif best_throughput >= 150_000:
            print("🔥 Very close to target! Approaching 200k rec/sec")
        elif best_throughput >= 100_000:
            print("✅ Significant improvement! Over 100k rec/sec")
        
        target_progress = (best_throughput / 200_000) * 100
        print(f"📈 Progress toward 200k target: {target_progress:.1f}%")

if __name__ == "__main__":
    main()
