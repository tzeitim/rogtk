#!/usr/bin/env python3
"""
Test hybrid with more segments to push toward 200k target
"""

import rogtk
import time
import sys

def test_segments(bam_path, limit, segments):
    print(f"🧪 Testing {limit:,} records with {segments} segments...")
    
    try:
        start = time.time()
        
        rogtk.bam_to_arrow_ipc_htslib_hybrid_segments(
            bam_path=bam_path,
            arrow_ipc_path=f"segments_test_{segments}.arrow",
            include_sequence=True,
            include_quality=True,
            limit=limit,
            num_segments=segments,
            batch_size=20_000,
            max_bgzf_threads=4,
            writing_threads=8
        )
        
        duration = time.time() - start
        throughput = limit / duration
        expected = 153_000 * segments
        efficiency = (throughput / expected) * 100
        per_segment = throughput / segments
        
        print(f"  ⚡ {throughput:,.0f} rec/sec total ({per_segment:,.0f} per segment)")
        print(f"  📊 {efficiency:.1f}% efficiency vs theoretical max")
        print(f"  ⏱️  {duration:.1f}s duration")
        
        return throughput
        
    except Exception as e:
        print(f"  ❌ Failed: {e}")
        return 0
    finally:
        # Cleanup
        import os
        output_file = f"segments_test_{segments}.arrow"
        if os.path.exists(output_file):
            os.remove(output_file)

def main():
    if len(sys.argv) > 1:
        bam_path = sys.argv[1]
    else:
        bam_path = "/home/projects/nyosef/pedro/projects/lt/datain/20250716/fracture/PBC15320_20250123_jurkat_1_1.bam"
    
    print("🚀 MULTI-SEGMENT SCALING TEST")
    print("=" * 50)
    print(f"BAM file: {bam_path}")
    print()
    
    limit = 2_000_000  # Use best performing scale from previous test
    print(f"Testing {limit:,} records with different segment counts:")
    print()
    
    segment_counts = [1, 2, 4, 8]
    results = []
    
    for segments in segment_counts:
        throughput = test_segments(bam_path, limit, segments)
        results.append((segments, throughput))
        print()
    
    print("📊 SEGMENT SCALING ANALYSIS:")
    print("-" * 60)
    print(f"{'Segments':<10} {'Throughput':<15} {'Per Segment':<15} {'Efficiency':<12}")
    print("-" * 60)
    
    for segments, throughput in results:
        if throughput > 0:
            per_segment = throughput / segments
            theoretical_max = 153_000 * segments
            efficiency = (throughput / theoretical_max) * 100
            
            print(f"{segments:<10} {throughput:>12,.0f} rec/sec {per_segment:>12,.0f} rec/sec {efficiency:>9.1f}%")
        else:
            print(f"{segments:<10} {'FAILED':>15} {'N/A':>15} {'N/A':>12}")
    
    print()
    
    # Analysis
    best_segments, best_throughput = max(results, key=lambda x: x[1])
    if best_throughput > 0:
        print(f"🎯 Best configuration: {best_segments} segments = {best_throughput:,.0f} rec/sec")
        
        target_progress = (best_throughput / 200_000) * 100
        print(f"📈 Progress toward 200k target: {target_progress:.1f}%")
        
        if best_throughput >= 200_000:
            print("🎉 TARGET ACHIEVED! Exceeded 200k records/sec!")
        elif best_throughput >= 150_000:
            print("🔥 Very close! Within 25% of 200k target")
        elif best_throughput >= 120_000:
            print("✅ Strong progress! Over 120k rec/sec")
        
        # Linear scaling analysis
        print(f"\n💡 Scaling Insights:")
        single_segment = next((t for s, t in results if s == 1), 0)
        if single_segment > 0:
            for segments, throughput in results:
                if segments > 1 and throughput > 0:
                    ideal_scaling = single_segment * segments
                    actual_efficiency = (throughput / ideal_scaling) * 100
                    print(f"   {segments} segments: {actual_efficiency:.1f}% of ideal linear scaling")

if __name__ == "__main__":
    main()