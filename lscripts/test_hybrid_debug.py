#!/usr/bin/env python3
"""
Debug version of hybrid test with detailed logging
"""

import rogtk
import time
import sys

def test_hybrid_debug(bam_path):
    print("🔍 HYBRID SEGMENTS DEBUG TEST")
    print("=" * 50)
    
    # Test with small limit first
    limit = 100_000  # 100k records for quick debugging
    
    print(f"Testing with {limit:,} records...")
    print(f"BAM file: {bam_path}")
    print()
    
    try:
        start = time.time()
        
        rogtk.bam_to_arrow_ipc_htslib_hybrid_segments(
            bam_path=bam_path,
            arrow_ipc_path="debug_output.arrow",
            include_sequence=True,
            include_quality=True,
            limit=limit,
            num_segments=2,
            batch_size=10_000,  # Smaller batches for debugging
            max_bgzf_threads=2,
            writing_threads=2
        )
        
        duration = time.time() - start
        throughput = limit / duration
        
        print(f"\n📊 RESULTS:")
        print(f"  Duration: {duration:.1f}s")
        print(f"  Throughput: {throughput:,.0f} records/sec")
        print(f"  Expected (2 segments): ~{153_000 * 2:,} rec/sec")
        print(f"  Efficiency: {(throughput / (153_000 * 2)) * 100:.1f}%")
        
        # Check if output file exists
        import os
        if os.path.exists("debug_output.arrow"):
            size_mb = os.path.getsize("debug_output.arrow") / (1024 * 1024)
            print(f"  Output file: {size_mb:.1f} MB")
            os.remove("debug_output.arrow")  # Cleanup
        else:
            print("  ⚠️  No output file created!")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        bam_path = sys.argv[1]
    else:
        bam_path = "/home/projects/nyosef/pedro/projects/lt/datain/20250716/fracture/PBC15320_20250123_jurkat_1_1.bam"
    
    test_hybrid_debug(bam_path)