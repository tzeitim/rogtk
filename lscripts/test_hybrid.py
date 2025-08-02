#!/usr/bin/env python3
"""
Quick test to verify hybrid segments function is available
"""

try:
    import rogtk
    print("✓ rogtk imported successfully")
    
    # Test if the hybrid function is available
    if hasattr(rogtk, 'bam_to_arrow_ipc_htslib_hybrid_segments'):
        print("✓ bam_to_arrow_ipc_htslib_hybrid_segments function is available")
        func = rogtk.bam_to_arrow_ipc_htslib_hybrid_segments
        print(f"✓ Function object: {func}")
        print("\n🎉 Hybrid segments implementation ready for testing!")
        print("\nTo test performance:")
        print("  python hybrid_segments_benchmark.py --bam-path /path/to/your/file.bam")
    else:
        print("❌ bam_to_arrow_ipc_htslib_hybrid_segments function not found")
        print("Available functions:")
        for attr in sorted(dir(rogtk)):
            if 'bam_to_arrow' in attr:
                print(f"  - {attr}")
                
except ImportError as e:
    print(f"❌ Failed to import rogtk: {e}")
except Exception as e:
    print(f"❌ Error: {e}")