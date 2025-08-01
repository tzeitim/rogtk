#!/usr/bin/env python3
"""
Verify Output Size Test 
=======================

Check if the original hybrid implementation is actually outputting
all the records it claims to be processing.
"""

import argparse
import time
import os
import sys
import pyarrow as pa
import pyarrow.ipc as ipc

try:
    import rogtk
except ImportError:
    print("ERROR: rogtk module not found")
    sys.exit(1)


def count_records_in_arrow_file(file_path):
    """Count actual records in Arrow IPC file"""
    try:
        with open(file_path, 'rb') as f:
            reader = ipc.RecordBatchFileReader(f)
            total_records = 0
            for i in range(reader.num_record_batches):
                batch = reader.get_record_batch(i)
                total_records += batch.num_rows
            return total_records
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return 0


def run_test_and_verify(name, function_name, bam_path, limit, **params):
    """Run test and verify actual output"""
    print(f"\n🧪 Testing: {name}")
    print("-" * 50)
    
    output_file = f"verify_{function_name}_{int(time.time())}.arrow.ipc"
    
    try:
        print(f"Running with {limit:,} record limit...", end=" ", flush=True)
        
        start_time = time.time()
        func = getattr(rogtk, function_name)
        func(
            bam_path=bam_path,
            arrow_ipc_path=output_file,
            include_sequence=True,
            include_quality=True,
            limit=limit,
            **params
        )
        end_time = time.time()
        
        duration = end_time - start_time
        reported_throughput = limit / duration
        
        print(f"✓ {duration:.1f}s")
        print(f"Reported throughput: {reported_throughput:,.0f} rec/sec")
        
        # Check actual output
        if os.path.exists(output_file):
            file_size_mb = os.path.getsize(output_file) / (1024 * 1024)
            actual_records = count_records_in_arrow_file(output_file)
            
            print(f"Output file size: {file_size_mb:.1f} MB")
            print(f"Actual records in output: {actual_records:,}")
            print(f"Expected records: {limit:,}")
            
            if actual_records > 0:
                completeness = (actual_records / limit) * 100
                actual_throughput = actual_records / duration
                
                print(f"Completeness: {completeness:.1f}%")
                print(f"ACTUAL throughput: {actual_throughput:,.0f} rec/sec")
                
                if completeness < 50:
                    print("🚨 MAJOR DATA LOSS! Less than 50% of records in output")
                elif completeness < 90:
                    print("⚠️  Significant data loss")
                elif completeness < 99:
                    print("⚠️  Minor data loss")
                else:
                    print("✅ Complete data output")
                    
                return actual_throughput, actual_records
            else:
                print("❌ No records found in output")
                return 0, 0
        else:
            print("❌ No output file created")
            return 0, 0
            
    except Exception as e:
        print(f"❌ Failed: {e}")
        return 0, 0
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def main():
    parser = argparse.ArgumentParser(description="Verify Output Size Test")
    parser.add_argument("--bam-path", required=True, help="Path to BAM file")
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("OUTPUT VERIFICATION TEST")
    print("=" * 70)
    print(f"BAM file: {args.bam_path}")
    print(f"Testing with smaller limit to verify behavior...")
    
    limit = 100_000  # Smaller limit for clearer comparison
    
    optimal_params = {
        'batch_size': 15000,
        'read_buffer_mb': 2048,
        'write_buffer_mb': 128,
        'max_bgzf_threads': 4,
        'writing_threads': 6
    }
    
    # Test 1: Original hybrid (suspected to have incomplete output)
    original_throughput, original_count = run_test_and_verify(
        "Original Hybrid Segments",
        "bam_to_arrow_ipc_htslib_hybrid_segments",
        args.bam_path,
        limit,
        **optimal_params
    )
    
    # Test 2: Single optimized reader (should be complete)
    single_throughput, single_count = run_test_and_verify(
        "Single Optimized Reader",
        "bam_to_arrow_ipc_htslib_optimized",
        args.bam_path,
        limit,
        **optimal_params
    )
    
    # Test 3: Minimal fix (should be complete but slower)
    fixed_throughput, fixed_count = run_test_and_verify(
        "Minimal Fix Hybrid",
        "bam_to_arrow_ipc_htslib_hybrid_minimal_fix",
        args.bam_path,
        limit,
        **optimal_params
    )
    
    print(f"\n" + "=" * 70)
    print("VERIFICATION ANALYSIS")
    print("=" * 70)
    
    print(f"Expected records: {limit:,}")
    print(f"Original hybrid:  {original_count:,} records → {original_throughput:,.0f} rec/sec")
    print(f"Single reader:    {single_count:,} records → {single_throughput:,.0f} rec/sec") 
    print(f"Minimal fix:      {fixed_count:,} records → {fixed_throughput:,.0f} rec/sec")
    
    if original_count < single_count:
        missing_pct = ((single_count - original_count) / single_count) * 100
        print(f"\n🚨 CONFIRMED: Original hybrid missing {missing_pct:.1f}% of data!")
        print(f"   This explains the artificially high throughput numbers")
    
    if fixed_count >= single_count and fixed_throughput > 0:
        print(f"\n✅ Minimal fix provides complete data")
        if fixed_throughput >= original_throughput * 0.5:  # Account for complete data
            print(f"   Performance is reasonable considering complete output")
        else:
            print(f"   Performance needs improvement")
    
    print(f"\n💡 CONCLUSION:")
    if original_count < single_count:
        print(f"   The original 'hybrid' was a false optimization")  
        print(f"   It threw away {((single_count - original_count) / single_count * 100):.0f}% of the data")
        print(f"   Real hybrid performance should be judged on COMPLETE output")


if __name__ == "__main__":
    main()