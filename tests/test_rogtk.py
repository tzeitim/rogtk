import pytest
from rogtk import (
    sum_as_string,  # Your basic example function
    parse_paired_fastqs,  # Your FASTQ parsing function
    merge_paired_fastqs,
    fastq_to_parquet,
    fracture_sequences,
)

def test_sum_as_string():
    assert sum_as_string(2, 3) == "5"
    assert sum_as_string(0, 0) == "0"
    assert sum_as_string(10, 20) == "30"

def test_cigar_basic_import():
    """Test that the CIGAR parsing functionality is importable"""
    from rogtk import oparse_cigar
    assert callable(oparse_cigar)

# Add placeholder tests for your main functions
# These will need real test data to be properly implemented
def test_parse_paired_fastqs_arguments():
    """Test that parse_paired_fastqs accepts the correct arguments"""
    with pytest.raises(TypeError):
        # Should fail because required arguments are missing
        parse_paired_fastqs()

def test_merge_paired_fastqs_arguments():
    """Test that merge_paired_fastqs accepts the correct arguments"""
    with pytest.raises(TypeError):
        # Should fail because required arguments are missing
        merge_paired_fastqs()

def test_fracture_sequences_arguments():
    """Test that fracture_sequences accepts the correct arguments"""
    with pytest.raises(TypeError):
        # Should fail because required arguments are missing
        fracture_sequences()
