from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from polars.plugins import register_plugin_function
from polars.type_aliases import IntoExpr
from .utils import *

# Try to import all functions including HTSlib if available
try:
    from rogtk.rogtk import (
        sum_as_string,
        oparse_cigar,
        merge_paired_fastqs,
        parse_paired_fastqs,
        fastq_to_parquet,
        fracture_fasta,
        fracture_sequences,
        bam_to_parquet,
        bam_to_arrow_ipc,
        bam_to_arrow_ipc_parallel,
        bam_to_arrow_ipc_gzp_parallel,
        bam_to_arrow_ipc_htslib_parallel,
        bam_to_arrow_ipc_htslib_multi_reader_parallel,
        bam_to_arrow_ipc_htslib_optimized,
        bam_to_arrow_ipc_htslib_mmap_parallel,
        bam_to_arrow_ipc_htslib_bgzf_blocks,
        bam_to_arrow_ipc_htslib_hybrid_segments,
        bam_to_arrow_ipc_htslib_hybrid_optimized,
        bam_to_arrow_ipc_htslib_hybrid_minimal_fix,
    )
    _HTSLIB_AVAILABLE = True
except ImportError:
    # Fallback: import without HTSlib function
    from rogtk.rogtk import (
        sum_as_string,
        oparse_cigar,
        merge_paired_fastqs,
        parse_paired_fastqs,
        fastq_to_parquet,
        fracture_fasta,
        fracture_sequences,
        bam_to_parquet,
        bam_to_arrow_ipc,
        bam_to_arrow_ipc_parallel,
        bam_to_arrow_ipc_gzp_parallel,
    )
    _HTSLIB_AVAILABLE = False
    bam_to_arrow_ipc_htslib_parallel = None
    bam_to_arrow_ipc_htslib_multi_reader_parallel = None
    bam_to_arrow_ipc_htslib_optimized = None
    bam_to_arrow_ipc_htslib_mmap_parallel = None

@pl.api.register_expr_namespace("dna")
class DnaNamespace:
    def __init__(self, expr: pl.Expr):
        self._expr = expr
    
    def reverse_complement(self) -> pl.Expr:
        """Generate reverse complement of DNA sequences."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="reverse_complement_series",
            args=self._expr,
            is_elementwise=True,
        )

#@pl.api.register_expr_namespace("cigar")
def parse_cigar(expr: IntoExpr, block_dels: bool=False) -> pl.Expr:
    """mr cigar"""
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="parse_cigar_series",
        args=expr,
        kwargs={"block_dels": block_dels},#ParseCigarKwargs
        is_elementwise=True,
    )
#@pl.api.register_expr_namespace("phred")
def phred_to_numeric_str(expr: IntoExpr, base: int=33) -> pl.Expr:
    """Convert a PHRED score string into it's numeric value"""
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="phred_to_numeric_series_str",
        args=expr,
        kwargs={"base": base},
        is_elementwise=True,
    )

#@pl.api.register_expr_namespace("phred")
def nn(expr: IntoExpr, base: int=33) -> pl.Expr:
    """Convert a PHRED score string into it's numeric value"""
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="nn",
        args=expr,
        kwargs={"base": base},
        is_elementwise=True,
    )

#@pl.api.register_expr_namespace("fracture")
def assemble_sequences(
    expr: IntoExpr,
    k: int = 10,
    min_coverage: int = 5,
    method: str = 'shortest_path',  
    start_anchor: str | None = None,
    end_anchor: str | None = None, 
    min_length: int | None = None,
    export_graphs: bool = False,
    only_largest: bool = False,
    auto_k: bool = False,
    prefix: str | None = None
) -> pl.Expr:
    """
    Assemble DNA sequences using a de Bruijn graph approach.
    
    Parameters
    ----------
    expr : IntoExpr
        Input expression containing DNA sequences
    k : int
        K-mer size for graph construction (used if auto_k=False)
    min_coverage : int
        Minimum k-mer coverage threshold
    method : str
        Assembly method ('compression' or 'shortest_path')
    start_anchor : str, optional
        Start sequence anchor for shortest_path method
    end_anchor : str, optional 
        End sequence anchor for shortest_path method
    export_graphs : bool, optional
        Whether to export graph visualization files
    ...
    """
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="assemble_sequences_expr",
        args=expr,
        kwargs={
            "k": k,
            "min_coverage": min_coverage,
            "method": method,
            "start_anchor": start_anchor,
            "end_anchor": end_anchor,
            "min_length": min_length,
            "export_graphs": export_graphs,
            "only_largest": only_largest,
            "auto_k": auto_k,
            "prefix": prefix
        },
        returns_scalar=True,
        is_elementwise=False,
    )

def assemble_sequences_with_anchors(
    expr: IntoExpr,
    start_anchor_col: IntoExpr,
    end_anchor_col: IntoExpr,
    k: int = 17,
    min_coverage: int = 25,
    method: str = 'shortest_path',
    min_length: int | None = None,
    export_graphs: bool = False,
    auto_k: bool = False,
    prefix: str | None = None
) -> pl.Expr:
    """
    Assemble DNA sequences with dynamic per-group anchor sequences from columns.

    This function enables single group_by assembly operations where each group
    can have different start/end anchors, avoiding slow Python loops over segment types.

    Parameters
    ----------
    expr : IntoExpr
        Input expression containing DNA sequences to assemble
    start_anchor_col : IntoExpr
        Column expression containing start anchor sequences (first value per group used)
    end_anchor_col : IntoExpr
        Column expression containing end anchor sequences (first value per group used)
    k : int, default 17
        K-mer size for graph construction
    min_coverage : int, default 25
        Minimum k-mer coverage threshold
    method : str, default 'shortest_path'
        Assembly method (only 'shortest_path' supported with dynamic anchors)
    min_length : int, optional
        Minimum contig length to return
    export_graphs : bool, default False
        Whether to export graph visualization files
    auto_k : bool, default False
        Whether to automatically adjust k-mer size
    prefix : str, optional
        Prefix for output files when export_graphs is True

    Returns
    -------
    pl.Expr
        Expression returning assembled consensus sequence(s)

    Examples
    --------
    >>> segments.group_by(['umi', 'start_meta', 'end_meta']).agg(
    ...     rogtk.assemble_sequences_with_anchors(
    ...         expr=pl.col('segment_seq'),
    ...         start_anchor_col=pl.first('start_meta_seq'),
    ...         end_anchor_col=pl.first('end_meta_seq'),
    ...         k=15,
    ...         min_coverage=20,
    ...     ).alias('consensus_seq')
    ... )
    """
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="assemble_sequences_with_anchors_expr",
        args=[expr, start_anchor_col, end_anchor_col],
        kwargs={
            "k": k,
            "min_coverage": min_coverage,
            "method": method,
            "start_anchor": None,  # Not used - anchors come from columns
            "end_anchor": None,    # Not used - anchors come from columns
            "min_length": min_length,
            "export_graphs": export_graphs,
            "only_largest": False,
            "auto_k": auto_k,
            "prefix": prefix
        },
        returns_scalar=True,
        is_elementwise=False,
    )

def sweep_assembly_params(
    expr: IntoExpr,
    k_start: int = 5,
    k_end: int = 32,
    k_step: int = 1,
    cov_start: int = 1,
    cov_end: int = 150,
    cov_step: int = 1,
    method: str = 'shortest_path',
    start_anchor: str | None = None,
    end_anchor: str | None = None,
    min_length: int | None = None,
    export_graphs: bool = False,
    prefix: str | None = None,
    auto_k: bool = False,
) -> pl.Expr:
    """
    Run sequence assembly across ranges of k-mer size and minimum coverage parameters.
    
    Parameters
    ----------
    ...
    method : str
        Assembly method ('compression' or 'shortest_path')
    start_anchor : str, optional
        Start sequence anchor for shortest_path method
    end_anchor : str, optional
        End sequence anchor for shortest_path method
    ...
    """
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="sweep_assembly_params_expr",
        args=expr,
        kwargs={
            "k_start": k_start,
            "k_end": k_end,
            "k_step": k_step,
            "cov_start": cov_start,
            "cov_end": cov_end,
            "cov_step": cov_step,
            "method": method,
            "start_anchor": start_anchor,
            "end_anchor": end_anchor,
            "min_length": min_length,
            "export_graphs": export_graphs,
            "prefix": prefix,
            "auto_k": auto_k
        },
        returns_scalar=True,
        is_elementwise=False,
    )

def optimize_assembly(
    expr: IntoExpr,
    method: str = 'shortest_path',
    start_anchor: str | None = None,
    end_anchor: str | None = None,
    start_k: int = 31,
    start_min_coverage: int = 1,
    min_length: int | None = None,
    export_graphs: bool = False,
    prefix: str | None = None,
    max_iterations: int | None = None,
    explore_k: bool | None = None,
    prioritize_length: bool | None = None,
) -> pl.Expr:
    if start_anchor is None or end_anchor is None:
        raise ValueError("Both start_anchor and e®d_anchor are required")
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="optimize_assembly_expr",
        args=expr,
        kwargs={
            "method": method,
            "start_anchor": start_anchor,
            "end_anchor": end_anchor,
            "start_k": start_k,
            "start_min_coverage": start_min_coverage,
            "min_length": min_length,
            "export_graphs": export_graphs,
            "prefix": prefix,
            "max_iterations": max_iterations,
            "explore_k": explore_k,
            "prioritize_length": prioritize_length
        },
        returns_scalar=True,
        is_elementwise=False,
    )

@pl.api.register_expr_namespace("hamming")
class HammingExpr:
    def __init__(self, expr: pl.Expr):
        self._expr = expr
    
    def distance(self, target: str) -> pl.Expr:
        """Calculate Hamming distance to target sequence."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="hamming_distance_expr",
            args=self._expr,
            kwargs={"target": target},
            is_elementwise=True,
        )
    
    def within(self, target: str, max_distance: int = 1) -> pl.Expr:
        """Check if sequence is within Hamming distance of target."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="hamming_within_expr", 
            args=self._expr,
            kwargs={"target": target, "max_distance": max_distance},
            is_elementwise=True,
        )

@pl.api.register_expr_namespace("fuzzy")
class FuzzyExpr:
    def __init__(self, expr: pl.Expr):
        self._expr = expr
    
    # Option 1: Pre-generated patterns (use with fuzzy_match_str)
    def replace(self, pattern: str, replacement: str, literal: bool = False) -> pl.Expr:
        """Replace fuzzy pattern matches in sequences (use with fuzzy_match_str)."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="fuzzy_replace_expr",
            args=self._expr,
            kwargs={"pattern": pattern, "replacement": replacement, "literal": literal},
            is_elementwise=True,
        )
    
    def contains(self, pattern: str, literal: bool = False) -> pl.Expr:
        """Check if sequence contains fuzzy pattern (use with fuzzy_match_str)."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="fuzzy_contains_expr",
            args=self._expr,
            kwargs={"pattern": pattern, "literal": literal},
            is_elementwise=True,
        )
    
    # Option 2: Native Rust pattern generation
    def match(self, target: str, wildcard: str = ".{0,1}", include_original: bool = True, max_length: int = 100) -> pl.Expr:
        """Check if sequence fuzzy matches target (pattern generated in Rust)."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="fuzzy_contains_native_expr",
            args=self._expr,
            kwargs={
                "target": target, 
                "wildcard": wildcard, 
                "include_original": include_original,
                "max_length": max_length
            },
            is_elementwise=True,
        )

    def replace_target(self, target: str, replacement: str, 
                      wildcard: str = ".{0,1}", include_original: bool = True, 
                      max_length: int = 100, replace_all: bool = False) -> pl.Expr:
        """Replace fuzzy matches of target (pattern generated in Rust)."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="fuzzy_replace_native_expr",
            args=self._expr,
            kwargs={
                "target": target,
                "replacement": replacement,
                "wildcard": wildcard, 
                "include_original": include_original,
                "max_length": max_length,
                "replace_all": replace_all
            },
            is_elementwise=True,
        )

@pl.api.register_expr_namespace("umi")
class UmiNamespace:
    def __init__(self, expr: pl.Expr):
        self._expr = expr
    
    def complexity_all(self) -> pl.Expr:
        """Calculate all UMI complexity scores and return as struct."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_complexity_all_expr",
            args=self._expr,
            is_elementwise=True,
        )
    
    def all_scores(self) -> pl.Expr:
        """Calculate all UMI complexity scores and return as struct (alias for complexity_all)."""
        return self.complexity_all()
    
    def shannon_entropy(self) -> pl.Expr:
        """Calculate Shannon entropy of UMI sequence."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_shannon_entropy_expr",
            args=self._expr,
            is_elementwise=True,
        )
    
    def linguistic_complexity(self) -> pl.Expr:
        """Calculate linguistic complexity (unique k-mers ratio) of UMI sequence."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_linguistic_complexity_expr",
            args=self._expr,
            is_elementwise=True,
        )
    
    def homopolymer_fraction(self) -> pl.Expr:
        """Calculate fraction of UMI sequence in homopolymer runs (3+ identical bases)."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_homopolymer_fraction_expr",
            args=self._expr,
            is_elementwise=True,
        )
    
    def dinucleotide_entropy(self) -> pl.Expr:
        """Calculate dinucleotide entropy of UMI sequence."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_dinucleotide_entropy_expr",
            args=self._expr,
            is_elementwise=True,
        )
    
    def combined_score(self) -> pl.Expr:
        """Calculate combined complexity score of UMI sequence."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_combined_score_expr",
            args=self._expr,
            is_elementwise=True,
        )
    
    def longest_homopolymer_run(self) -> pl.Expr:
        """Calculate longest homopolymer run length in UMI sequence."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_longest_homopolymer_expr",
            args=self._expr,
            is_elementwise=True,
        )
    
    def dust_score(self) -> pl.Expr:
        """Calculate DUST score (low-complexity measure) for UMI sequence."""
        return register_plugin_function(
            plugin_path=Path(__file__).parent,
            function_name="umi_dust_score_expr",
            args=self._expr,
            is_elementwise=True,
        )

def umi_complexity_scores(expr: IntoExpr) -> pl.Expr:
    """
    Calculate all UMI complexity scores and return them as separate columns.
    
    This function returns a struct containing:
    - shannon_entropy: Shannon entropy of nucleotide frequencies
    - linguistic_complexity: Ratio of unique k-mers to total possible k-mers  
    - homopolymer_fraction: Fraction of sequence in homopolymer runs (3+ bases)
    - dinucleotide_entropy: Entropy of dinucleotide frequencies
    - longest_homopolymer_run: Length of longest homopolymer run
    - dust_score: DUST score for low-complexity detection
    - combined_score: Weighted combination of all metrics
    
    Parameters
    ----------
    expr : IntoExpr
        Input expression containing UMI sequences
        
    Returns
    -------
    pl.Expr
        Struct expression with all complexity scores
        
    Examples
    --------
    >>> df.with_columns(umi_complexity_scores(pl.col("umi")))
    >>> df.unnest("umi_complexity")  # To expand struct into separate columns
    """
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="umi_complexity_all_expr",
        args=expr,
        is_elementwise=True,
    )
