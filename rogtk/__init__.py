from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from polars.plugins import register_plugin_function
from polars.type_aliases import IntoExpr
from .utils import *

from rogtk.rogtk import (
    sum_as_string,
    oparse_cigar,
    merge_paired_fastqs,
    parse_paired_fastqs,
    fastq_to_parquet,
    fracture_fasta,
    fracture_sequences,
    bam_to_parquet
)

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
    
    def replace_target(self, target: str, replacement: str = "[FUZZY_MATCH]", wildcard: str = ".{0,1}", include_original: bool = True, max_length: int = 100) -> pl.Expr:
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
                "max_length": max_length
            },
            is_elementwise=True,
        )
