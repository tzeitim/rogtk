from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from polars.plugins import register_plugin_function
from polars.type_aliases import IntoExpr
from .utils import *

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
    k: int = 31,
    min_coverage: int = 1,
    export_graphs: bool = False,
    only_largest: bool = False,
    min_length: int | None = None,
    auto_k: bool = True,
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
    export_graphs : bool, optional
        Whether to export graph visualization files
    only_largest : bool, optional
        Return only the largest contig
    min_length : int, optional
        Minimum contig length to return
    auto_k : bool, optional
        Automatically estimate optimal k-mer size
    prefix : str, optional
        Prefix for output files
        
    Returns
    -------
    pl.Expr
        Expression containing assembled contigs
    """
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="assemble_sequences_expr",
        args=expr,
        kwargs={
            "k": k,
            "min_coverage": min_coverage,
            "export_graphs": export_graphs,
            "only_largest": only_largest,
            "min_length": min_length,
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
    export_graphs: bool = False,
    prefix: str | None = None
) -> pl.Expr:
    """
    Run sequence assembly across ranges of k-mer size and minimum coverage parameters.
    
    Parameters
    ----------
    expr : IntoExpr
        Input expression containing DNA sequences
    k_start : int
        Starting k-mer size
    k_end : int
        Ending k-mer size (inclusive)
    k_step : int 
        Step size for k-mer parameter
    cov_start : int
        Starting minimum coverage
    cov_end : int
        Ending minimum coverage (inclusive)
    cov_step : int
        Step size for coverage parameter
    export_graphs : bool
        Whether to export assembly graphs
    prefix : str, optional
        Prefix for output files
        
    Returns
    -------
    pl.Expr
        Expression containing list of [k, min_coverage, max_contig_length] for each parameter combination
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
            "export_graphs": export_graphs,
            "prefix": prefix
        },
        returns_scalar=True,
        is_elementwise=False,
    )


def optimize_assembly(
    expr: IntoExpr,
    start_k: int,
    start_min_coverage: int,
    start_anchor: str,
    end_anchor: str,
    max_iterations: int | None = None,
    explore_k: bool | None = None
) -> pl.Expr:
    """
    Optimize assembly parameters by exploring parameter space until anchors are found.
    
    Parameters
    ----------
    expr : IntoExpr
        Input expression containing sequences to assemble
    start_k : int
        Initial k-mer size
    start_min_coverage : int
        Initial minimum coverage threshold
    start_anchor : str
        Sequence that should appear at start of contig
    end_anchor : str
        Sequence that should appear at end of contig
    max_iterations : int, optional
        Maximum optimization iterations, defaults to 50
    explore_k : bool, optional
        Whether to explore k-mer size variations, defaults to False
        
    Returns
    -------
    pl.Expr
        Expression containing optimized assembly result with fields:
        - contig: str
        - k: int
        - min_coverage: int
        - length: int
    """
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="optimize_assembly_expr",
        args=expr,
        kwargs={
            "start_k": start_k,
            "start_min_coverage": start_min_coverage,
            "start_anchor": start_anchor,
            "end_anchor": end_anchor,
            "max_iterations": max_iterations,
            "explore_k": explore_k
        },
        returns_scalar=True,
        is_elementwise=False,
    )
