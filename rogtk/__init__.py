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
    method: str = 'compression',  # Changed from enum to str
    start_anchor: str | None = None,  # New parameter
    end_anchor: str | None = None,    # New parameter 
    export_graphs: bool = False,
    only_largest: bool = False,
    min_length: int | None = None,
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
    method: str = 'compression',  # Changed from enum to str
    start_anchor: str | None = None,  # New parameter
    end_anchor: str | None = None,    # New parameter
    export_graphs: bool = False,
    prefix: str | None = None
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
    explore_k: bool | None = None,
    prioritize_length: bool | None = None,
) -> pl.Expr:
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
            "explore_k": explore_k,
            "prioritize_length": prioritize_length
        },
        returns_scalar=True,
        is_elementwise=False,
    )
