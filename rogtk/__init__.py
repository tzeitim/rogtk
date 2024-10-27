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

#def noop(expr: IntoExpr) -> pl.Expr:
#    expr = parse_into_expr(expr)
#    return expr.register_plugin(
#        lib=lib,
#        symbol="noop",
#        is_elementwise=True,
#    )
#
