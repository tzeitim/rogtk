from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from polars.plugins import register_plugin_function
from polars.type_aliases import IntoExpr

@pl.api.register_expr_namespace("cigar")
def parse_cigar(expr: IntoExpr) -> pl.Expr:
    """mr cigar"""
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="parse_cigar_series",
        args=expr,
        is_elementwise=True,
    )

#@pl.api.register_expr_namespace("language")
def pig_latinnify(expr: IntoExpr) -> pl.Expr:
    """Pig-latinnify expression."""
    return register_plugin_function(
        plugin_path=Path(__file__).parent,
        function_name="pig_latinnify",
        args=expr,
        is_elementwise=True,
    )
