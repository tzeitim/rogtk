import rogtk as rr
import polars as pl

def phred_to_numeric(df, col_name):
    """
    Transforms a phred score string into a list ints
    """
    return(
            df
            .lazy()
            .with_columns(
                rr.phred_to_numeric_str(col_name)
                  .str.split('|')
                  .list.eval(pl.element().cast(pl.UInt8)))
            .collect()
        )

