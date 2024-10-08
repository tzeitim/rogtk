import polars as pl
import rogtk

POLARS_VERBOSE=1

df = pl.DataFrame({'dense': [[0, 9], [8, 6, 0, 9], None, [3, 3]]})
df = pl.DataFrame({'a': [['a', 'b'], ['c', 'd', 'e', 'f'], None, ['g', 'h']]})
df = pl.DataFrame({'a': ['!!!!!', 'IIIIIII']})

print(df.with_columns(rogtk.phred_to_numeric_str('a')))
