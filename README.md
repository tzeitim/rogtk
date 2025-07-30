rust ogtk

## Installation

**Standard installation:**
```bash
pip install -e src/rogtk/
# or
maturin develop
# HTSLIB support
LIBCLANG_PATH=$CONDA_PREFIX/lib maturin develop --features htslib --release
```

**With HTSlib support (for parallel BAM processing):**
```bash
# Requires libclang
mamba install libclang
# Build with HTSlib feature
LIBCLANG_PATH=$CONDA_PREFIX/lib maturin develop --features htslib --release
```


Example
```
import rogtk.rogtk as rogtk
rogtk.sum_as_string(1, 2)
```

Useful references

[setuptools](https://setuptools-rust.readthedocs.io/en/latest/reference.html) 
Here one can control 'release' or 'debug' options


- [bio::alignment](https://docs.rs/bio/0.24.0/bio/alignment/pairwise/index.html)
- [seal alignment](https://github.com/regexident/rust-seal)


Notes:

# compilations not updating

```
The issue was that even after running maturin develop (or build with or
without --release) there were some .so files for the python package that were
not overwritten. That's why starting from scratch worked (once) since those
original files were gone. After I delete it (together with the __pycache__ dir)
the python side is up to date. 

```
