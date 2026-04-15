# Installation

## From PyPI

```bash
pip install molcrafts-molpack
```

The PyPI package `molcrafts-molpack` installs a Python module named
`molpack`:

```python
import molpack
```

Pre-built wheels are published for CPython 3.9–3.12 on Linux (x86_64,
aarch64), macOS (x86_64, arm64), and Windows (x86_64).

## Optional: `molcrafts-molrs`

`molpack` does not parse file formats. The canonical companion for PDB
/ XYZ / LAMMPS / SMILES is `molcrafts-molrs`:

```bash
pip install molcrafts-molrs
```

The two packages are independent — `molpack` operates on numpy arrays
only. Installing `molrs` is recommended but not required.

## Building from source

Requires a Rust toolchain (1.85+) and `maturin`:

```bash
git clone https://github.com/MolCrafts/molpack
cd molpack/python
pip install maturin
maturin develop --release
```

This builds against the local `molpack` Rust crate under `../`.

## Verification

```python
import molpack
import numpy as np

target = molpack.Target.from_coords(
    np.zeros((1, 3)),
    np.array([1.0]),
    count=1,
)
print(target)  # Target(natoms=1, count=1, name=None)
```
