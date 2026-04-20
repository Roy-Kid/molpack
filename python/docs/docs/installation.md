# Installation

## From PyPI

```bash
pip install molcrafts-molpack
```

The PyPI package `molcrafts-molpack` installs a Python module named `molpack`:

```python
import molpack
```

Pre-built wheels are published for CPython 3.12–3.13 on Linux (x86_64,
aarch64), macOS (x86_64, arm64), and Windows (x86_64).

## With `molcrafts-molrs` (recommended)

`molpack` accepts `molrs.Frame` objects directly — the recommended way
to load PDB and XYZ files:

```bash
pip install molcrafts-molpack molcrafts-molrs
```

```python
import molrs
from molpack import InsideBoxRestraint, Molpack, Target

frame = molrs.read_pdb("water.pdb")
water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)
result = Molpack().with_seed(42).pack([water], max_loops=200)
```

`molpack` also accepts plain Python dicts-of-dicts, so `molcrafts-molrs`
is not strictly required.

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
from molpack import Molpack, Target

frame = {
    "atoms": {
        "x": [0.0], "y": [0.0], "z": [0.0],
        "element": ["O"],
    }
}
target = Target(frame, count=1).with_name("mol")
print(target)  # Target(natoms=1, count=1, name=Some("mol"))
```
