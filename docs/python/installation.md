# Installation

## From PyPI

```bash
pip install molcrafts-molpack
```

The PyPI package `molcrafts-molpack` installs a Python module named `molpack`:

```python
import molpack
```

Pre-built wheels are published for CPython **3.12** and **3.13** on Linux
(manylinux x86-64) and macOS (universal2). Other platforms fall back to the
sdist and need a Rust toolchain to build locally.

## Frame I/O

`molcrafts-molpack` installs `molcrafts-molrs` as a dependency. The
importable `molrs` module provides the frame type plus PDB and XYZ readers:

```bash
pip install molcrafts-molpack
```

```python
import molrs
from molpack import InsideBoxRestraint, Molpack, Target

frame = molrs.io.read_pdb("water.pdb")
water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)
frame = Molpack().with_seed(42).pack([water], max_loops=200)
```

`Target` takes a `molrs.Frame` or `molpy.Frame`, resolved zero-copy through
its FFI capsule.

## Building from source

Requires a Rust toolchain (1.91+) and `maturin`:

```bash
git clone https://github.com/MolCrafts/molpack
cd molpack/python
pip install maturin
maturin develop --release
```

This builds against the local `molpack` Rust crate under `../`.

## Verification

```python
import molrs
from molpack import Target

frame = molrs.Frame.from_dict({
    "blocks": {"atoms": {"x": [0.0], "y": [0.0], "z": [0.0], "element": ["O"]}},
    "meta": {},
})
target = Target(frame, count=1).with_name("mol")
print(target)  # Target(natoms=1, count=1, name=Some("mol"))
```
