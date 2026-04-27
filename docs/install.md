# Install

`molpack` ships in three forms — the CLI binary, the Rust crate, and the
Python binding. Pick the surface that matches your workflow.

## CLI

For Packmol-style `.inp` scripts:

```bash
cargo install molcrafts-molpack --features cli
molpack --help
```

## Rust crate

For programmatic use from Rust code:

```bash
cargo add molcrafts-molpack
```

Add the `cli` feature in your workspace if you also want the binary
in-tree.

## Python binding

For notebooks and pipelines:

```bash
pip install molcrafts-molpack molcrafts-molrs
```

`molcrafts-molrs` provides `molrs.Frame` plus PDB / XYZ readers — the
Python wheel itself is I/O-free.

```python
import molpack
print(molpack.Molpack)
```

Pre-built wheels are published for CPython 3.12 and 3.13 on Linux,
macOS, and Windows.

## Build from source

When you are modifying the crate or Python binding:

```bash
git clone https://github.com/MolCrafts/molpack
cd molpack

# Rust library + CLI binary
cargo build
cargo build --features cli --bin molpack

# Python binding
cd python
pip install maturin
maturin develop --release
```

The Python package builds against the local Rust crate under `../`, so
this is the right path for contributor workflows.

## Next

After verifying the install, run a first pack from
[Getting Started](getting_started.md).
