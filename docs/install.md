# Install

`molpack` ships three surfaces — the CLI binary, the Rust crate, and the
Python binding. Pick the one that matches your workflow; they share the same
engine and packing model.

## CLI

For Packmol-style `.inp` scripts:

```bash
cargo install molcrafts-molpack --features cli
molpack --help
```

!!! tip "Path resolution"
    In file-arg mode (`molpack job.inp`), relative paths inside the script are
    resolved against the **script directory**. In stdin mode they resolve
    against the current working directory.

## Rust crate

For programmatic use from Rust:

```bash
cargo add molcrafts-molpack
```

Optional features (crate defaults to none enabled):

| Feature | Purpose |
|---|---|
| `cli` | `molpack` binary + clap (implies `io`) |
| `io` | Template reading / output writing through the molrs reader and writer of each file's format (`script::StructureFormat`), and `XYZHandler` |
| `rayon` | Parallel objective evaluation |

A force-field optimizer bound through `with_optimizer` (e.g. molrs's `Lbfgs`)
needs molrs's `ff` feature; enable it on your own `molcrafts-molrs` dependency.

```toml
# Cargo.toml — common combinations
molcrafts-molpack = { version = "0.4", features = ["io", "rayon"] }
```

## Python binding

For notebooks and pipelines (Python 3.12+):

```bash
pip install molcrafts-molpack
```

`molcrafts-molrs` is installed as a dependency and provides `molrs.core.Frame` plus
PDB / XYZ readers. The wheel itself is I/O-free — pass frames in, get frames
out.

```python
import molpack
print(molpack.GenCanPack)
```

!!! note "Pre-built wheels"
    Wheels are published for CPython **3.12** and **3.13** on Linux
    (manylinux x86-64) and macOS (universal2). Other platforms fall back to the
    sdist and need a Rust toolchain to build.

## Build from source

When you are modifying the crate or Python binding, check out **molrs** as a
sibling (path deps resolve `../molrs/molrs`). molpack 0.4 builds on the molrs
**0.16** line — the `v0.16.0` tag or a later 0.16 commit:

```bash
# sibling layout
# workspace/
# ├── molrs/
# └── molpack/

git clone https://github.com/MolCrafts/molpack
cd molpack

# Rust library + CLI
cargo build --features cli

# Python wheel (editable)
cd python && maturin develop --release
```

See [Development](development/) for tests, hooks, and contribution workflow.
