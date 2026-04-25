# Install

`molpack` is one engine with three entry points. The correct installation path
depends on how you want to describe packing jobs, not on three different
feature sets. This chapter explains which surface to choose and how to verify
that it works.

## What You'll Learn

After this chapter, you should be able to:

- choose between the CLI, Rust crate, and Python binding;
- install the surface that matches your workflow;
- verify that the installed toolchain is usable before writing real jobs.

## 1. Install the Command-Line Interface

Choose the CLI when your workflow already centers on Packmol-style `.inp`
scripts and structure files.

```bash
cargo install molcrafts-molpack --features cli
```

Verification:

```bash
molpack --help
```

This path preserves the classic file-driven flow while extending supported input
formats beyond Packmol's traditional PDB and XYZ coverage.

## 2. Install the Rust Crate

Choose the Rust crate when packing is part of a compiled program or when you
need typed access to targets, handlers, regions, or custom restraints.

```bash
cargo add molcrafts-molpack
```

If you need the binary in the same project, enable the `cli` feature in your
workspace or install it separately through `cargo install`.

## 3. Install the Python Binding

Choose Python when the packing step belongs in a notebook, data pipeline, or
general scripting environment.

```bash
pip install molcrafts-molpack
```

Most real workflows should also install `molcrafts-molrs`, because that gives
you `molrs.Frame` values and reader/writer helpers for PDB and XYZ:

```bash
pip install molcrafts-molpack molcrafts-molrs
```

Verification:

```python
import molpack
print(molpack.Molpack)
```

Pre-built wheels are published for CPython 3.12 and 3.13 on Linux, macOS, and
Windows targets covered by the release pipeline.

## 4. Build from Source

Build from source when you are modifying the crate itself, changing the Python
binding, or validating unpublished local changes.

```bash
git clone https://github.com/MolCrafts/molpack
cd molpack
```

Rust library and CLI:

```bash
cargo build
cargo build --features cli --bin molpack
```

Python binding:

```bash
cd python
pip install maturin
maturin develop --release
```

The Python package builds against the local Rust crate under `../`, so it is
the right path for contributor workflows and binding changes.

## 5. Choose the Right Path

If you are unsure, use this rule of thumb:

- start with the CLI if you are replacing Packmol scripts;
- start with Rust if you need new geometry or integration into a compiled tool;
- start with Python if you want notebook- or pipeline-friendly orchestration
  around existing frames.

## Takeaway

Installation in `molpack` is not about unlocking features. It is about choosing
the surface that fits your workflow. Once the tool is installed and verified,
the next useful step is a first complete packing run in
[Getting Started](getting_started.md).
