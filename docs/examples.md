# Examples

The examples in this repository are organized around complete workloads rather
than isolated API fragments. That choice is deliberate: packing only becomes
intelligible when target counts, restraint geometry, and iteration budgets are
studied together.

## What You'll Learn

After this chapter, you should be able to:

- choose the example that matches the physical setup you care about;
- run the Rust or Python version of a canonical workload;
- use examples as migration guides from Packmol input decks.

## 1. Canonical Workloads

The repository carries the five canonical Packmol-style workloads in both
Rust-oriented and Python-oriented forms where that makes sense:

| Workload | Rust example or input | Python example | What it teaches |
|---|---|---|---|
| Water cube | `examples/pack_mixture/` is the closest minimal Rust entry point; the crate docs also use a one-species box | `python/examples/pack_water_cube.py` | the smallest useful pack job |
| Mixture | `examples/pack_mixture` | `python/examples/pack_mixture.py` | two species sharing one restraint volume |
| Bilayer | `examples/pack_bilayer` | `python/examples/pack_bilayer.py` | atom-subset restraints for orienting molecules |
| Interface | `examples/pack_interface` | `python/examples/pack_interface.py` | fixed reference structures plus mobile solvents |
| Spherical | `examples/pack_spherical` | `python/examples/pack_spherical.py` | nested regions and shell-like packing |
| Solvated protein | `examples/pack_solvprotein` | `python/examples/pack_solvprotein.py` | heterogeneous systems with one fixed macromolecule |

The Python examples are regression-aligned with the Rust engine, so they are a
good way to learn the binding without learning a second algorithm.

## 2. Run the Rust Examples

The Rust examples mirror the canonical Packmol workloads and compile against the
crate in this repository:

```bash
cargo run --release --example pack_mixture
cargo run --release --example pack_bilayer
cargo run --release --example pack_interface
cargo run --release --example pack_spherical
cargo run --release --example pack_solvprotein
```

The CLI-facing `.inp` scripts are also kept in the example directories, so you
can compare the programmatic and script-driven representations directly.

## 3. Run the Python Examples

The Python examples live under `python/examples/`:

```bash
cd python
pip install -e .
python examples/pack_water_cube.py
python examples/pack_mixture.py
```

Any example that reads real structure files also needs `molcrafts-molrs`.

## 4. Start with the Right Example

Pick the workload that exercises the constraint you care about:

- choose the water cube when you only need to learn the shape of the API;
- choose the mixture when you care about multi-species packing;
- choose the bilayer when you need atom-subset restraints;
- choose the interface or solvated-protein examples when fixed placements are
  part of the model;
- choose the spherical example when your geometry is easier to describe as
  spheres than as boxes.

## 5. Use Examples as Migration Guides

If you are migrating from Packmol, the examples are more useful than a raw API
index because they show the same physical setup in the language of this crate.
Follow them in this order:

1. read the corresponding `.inp` script or Python example;
2. identify the target count, restraint region, and iteration budget;
3. map those parts into your own system;
4. only then drill into the [API Reference](api/index.md) for exact type
   signatures.

## Takeaway

Examples in `molpack` are not decorations. They are executable case studies.
Choose the workload that matches your geometry, run it unchanged, and only then
generalize to your own system.
