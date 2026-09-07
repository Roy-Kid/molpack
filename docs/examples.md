# Examples

Five canonical Packmol-equivalent workloads ship in `examples/`. Each
exercises a different combination of restraints, fixed placements, and
target counts, and all five are covered by the regression suite
(`tests/examples_batch.rs`). A sixth program, `examples/pack_peo`, is a
measurement harness rather than a Packmol workload: it evaluates the
chain-growth solver against the rigid-body path on a PEO melt
(`cargo run --release --example pack_peo --features io -- grow 200 25 1.0 42`).
Polymer PEO scenes are Python-only (`python/examples/pack_peo_*.py`):
monomers from molrs SMILES + conformer, architecture from molpy
`PolymerBuilder`, then molpack packing. Linear melt and mixed
linear+star share `LatticeGrow` @ 2.0 Å then `GenCanPack.with_restart`;
the ring path is `pack_peo_topo.py ring` (named grower reject, rigid
`GenCanPack`); STL-confined linear PEO is `pack_peo_stl.py` (`StlRegion`
masks lattice sites outside the mesh, then `GenCanPack.with_restart`).

| Workload | Rust example | Python example | Molecules | Restraints | Demonstrates |
|---|---|---|---|---|---|
| Mixture | `examples/pack_mixture` | `python/examples/pack_mixture.py` | 1000 water + 400 urea | one `inside box` per type | the simplest two-component fill in a cube |
| Bilayer | `examples/pack_bilayer` | `python/examples/pack_bilayer.py` | water slabs + palmitoil lipids | per-atom `above`/`below plane` | orienting the two leaflets of a membrane |
| Interface | `examples/pack_interface` | `python/examples/pack_interface.py` | water + chloroform + 1 fixed t3 | `inside box` per solvent + one `fixed` molecule | a liquid–liquid interface around a fixed structure |
| Spherical | `examples/pack_spherical` | `python/examples/pack_spherical.py` | concentric lipid + water shells | per-atom radial `inside`/`outside sphere` | nested shell packing (largest / slowest case) |
| Solvated protein | `examples/pack_solvprotein` | `python/examples/pack_solvprotein.py` | 1 fixed protein + water + Na⁺ + Cl⁻ | `inside sphere` solvent around a `fixed` solute | dense solvation that relies on `avoid_overlap` |

## Run the Rust examples

The examples need the `io` feature so they can read the bundled structure
files:

```bash
cargo run --release --example pack_mixture     --features io
cargo run --release --example pack_bilayer     --features io
cargo run --release --example pack_interface   --features io
cargo run --release --example pack_spherical   --features io
cargo run --release --example pack_solvprotein --features io
```

Each directory also contains the matching `.inp` script
(`mixture.inp`, `bilayer-comment.inp`, `interface.inp`,
`spherical-comment.inp`, `solvprotein.inp`), so the CLI form and the
programmatic form can be compared side by side:

```bash
cargo run --release --features cli --bin molpack -- examples/pack_mixture/mixture.inp
```

Optional progress / trajectory dumps are gated behind environment
variables:

```bash
MOLRS_PACK_EXAMPLE_PROGRESS=1 cargo run --release --example pack_mixture --features io
MOLRS_PACK_EXAMPLE_XYZ=1      cargo run --release --example pack_mixture --features io
```

## Run the Python examples

The Python examples mirror the Rust ones and load their structure files
through the `molrs` dependency installed with molpack:

```bash
cd python
maturin develop --release
python examples/pack_mixture.py
python examples/pack_bilayer.py
python examples/pack_interface.py
python examples/pack_spherical.py
python examples/pack_solvprotein.py
python examples/pack_peo_linear.py 8 8 0.5 42
python examples/pack_peo_mix.py 4 2 4 4 0.5 42
python examples/pack_peo_topo.py star 4 8 0.5 42
python examples/pack_peo_stl.py 4 4 30 42
```

`python/examples/pack_water_cube.py` is a minimal standalone starter
(single-species cube) for a first look at the Python API.
