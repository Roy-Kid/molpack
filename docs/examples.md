# Examples

Five canonical Packmol-equivalent workloads ship in `examples/`. Each
exercises a different combination of restraints, fixed placements, and
target counts.

| Workload | Rust example | Python example | Demonstrates |
|---|---|---|---|
| Mixture | `examples/pack_mixture` | `python/examples/pack_mixture.py` | two species in one box |
| Bilayer | `examples/pack_bilayer` | `python/examples/pack_bilayer.py` | atom-subset restraints to orient lipids |
| Interface | `examples/pack_interface` | `python/examples/pack_interface.py` | fixed structures plus mobile solvents |
| Spherical | `examples/pack_spherical` | `python/examples/pack_spherical.py` | nested regions, shell packing |
| Solvated protein | `examples/pack_solvprotein` | `python/examples/pack_solvprotein.py` | one fixed macromolecule, heterogeneous solvent |

## Run the Rust examples

```bash
cargo run --release --example pack_mixture
cargo run --release --example pack_bilayer
cargo run --release --example pack_interface
cargo run --release --example pack_spherical
cargo run --release --example pack_solvprotein
```

Each directory also contains the matching `.inp` script, so the CLI
form (`molpack examples/pack_mixture/mixture.inp`) and the programmatic
form can be compared side by side.

Two heavier ad-hoc demonstrations live alongside the canonical set
(not part of the regression suite):

```bash
cargo run --release --example mc_fold_chain
cargo run --release --example pack_polymer_vesicle
```

Optional progress / trajectory dumps:

```bash
MOLRS_PACK_EXAMPLE_PROGRESS=1 cargo run --release --example pack_mixture
MOLRS_PACK_EXAMPLE_XYZ=1      cargo run --release --example pack_mixture
```

## Run the Python examples

```bash
cd python
pip install -e .
python examples/pack_mixture.py
python examples/pack_bilayer.py
```

Examples that load real structure files also need `molcrafts-molrs`.
