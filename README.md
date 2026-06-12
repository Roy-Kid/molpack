# molcrafts-molpack

[![CI](https://github.com/MolCrafts/molpack/actions/workflows/ci.yml/badge.svg)](https://github.com/MolCrafts/molpack/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/molcrafts-molpack.svg)](https://crates.io/crates/molcrafts-molpack)
[![PyPI](https://img.shields.io/pypi/v/molcrafts-molpack.svg)](https://pypi.org/project/molcrafts-molpack/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![ty](https://img.shields.io/badge/type--checked-ty-informational)](https://github.com/astral-sh/ty)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)

Packmol-grade molecular packing in pure Rust, with Python bindings.
Part of the [molrs](https://github.com/MolCrafts/molrs) toolkit.

## Install

**CLI**

```bash
cargo install molcrafts-molpack --features cli
```

**Rust library**

```bash
cargo add molcrafts-molpack
```

**Python**

```bash
pip install molcrafts-molpack molcrafts-molrs
```

## CLI

The `molpack` binary accepts the same `.inp` script format as Packmol, making it a drop-in replacement. Relative file paths in the script are resolved against the script's own directory (file-arg mode) or the current directory (stdin mode).

```bash
# File argument (paths resolved relative to the .inp file's directory)
molpack mixture.inp

# Stdin — compatible with Packmol usage
molpack < mixture.inp
cat mixture.inp | molpack
```

**Supported `.inp` keywords**

| Keyword | Description |
|---|---|
| `tolerance <f>` | Minimum atom–atom distance (Å). Default 2.0 |
| `seed <n>` | Random seed for reproducibility |
| `filetype <fmt>` | Input format for all structure files |
| `output <path>` | Output file path (format inferred from extension) |
| `nloop <n>` | Maximum outer-loop iterations. Default 400 |
| `pbc x y z` / `pbc x0 y0 z0 x1 y1 z1` | Periodic box: side lengths from origin, or explicit min+max |
| `avoid_overlap <no\|false\|0>` | Reject initial placements that land inside a fixed molecule (default on; faithful to Packmol's `initial.f90`) |
| `structure <file>` … `end structure` | One block per molecule type |
| `number <n>` | Copies to pack |
| `inside\|outside box x0 y0 z0 x1 y1 z1` | Axis-aligned box restraint |
| `inside\|outside cube x0 y0 z0 d` | Axis-aligned cube (origin corner + side) |
| `inside\|outside sphere cx cy cz r` | Sphere restraint |
| `inside\|outside ellipsoid a1 a2 a3 b1 b2 b3 c` | Ellipsoid (center, semi-axes, exponent) |
| `inside\|outside cylinder a1 a2 a3 d1 d2 d3 r l` | Finite cylinder (center, axis, radius, length) |
| `over\|above plane nx ny nz d` | Half-space (above) restraint |
| `below plane nx ny nz d` | Half-space (below) restraint |
| `center` | Center molecule at origin before packing |
| `fixed x y z ex ey ez` | Fix molecule at position + Euler angles |
| `atoms i j …` … `end atoms` | Per-atom-subset restraints |

All 12 box/cube/sphere/ellipsoid/cylinder/plane restraints are reachable from
both whole-molecule and `atoms … end atoms` blocks. The Gaussian-surface
restraints (`AboveGaussianRestraint` / `BelowGaussianRestraint`) are available
through the Rust/Python API but not yet exposed as `.inp` keywords.

**Extended format support** — beyond Packmol's PDB/XYZ, molpack also reads SDF/MOL, LAMMPS dump, and LAMMPS data files. Set `filetype` to `sdf`, `lammps_dump`, or `lammps_data`, or use the matching file extension:

| Format | Read | Write | Extension / filetype |
|---|---|---|---|
| PDB | ✓ | ✓ | `.pdb` / `pdb` |
| XYZ | ✓ | ✓ | `.xyz` / `xyz` |
| SDF / MOL | ✓ | — | `.sdf`, `.mol` / `sdf` |
| LAMMPS dump | ✓ | ✓ | `.lammpstrj` / `lammps_dump` |
| LAMMPS data | ✓ | — | `.data` / `lammps_data` |

## Quick start

**Rust**

```rust
use molpack::{InsideBoxRestraint, Molpack, Target};

let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
let radii = [1.52, 1.20, 1.20];

let target = Target::from_coords(&positions, &radii, 100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3]));

// `pack` returns the packed, topology-complete `molrs::Frame`.
// Every tuning knob has a Packmol-matching default, so `new().pack(...)`
// is a complete call; `200` is the outer-loop budget.
let frame = Molpack::new().pack(&[target], 200)?;

// For full diagnostics, use `pack_with_report` → `PackResult`
// (`frame`, `fdist`, `frest`, `converged`).
let report = Molpack::new().pack_with_report(&[target], 200)?;
```

**Python**

```python
import molrs
from molpack import InsideBoxRestraint, Molpack, Target

frame = molrs.read_pdb("water.pdb")

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
)
frame = Molpack().pack([water], max_loops=200)
```

## Examples

Five canonical workloads ship in `examples/` (they need the `io` feature
to read the bundled structure files):

```bash
cargo run --release --example pack_mixture     --features io   # 1000 water + 400 urea in a cube
cargo run --release --example pack_bilayer     --features io   # membrane leaflets via per-atom plane restraints
cargo run --release --example pack_interface   --features io   # water + chloroform around a fixed molecule
cargo run --release --example pack_spherical   --features io   # concentric lipid/water shells (largest case)
cargo run --release --example pack_solvprotein --features io   # fixed protein solvated in a sphere (avoid_overlap)
```

The same workloads run through the CLI from their bundled `.inp` scripts,
e.g. `cargo run --release --features cli --bin molpack -- examples/pack_mixture/mixture.inp`.
Python equivalents are in [`python/examples/`](./python/examples/).

## Testing

```bash
cargo test                                                  # unit + integration
cargo test --release --test examples_batch -- --ignored     # Packmol regression (all 5 workloads)
cd python && maturin develop --release && pytest            # Python wheel
```

## Documentation

- **Guide** — the Markdown chapters under [`docs/`](./docs/): install,
  getting started, concepts, examples, Packmol parity, architecture, and
  extending.
- **Rust API** — `cargo doc --open`, or [docs.rs](https://docs.rs/molcrafts-molpack).
  The four long-form chapters (getting started, concepts, architecture,
  extending) are also embedded in the rustdoc as
  `molpack::getting_started`, `molpack::concepts`,
  `molpack::architecture`, and `molpack::extending`.
- **Python** — [`python/docs/`](./python/docs/).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Bugs and feature requests via [GitHub Issues](https://github.com/MolCrafts/molpack/issues).

## License

BSD-3-Clause

## References

- Martínez, L.; Andrade, R.; Birgin, E. G.; Martínez, J. M.
  **PACKMOL: A package for building initial configurations for molecular dynamics simulations.**
  *J. Comput. Chem.* **2009**, *30* (13), 2157–2164.
  https://doi.org/10.1002/jcc.21224
