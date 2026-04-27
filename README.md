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
| `avoid_overlap <yes\|no>` | Parsed for compatibility |
| `structure <file>` … `end structure` | One block per molecule type |
| `number <n>` | Copies to pack |
| `inside box x0 y0 z0 x1 y1 z1` | Axis-aligned box restraint |
| `inside sphere cx cy cz r` | Sphere restraint |
| `outside sphere cx cy cz r` | Exclusion sphere |
| `over plane nx ny nz d` | Half-space (above) restraint |
| `below plane nx ny nz d` | Half-space (below) restraint |
| `center` | Center molecule at origin before packing |
| `fixed x y z ex ey ez` | Fix molecule at position + Euler angles |
| `atoms i j …` … `end atoms` | Per-atom-subset restraints |

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

let result = Molpack::new()
    .with_tolerance(2.0)
    .with_seed(42)
    .pack(&[target], 200)?;
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
result = (
    Molpack()
    .with_tolerance(2.0)
    .with_seed(42)
    .pack([water], max_loops=200)
)
```

## Examples

```bash
cargo run --release --example pack_mixture
cargo run --release --example pack_bilayer
cargo run --release --example pack_spherical
cargo run --release --example pack_interface
cargo run --release --example pack_solvprotein
```

Python equivalents are in [`python/examples/`](./python/examples/).

## Testing

```bash
cargo test                                                               # unit + integration
cargo test --release --test examples_batch -- --ignored                 # Packmol regression
cargo bench --bench pack_end_to_end -- mixture                          # e2e benchmark
```

## Documentation

- **Site** — [molcrafts.github.io/molpack](https://molcrafts.github.io/molpack/)
- **Rust API** — `cargo doc --open` or [docs.rs](https://docs.rs/molcrafts-molpack)
- **Python** — [`python/docs/`](./python/docs/)

The site has four top-level tabs: **Home**, **Get Started** (install,
getting started), **Guide** (concepts, examples, Packmol parity), and
**Internals** (architecture, extending). The same long-form chapters
are embedded in the rustdoc as `molpack::getting_started`,
`molpack::concepts`, `molpack::architecture`, and `molpack::extending`.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Bugs and feature requests via [GitHub Issues](https://github.com/MolCrafts/molpack/issues).

## License

BSD-3-Clause

## References

- Martínez, L.; Andrade, R.; Birgin, E. G.; Martínez, J. M.
  **PACKMOL: A package for building initial configurations for molecular dynamics simulations.**
  *J. Comput. Chem.* **2009**, *30* (13), 2157–2164.
  https://doi.org/10.1002/jcc.21224
