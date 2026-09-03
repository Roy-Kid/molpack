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
pip install molcrafts-molpack
```

## CLI

The `molpack` binary accepts Packmol-style `.inp` scripts. Use the official
[Packmol user guide](https://m3g.github.io/packmol/userguide.shtml) as the
script-language reference; molpack documents only invocation behavior and
extensions here. Relative file paths in the script are resolved against the
script's own directory (file-arg mode) or the current directory (stdin mode).

```bash
# File argument (paths resolved relative to the .inp file's directory)
molpack mixture.inp

# Stdin — compatible with Packmol usage
molpack < mixture.inp
cat mixture.inp | molpack
```

**molpack script additions**

| Addition | Description |
|---|---|
| `avoid_overlap <no\|false\|0>` | Disable the default fixed-solute initial-placement guard. Leave it on unless you need to reproduce a less guarded initialization. |
| `filetype sdf` | Read SDF/MOL inputs. Read-only. |
| `filetype lammps_dump` | Read LAMMPS dump inputs and write `.lammpstrj` outputs. |
| `filetype lammps_data` | Read LAMMPS data inputs. Read-only. |

Unknown top-level keywords are rejected instead of ignored. Output format is
inferred from the `output` extension.

## Quick start

**Rust**

```rust
use molpack::{GenCanPack, InsideBoxRestraint, PackEngine, Target};

let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
let radii = [1.52, 1.20, 1.20];

let target = Target::from_coords(&positions, &radii, 100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3]));

// You pick the algorithm by picking the entry: `GenCanPack` for rigid-body
// packing, `CbmcGrow` for chain growth. Both share the `PackEngine` builders
// and the single terminal verb `run`, which consumes the entry — one engine,
// one run. Every tuning knob has a Packmol-matching default, so
// `GenCanPack::new().run(...)` is a complete call; `200` is the outer-loop
// budget.
//
// `run` returns a `PackResult`: the packed, topology-complete `molrs::Frame`
// in `.frame`, plus `fdist`, `frest`, `converged`, `softened`.
let result = GenCanPack::new().run(&[target], 200)?;
let frame = result.frame;
```

**Python**

```python
import molrs
from molpack import GenCanPack, InsideBoxRestraint, Target

frame = molrs.io.read_pdb("water.pdb")

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
)
result = GenCanPack().run([water], max_loops=200)
frame = result.frame
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

Two measurement harnesses also live under `examples/`: `pack_peo`
(PEO-melt evaluation of the chain-growth solver against the rigid-body
path, needs `--features io`) and `mt_scaling` (parallel
speed-up-vs-size sweep, needs `--features rayon`).

## Testing

```bash
cargo test                                                  # unit + integration
cargo test --release --test examples_batch -- --ignored     # Packmol regression (all 5 workloads)
cargo bench --benches                                       # criterion regression benches (no io)
cd python && maturin develop --release && pytest            # Python wheel
```

## Documentation

- **Guide** — [`docs/`](./docs/), built with `zensical build`. Chapters:
  install, getting started, concepts, examples, Packmol parity, architecture,
  and extending.
- **Rust API** — `cargo doc --open`, or [docs.rs](https://docs.rs/molcrafts-molpack).
  The four long-form chapters (getting started, concepts, architecture,
  extending) are also embedded in the rustdoc as
  `molpack::getting_started`, `molpack::concepts`,
  `molpack::architecture`, and `molpack::extending`.
- **Python** — binding docs live under `docs/python/`, published as the
  **Python** section of the site.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Bugs and feature requests via [GitHub Issues](https://github.com/MolCrafts/molpack/issues).

## License

BSD-3-Clause

## References

- Martínez, L.; Andrade, R.; Birgin, E. G.; Martínez, J. M.
  **PACKMOL: A package for building initial configurations for molecular dynamics simulations.**
  *J. Comput. Chem.* **2009**, *30* (13), 2157–2164.
  https://doi.org/10.1002/jcc.21224
