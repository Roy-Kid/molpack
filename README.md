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

**Rust**

```toml
[dependencies]
molcrafts-molpack = "0.1"
```

**Python**

```bash
pip install molcrafts-molpack
```

## Quick start

**Rust**

```rust
use molpack::{InsideBoxRestraint, Molpack, Target};

let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
let radii = [1.52, 1.20, 1.20];

let target = Target::from_coords(&positions, &radii, 100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3]));

let result = Molpack::new()
    .tolerance(2.0)
    .pack(&[target], 200, Some(42))?;
```

**Python**

```python
import numpy as np
from molpack import InsideBox, Packer, Target

target = (
    Target.from_coords(positions, radii, count=100, elements=["O", "H", "H"])
    .with_name("water")
    .with_constraint(InsideBox([0, 0, 0], [40, 40, 40]))
)
result = Packer(tolerance=2.0).pack([target], max_loops=200, seed=42)
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

- **Rust** — `cargo doc --open` or [docs.rs](https://docs.rs/molcrafts-molpack)
- **Python** — [`python/docs/`](./python/docs/)
- **Concepts, architecture, extending** — rustdoc chapters (`molpack::concepts`, `molpack::architecture`, `molpack::extending`)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Bugs and feature requests via [GitHub Issues](https://github.com/MolCrafts/molpack/issues).

## License

BSD-3-Clause
