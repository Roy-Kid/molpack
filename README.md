# molcrafts-molpack

[![CI](https://github.com/MolCrafts/molpack/actions/workflows/ci.yml/badge.svg)](https://github.com/MolCrafts/molpack/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/molcrafts-molpack.svg)](https://crates.io/crates/molcrafts-molpack)
[![PyPI](https://img.shields.io/pypi/v/molcrafts-molpack.svg)](https://pypi.org/project/molcrafts-molpack/)
[![docs.rs](https://img.shields.io/docsrs/molcrafts-molpack)](https://docs.rs/molcrafts-molpack)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)

Packmol-grade molecular packing in pure Rust, with Python bindings.
Part of the [molrs](https://github.com/MolCrafts/molrs) toolkit.

## What it does

Given several molecule types with copy counts and geometric restraints, produce a non-overlapping packing by minimizing a penalty objective (GENCAN quasi-Newton with bound constraints, three-phase schedule: per-type init → geometric pre-fit → main loop with `movebad` heuristic). Output matches Packmol within tolerance; byte-identical for the 5 canonical examples checked into this repo.

## Core API

| Type | Role |
|---|---|
| `Molpack` | Packing configuration builder — `new()`, `tolerance`, `precision`, `pbc`, `add_handler`, `add_restraint`, `pack(&targets, max_loops, seed)` |
| `Target` | One molecule type — `Target::new(frame, count)` / `Target::from_coords(coords, radii, count)`, `with_name`, `with_restraint`, `with_restraint_for_atoms`, `with_relaxer`, `fixed_at`, `constrain_rotation_x/y/z` |
| `Restraint` trait | Soft-penalty contract: `f(x, scale, scale2) -> F`, `fg(x, scale, scale2, g) -> F`. `Send + Sync + Debug`. Accumulate gradient into `g` with `+=`. |
| 15 concrete `*Restraint` structs | Packmol kinds 2–15 as independent public types, each holding its own semantically-named fields |
| `Region` trait + `And`/`Or`/`Not` + `FromRegion<R>` | Geometric predicates with signed-distance; combinators compose; `FromRegion` lifts any Region to a quadratic-exterior-penalty Restraint |
| `Handler` trait | Observer callbacks — `on_start`, `on_initial`, `on_step`, `on_phase_start`, `on_phase_end`, `on_inner_iter`, `on_finish`, `should_stop` (all default no-op except `on_step`). Built-ins: `NullHandler`, `ProgressHandler`, `EarlyStopHandler`, `XYZHandler`. |
| `Relaxer` trait | Per-target reference-geometry modification between outer loops. Built-in: `TorsionMcRelaxer`. |

## Quick start

```rust
use molpack::{InsideBoxRestraint, Molpack, Target};

let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
let radii     = [1.52, 1.20, 1.20];

let target = Target::from_coords(&positions, &radii, 100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0, 40.0, 40.0]));

let result = Molpack::new()
    .tolerance(2.0)
    .precision(0.01)
    .pack(&[target], 200, Some(42))?;

println!("converged={}, natoms={}", result.converged, result.natoms());
# Ok::<(), molpack::PackError>(())
```

## Restraint scopes

Restraints attach at three granularities — all use the same `Restraint` trait:

```rust
use molpack::{InsideBoxRestraint, InsideSphereRestraint, BelowPlaneRestraint, Molpack, Target};

let target_a = Target::from_coords(&pos, &rad, 100)
    // per-target, all atoms:
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3]))
    // per-target, atom subset (1-based indexing per Packmol):
    .with_restraint_for_atoms(&[1, 2], BelowPlaneRestraint::new([0.0, 0.0, 1.0], 5.0));

let result = Molpack::new()
    // global — broadcast to every target ≡ each target.with_restraint(r.clone()):
    .add_restraint(InsideSphereRestraint::new([20.0; 3], 30.0))
    .pack(&[target_a], 200, Some(42))?;
# Ok::<(), molpack::PackError>(())
```

Scope equivalence: `Molpack::add_restraint(r) ≡ for t in targets { t.with_restraint(r.clone()) }` — no separate global-storage path inside the packer.

## Custom restraints (user plugins)

Any type that `impl Restraint` composes with built-ins identically — direction-3 design, no second-class citizenship:

```rust
use molrs::types::F;
use molpack::{InsideBoxRestraint, Restraint, Target};

#[derive(Debug, Clone, Copy)]
struct AboveZ { z_min: F }

impl Restraint for AboveZ {
    fn f(&self, pos: &[F; 3], _scale: F, scale2: F) -> F {
        let v = (self.z_min - pos[2]).max(0.0);
        scale2 * v * v
    }
    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let v = (self.z_min - pos[2]).max(0.0);
        if v > 0.0 { g[2] += -2.0 * scale2 * v; }
        self.f(pos, scale, scale2)
    }
}

let target = Target::from_coords(&pos, &rad, 10)
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3]))
    .with_restraint(AboveZ { z_min: 10.0 });
```

## Region composition

For geometric intersection / union that's more ergonomic than multiple `with_restraint` calls, use `Region` + `FromRegion`:

```rust
use molpack::{FromRegion, InsideSphereRegion, RegionExt, Target};

// Spherical shell: inside outer ∧ NOT inside inner
let shell = InsideSphereRegion::new([0.0; 3], 10.0)
    .and(InsideSphereRegion::new([0.0; 3], 5.0).not());

let target = Target::from_coords(&pos, &rad, 100)
    .with_restraint(FromRegion(shell));
```

`Region` methods: `contains`, `signed_distance`, `signed_distance_grad` (default finite-difference; concrete regions override analytically), `bounding_box`. Combinators `And` / `Or` / `Not` are zero-cost type algebra (analytic chain-rule gradients).

## Examples

5 canonical Packmol workloads reproduced under `examples/`:

```bash
cargo run -p molcrafts-molpack --release --example pack_mixture
cargo run -p molcrafts-molpack --release --example pack_bilayer
cargo run -p molcrafts-molpack --release --example pack_spherical
cargo run -p molcrafts-molpack --release --example pack_interface
cargo run -p molcrafts-molpack --release --example pack_solvprotein
```

Plus 2 demonstrations not wired into the regression suite (heavier, ad-hoc):

```bash
cargo run -p molcrafts-molpack --release --example mc_fold_chain        # pure torsion-MC fold (no packing)
cargo run -p molcrafts-molpack --release --example pack_polymer_vesicle # PE chains + vesicle
```

## Testing and benchmarks

```bash
# Correctness
cargo test -p molcrafts-molpack                                    # unit + integration
cargo test -p molcrafts-molpack --release --test examples_batch -- --ignored   # Packmol-equivalence regression (all 5 canonical cases)

# Performance
cargo bench -p molcrafts-molpack --bench pack_end_to_end -- mixture    # catastrophic-regression alarm (e2e, one workload)
cargo bench -p molcrafts-molpack --bench run_phase                     # hot-path microbench
cargo bench -p molcrafts-molpack --bench run_iteration
cargo bench -p molcrafts-molpack --bench evaluate_unscaled
cargo bench -p molcrafts-molpack --bench objective_dispatch
```

## Documentation

Primary documentation is delivered as rustdoc chapters; run `cargo doc
--open -p molcrafts-molpack` or visit [docs.rs](https://docs.rs/molcrafts-molpack).
The chapters live under `docs/` and are included into the rustdoc output
as top-level modules:

- `molpack::getting_started` — install, hello-world packing, the
  three restraint scopes, handlers, relaxers, running the canonical
  examples.
- `molpack::concepts` — every abstraction defined in one place
  (Restraint / Region / Relaxer / Handler / Objective / Target /
  Molpack / PackContext), the scope equivalence law, the two-scale
  contract, the direction-3 extension pattern.
- `molpack::architecture` — module map, dependency graph, full
  `pack()` lifecycle diagram, hot-path `evaluate()` walkthrough,
  invariants, design decisions.
- `molpack::extending` — tutorials for writing your own
  `Restraint` / `Region` / `Handler` / `Relaxer`, testing and
  benchmarking discipline, common pitfalls, contributing flow.

Python binding:

```bash
pip install molcrafts-molpack
```

```python
import numpy as np
from molpack import Target, Packer, InsideBox

target = (
    Target.from_coords(positions, radii, count=100, elements=["O", "H", "H"])
    .with_name("water")
    .with_constraint(InsideBox([0, 0, 0], [40, 40, 40]))
)
result = Packer(tolerance=2.0).pack([target], max_loops=200, seed=42)
```

The Python docs are a [Zensical](https://zensical.org/) site under
[`python/docs/`](./python/docs/); Python ports of the five canonical
Packmol workloads live under [`python/examples/`](./python/examples/).

Reference material outside rustdoc:

- [`docs/packmol_alignment.md`](./docs/packmol_alignment.md) —
  Packmol kind-number ↔ Rust struct mapping with Fortran pointers.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for dev-environment setup, test instructions,
and the PR checklist. Bug reports and feature requests go through
[GitHub Issues](https://github.com/MolCrafts/molpack/issues).

## License

BSD-3-Clause
