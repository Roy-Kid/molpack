# Rust

Use the Rust API when packing is part of a Rust program, when you need
structured convergence diagnostics, or when you are extending molpack itself.

```rust
use std::sync::Arc;
use molpack::{GenCanPack, PackEngine, RegionRestraint, Target};
use molrs::core::Cuboid;
use ndarray::array;

let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
let radii = [1.52, 1.20, 1.20];

let water = Target::from_coords(&positions, &radii, 100)
    .with_name("water")
    .with_restraint(RegionRestraint(Arc::new(Cuboid::new(array![0.0, 0.0, 0.0], array![40.0, 40.0, 40.0]))));

let result = GenCanPack::new().with_seed(42).run(&[water], 200)?;
let frame = result.frame;
```

The shared builders (`with_seed`, `with_tolerance`, handlers, boxes, …) and
the terminal `run` come from the `PackEngine` trait, so it has to be in scope.
`GenCanPack` is the rigid-body entry; `CbmcGrow` is the chain-growth one.

Each entry is a **single-stage preset**: calling `.run(...)` on `GenCanPack`
or `CbmcGrow` drives exactly one packing algorithm end to end (internally,
`Pipeline::single(self).run(...)`). When a pack needs more than one algorithm
in sequence — grow a chain, then push it apart with rigid-body descent —
compose stages directly with `Pipeline` instead of chaining separate runs:

```rust
use molpack::{CbmcGrow, GenCanPack, PackEngine, Pipeline};

let result = Pipeline::new()
    .with_stage(CbmcGrow::new(prior))
    .with_stage(GenCanPack::new())
    .run(&targets, max_loops)?;
```

`Pipeline` drives every stage through the same lifecycle a preset uses,
continuing from the first stage's placements rather than re-placing from
scratch. See [Composing stages](../extending.md#composing-stages) for the
full walkthrough, including the shared-settings rule and the two stage
combinators (`with_repeat`, `with_guarded`).

## Install

```bash
cargo add molcrafts-molpack
```

Feature flags:

| Feature | Enables |
|---|---|
| `io` | Template reading and output writing through the molrs reader and writer of each file's format (`script::StructureFormat`), and `XYZHandler`. |
| `cli` | The `molpack` binary plus `io`. |
| `rayon` | Parallel objective evaluation. |

molpack has no `ff` feature: a force-field optimizer bound through
`with_optimizer` (molrs's `Lbfgs`) needs molrs's `ff` feature on your own
`molcrafts-molrs` dependency.

## Pages

- [Quickstart](getting-started.md) walks through a first target and run.
- [Restraints and PBC](restraints-and-pbc.md) explains target-level,
  atom-subset, global, and periodic restraints.
- [Handlers and Optimizers](handlers-optimizers.md) covers progress output,
  observers, early stop, trajectory dumping, and in-loop conformation sampling.
- [Examples](examples.md) lists the checked-in Rust workloads.
