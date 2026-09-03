# Rust

Use the Rust API when packing is part of a Rust program, when you need
structured convergence diagnostics, or when you are extending molpack itself.

```rust
use molpack::{GenCanPack, InsideBoxRestraint, PackEngine, Target};

let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
let radii = [1.52, 1.20, 1.20];

let water = Target::from_coords(&positions, &radii, 100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [40.0, 40.0, 40.0],
        [false, false, false],
    ));

let result = GenCanPack::new().with_seed(42).run(&[water], 200)?;
let frame = result.frame;
```

The shared builders (`with_seed`, `with_tolerance`, handlers, boxes, …) and
the terminal `run` come from the `PackEngine` trait, so it has to be in scope.
`GenCanPack` is the rigid-body entry; `CbmcGrow` is the chain-growth one.

## Install

```bash
cargo add molcrafts-molpack
```

Feature flags:

| Feature | Enables |
|---|---|
| `io` | PDB, XYZ, SDF, and LAMMPS readers through `molrs_io`. |
| `cli` | The `molpack` binary plus `io`. |
| `rayon` | Parallel objective evaluation. |
| `ff` | Force-field-backed relaxers. |

## Pages

- [Quickstart](getting-started.md) walks through a first target and run.
- [Restraints and PBC](restraints-and-pbc.md) explains target-level,
  atom-subset, global, and periodic restraints.
- [Handlers and Relaxers](handlers-relaxers.md) covers progress output,
  observers, early stop, trajectory dumping, and in-loop conformation sampling.
- [Examples](examples.md) lists the checked-in Rust workloads.
