# Quickstart

A Rust packing job has three parts:

1. Build one `Target` per molecule species.
2. Attach at least one spatial restraint to each mobile target, or use an
   engine-level global restraint.
3. Pick an engine entry and call `run(&targets, max_loops)` on it —
   `GenCanPack` for rigid-body packing, `CbmcGrow` for chain growth.

## One molecule type in a box

```rust
use std::sync::Arc;
use molpack::{GenCanPack, PackEngine, RegionRestraint, Target};
use molrs::core::Cuboid;
use ndarray::array;

let water_positions = [
    [0.0, 0.0, 0.0],
    [0.96, 0.0, 0.0],
    [-0.24, 0.93, 0.0],
];
let water_radii = [1.52, 1.20, 1.20];

let water = Target::from_coords(&water_positions, &water_radii, 100)
    .with_name("water")
    .with_restraint(RegionRestraint(Arc::new(Cuboid::new(array![0.0, 0.0, 0.0], array![40.0, 40.0, 40.0]))));

let result = GenCanPack::new()
    .with_tolerance(2.0)
    .with_seed(42)
    .run(&[water], 200)?;

let natoms = result.natoms();
println!("packed {natoms} atoms");
```

`run()` returns a `State`. The packed `molrs::core::Frame` is its `frame`
field, alongside the convergence diagnostics:

```rust
let result = GenCanPack::new().with_seed(42).run(&targets, 200)?;
println!("converged={} fdist={} frest={}", result.converged, result.fdist, result.frest);
let frame = result.frame;
```

## Builder defaults

Every tuning knob except `max_loops` has a Packmol-compatible default. Set a
builder value only when you need to change the default:

```rust
let engine = GenCanPack::new()
    .with_tolerance(2.0)
    .with_precision(0.01)
    .with_inner_iterations(20)
    .with_seed(42);
```

`max_loops` is positional because the right iteration budget depends on system
size and packing difficulty.

## One engine, one run

`run` takes the entry **by value**, so an engine is consumed by the run it
performs. Build a fresh `GenCanPack` (or `CbmcGrow`) for each pack; a second
`run` on the same value does not compile. This is what makes it impossible to
lose an engine's handler set on a repeat call.

## Targets are snapshots

`Target` is a builder value. `run()` snapshots the target slice at call time;
mutating or rebuilding a target after that does not affect an already running
pack.

## Chaining two engines

Mixed rigid + grown packs are not a single call, and non-convergence never
falls back silently. Two explicit shapes:

```rust
use molpack::{CbmcGrow, GenCanPack, PackEngine, Target};

let grown = CbmcGrow::new(prior)
    .with_density(0.9)
    .with_seed(42)
    .run(&[chain.clone()], 60)?;

// Push-off: continue the SAME free targets on the grown state. The seeded
// run skips `initial()` and pushes remaining contacts apart by rigid-body
// descent; the cell travels with the seed.
let pushed = GenCanPack::new()
    .with_restart(&grown)
    .with_seed(42)
    .run(&[chain], 60)?;

// Fixed matrix: freeze the first result, pack new species around it.
let result = GenCanPack::new()
    .with_seed(42)
    .run(&[Target::fixed_from(&pushed.frame), solvent], 200)?;
```

`GenCanPack::with_restart(&result)` carries the placement solution over
verbatim (bitwise — no frame round-trip); `Target::fixed_from(&result.frame)`
wraps that frame as one fixed target with its coordinates kept
verbatim.
