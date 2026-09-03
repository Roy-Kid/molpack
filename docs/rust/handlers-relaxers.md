# Handlers and Relaxers

Handlers observe a packing run. Relaxers modify a target's reference geometry
between optimizer iterations.

## Screen output

Enable LAMMPS-style progress output through the builder:

```rust
use molpack::{GenCanPack, LogLevel, PackEngine};

let engine = GenCanPack::new()
    .with_log_level(LogLevel::Progress)
    .with_log_frequency(10);
```

The CLI enables screen output by default; library callers stay quiet unless you
opt in.

## Built-in handlers

```rust
use molpack::{EarlyStopHandler, GenCanPack, PackEngine, XYZHandler};

let engine = GenCanPack::new()
    .with_handler(Box::new(XYZHandler::new("traj.xyz", 10)))
    .with_handler(Box::new(EarlyStopHandler::new(1e-4)));
```

`with_handler` is a `PackEngine` builder, so it works the same on `CbmcGrow`.
Use handlers for progress logs, trajectory snapshots, custom observation, and
early stop. Handler callbacks receive an immutable `PackContext` view; they do
not mutate engine state.

## Custom handlers

Implement the `Handler` trait when you need structured events from a run:

```rust
use molpack::{Handler, PackContext, StepInfo};

#[derive(Debug)]
struct WatchFdist;

impl Handler for WatchFdist {
    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        eprintln!("phase={} loop={} fdist={}", info.phase, info.loop_idx, info.fdist);
    }
}
```

See [Extending](../extending.md) for a full custom-handler walkthrough.

## Relaxers

Relaxers update a molecule's reference geometry during packing. They are useful
for flexible molecules that need to sample torsions while being placed.

```rust
use molpack::{InsideSphereRestraint, Target, TorsionMcRelaxer};

let target = Target::from_coords(positions, radii, 1)
    .with_restraint(InsideSphereRestraint::new([0.0; 3], 20.0))
    .with_relaxer(
        TorsionMcRelaxer::new(&graph)
            .with_temperature(0.5)
            .with_steps(20),
    );
```

Relaxers require `count == 1` because every copy of a target shares one
reference geometry.
