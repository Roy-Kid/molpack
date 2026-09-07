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

Every `StepInfo` also says which packing algorithm emitted it. A **stage** is
one algorithm behind molpack's packing seam (the `Stage` trait), and
`info.stage` is a `StageInfo` carrying `index` (0-based position of the stage
in the run), `total` (how many stages the run has), and `name` (the stage's own
name, `"gencan"` for the rigid-body path). One engine entry drives one stage,
so a plain `GenCanPack` or `CbmcGrow` run reports `index = 0`, `total = 1`.

Two further callbacks bracket a whole stage, the way `on_phase_start` /
`on_phase_end` bracket one GENCAN phase:

```rust
use molpack::handler::StageInfo;
use molpack::{Handler, PackContext, StageOutcome, StepInfo};

struct WatchStages;

impl Handler for WatchStages {
    fn on_step(&mut self, _info: &StepInfo, _sys: &PackContext) {}

    fn on_stage_start(&mut self, info: &StageInfo) {
        eprintln!("stage {}/{} ({}) starting", info.index + 1, info.total, info.name);
    }

    fn on_stage_end(&mut self, info: &StageInfo, outcome: &StageOutcome, sys: &PackContext) {
        eprintln!(
            "stage {} converged={} degraded={} fdist={} frest={}",
            info.name, outcome.converged, outcome.degraded, sys.fdist, sys.frest,
        );
    }
}
```

Both have default no-op bodies, and a run driven by a single engine entry calls
neither; they are the seam a caller that chains stages itself brackets each
stage with. Note where the numbers come from:
`StageOutcome` reports only what the stage alone knows (`converged`,
`degraded`), while the violation maxima `fdist` and `frest` are read off the
post-stage `PackContext`, so every algorithm is judged by the same objective.

See [Extending](../extending.md) for a full custom-handler walkthrough and for
writing a stage of your own.

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
