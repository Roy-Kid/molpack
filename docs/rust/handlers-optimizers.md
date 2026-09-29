# Handlers and Optimizers

Handlers observe a packing run. In-loop optimizers modify a copy's reference
geometry between GENCAN iterations.

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
    // Stricter than the default early stop: 20 % per 5 loops.
    .with_early_stop(EarlyStopHandler::new(20.0).with_patience(5));
```

Packmol runs every phase until it converges or reaches `nloop`. `GenCanPack`
additionally ends a phase once, at the user's radii (`radscale == 1`),
Packmol's `bestf` has improved by less than 10 % (its movebad threshold) over
10 loops — measured with Packmol's `fimprov` formula. A pack that stalls there
is almost always too dense, and a `converged == false` in minutes beats
running to `max_loops`. `with_early_stop(None)` restores Packmol's behaviour.

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

## In-loop optimizers

An in-loop optimizer reshapes a copy's reference geometry while it is being
placed. It is useful for flexible molecules that need to sample torsions to fit
their restraints. Bind one to the engine, naming the targets it applies to:

```rust
use molpack::{GenCanPack, OptimizeSelect, PackEngine, TorsionMcOptimizer};

let engine = GenCanPack::new().with_optimizer(
    OptimizeSelect::per_copy(["chain"]).with_environment(8.0),
    TorsionMcOptimizer::new(&graph)
        .with_temperature(0.5)
        .with_steps(20),
);
```

Each copy is relaxed on its own, so copies of one target start identical and
then diverge. After each call molpack re-evaluates the packing objective and
reverts the conformer if it got worse. `OptimizeSelect::joint` relaxes all
selected copies as one group. Any molrs `Optimizer` fits the same slot; binding
a force-field one such as `LBFGS` needs molrs's `ff` module (molpack's `ff`
feature forwards it).
