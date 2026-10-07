# Packer

`GencanPack` drives the GENCAN-based three-phase optimizer. It is one of
two engines — the other, `CbmcGrow`, grows chains instead of
placing rigid bodies; see [Chain growth](growth.md). You pick the
algorithm by picking the engine, and both share the builder names below.
All tuning is through `with_*` builder methods — the `GencanPack`
constructor takes no arguments.

## Constructor

```python
from molpack import GencanPack

packer = GencanPack()
```

All defaults match Packmol's reference behaviour. Override any of
them via the builders below.

## Builder methods

Every builder returns a **new** `GencanPack`:

```python
packer = (
    GencanPack()
    .with_tolerance(2.0)            # minimum allowed pairwise distance (Å)
    .with_precision(0.01)           # convergence threshold on fdist and frest
    .with_inner_iterations(20)      # GENCAN inner-loop cap (Packmol `maxit`)
    .with_init_passes(0)            # initial compaction passes (Packmol `nloop0`; 0 = auto)
    .with_init_box_half_size(1000)  # hard bound on init placement (Packmol `sidemax`)
    .with_perturb(0.05, False, True)  # stall heuristic: fraction, random pick, on/off
    .with_seed(42)                  # deterministic RNG
    .with_parallel_eval(False)      # rayon-backed pair-kernel eval (opt-in)
    .with_log_level("progress")     # LAMMPS-style screen output
)
```

Runs are silent by default (`"quiet"`); `.with_log_level(level)` picks
`"quiet"`, `"summary"`, `"progress"` (LAMMPS-style thermo lines) or
`"verbose"`, and `.with_log_frequency(n)` prints every *n*-th loop.

A few more builders cover specific needs:

```python
packer = (
    packer
    .with_periodic_box([0, 0, 0], [30, 30, 30])  # fully-periodic cell (Packmol `pbc`)
    .with_density(0.9)                            # or: size a cubic cell from g/cm³
    .with_avoid_overlap(True)                     # reject init placements onto a fixed molecule
)
```

`with_tolerance`, `with_precision`, `with_seed`, `with_periodic_box`,
`with_density`, `with_parallel_eval`, `with_log_level`, `with_callback`,
and `with_global_restraint` are the shared engine builders — they exist
on `CbmcGrow` too. The rest are GENCAN-only.

## Global restraints

Attach a single restraint to every target in a pack:

```python
packer = packer.with_global_restraint(
    molrs.core.Cuboid([0, 0, 0], [40, 40, 40])
)
```

Semantically equivalent to calling `.with_restraint(r)` on every
target.

## Callbacks

Attach any object implementing some subset of `on_start(ntotat, ntotmol)`,
`on_step(step, sys) -> bool | None` (`sys` is a `PackSystemView`, valid
only inside the call), `on_finish()`:

```python
class MyCallback:
    def on_step(self, step, sys):
        print(f"phase={step.phase} loop={step.loop_idx} fdist={step.fdist:.3f}")
        return None  # or True to request early stop

packer = packer.with_callback(MyCallback())
```

Returning `True` from `on_step` halts the run at the next check. See
the `Callback` Protocol in `molpack`.

## Periodic boundaries

PBC is declared on the engine via `.with_periodic_box(min, max)`; a region
only confines. See
[Periodic boundaries](periodic-boundaries.md).

## Running

```python
result = packer.run(targets, max_loops=200)
```

- `targets`   — list of `Target` objects (must be non-empty).
- `max_loops` — per-phase outer-iteration budget.

Raises one of the typed `PackError` subclasses on failure
(`NoTargetsError`, `InvalidPbcBoxError`, …).

`run()` is the engine's only terminal verb and it consumes the engine —
one engine, one run. Calling `run()` twice on the same object raises
`RuntimeError`; build a fresh `GencanPack` for the next pack.

## State

```python
result.positions   # (N, 3) float64 ndarray — packed coordinates
result.frame       # molrs.core.Frame — topology-complete packed frame
result.elements    # list[str]  — one entry per atom
result.natoms      # int
result.converged   # bool — True iff both fdist and frest < precision
result.fdist       # float — final distance-violation sum
result.frest       # float — final restraint-violation sum
result.degraded    # int — growth-only; always 0 on the GencanPack path
result.intra       # IntraResidual — same-copy scored / exempted minima (Å)
```

Inspect convergence:

```python
if not result.converged:
    print(f"not converged: fdist={result.fdist:.4f} frest={result.frest:.4f}")
```

`State.frame` is the packed frame. Pass it to a writer of your
choice (e.g. `molrs.io.write_pdb`). molpack does **not** provide
writers.

## Composing stages

`GencanPack`, `CbmcGrow`, and `LatticeGrow` are each a **single-stage
preset**: calling `.run(...)` on one of them drives exactly one packing
algorithm end to end. When a pack needs more than one algorithm in the same
run — grow a chain, then compact it with rigid-body descent — compose the
stage objects with `Pipeline` instead of writing two separate `run()` calls:

```python
from molpack import CbmcGrow, GencanPack, Pipeline, TorsionPrior

prior = TorsionPrior.uniform()  # see the Chain growth guide for a real prior
result = (
    Pipeline([CbmcGrow(prior), GencanPack()])
    .with_seed(7)
    .run(targets, max_loops=200)
)
```

`Pipeline([stage, ...])` accepts a list of stage objects at construction, and
`.with_stage(x)` appends one more; both accept `GencanPack`, `CbmcGrow`, and
`LatticeGrow` instances. The pipeline runs every stage in one lifecycle,
continuing from the previous stage's placements rather than starting over —
the same continuation `.with_restart` gives you across two separate runs (see
[Staging a mixed pack](growth.md#staging-a-mixed-pack)), but in one call and
one `State`.

### Shared settings live on the pipeline, not the stages

`with_tolerance`, `with_precision`, `with_seed`, `with_periodic_box`,
`with_density`, and the rest of the shared builders are one ruler for the
whole run — declare them on the `Pipeline`, never on a stage object that goes
into one. A stage object that still carries a non-default shared setting when
it enters a `Pipeline` raises `ValueError`, naming both the offending stage
and the setting, rather than silently picking one of two conflicting values:

```python
Pipeline([CbmcGrow(prior), GencanPack().with_seed(7)]).run(targets, max_loops=200)
# ValueError: preset `gencan` carries a non-default `seed` inside a pipeline;
#             set `seed` on the Pipeline instead (shared settings are one ruler)
```

### Callbacks travel with their stage

A stage object's own `.with_callback(...)` callback is **not** dropped when
that stage is composed into a `Pipeline` — it is adopted into the pipeline's
callback set and keeps firing for every stage in the run, not only the one it
was attached to:

```python
class CountSteps:
    def __init__(self):
        self.count = 0
    def on_step(self, step, sys):
        self.count += 1

counter = CountSteps()
result = Pipeline([GencanPack().with_callback(counter)]).run(targets, max_loops=200)
assert counter.count > 0
```

### `StepReport.stage` names the running algorithm

Inside a multi-stage pipeline, `step.stage` on every `StepReport` a callback
receives tells you which stage emitted that step:

```python
class WatchStages:
    def on_step(self, step, sys):
        s = step.stage
        print(f"stage {s.index + 1}/{s.total} ({s.name}) loop={step.loop_idx}")
        return None

Pipeline([CbmcGrow(prior), GencanPack()]).with_callback(WatchStages()).run(targets, max_loops=200)
```

`step.stage.index` is 0-based and increases monotonically over the run,
`step.stage.total` is the number of stages the pipeline holds, and
`step.stage.name` is the stage's own name — `"gencan"`, `"growth"`, or
`"lattice"`. A single-engine run (`GencanPack().run(...)` directly, with no
`Pipeline`) reports the same triple, with `index = 0` and `total = 1`.

### Composition errors

- An empty pipeline (no stages by the time `.run(...)` is called) raises
  `ValueError`.
- A stage-ordering error — a stage whose `requires` the previous stage did
  not `guarantee` — raises `ValueError` naming the stage.
- A preset entering the pipeline with a non-default shared setting (above)
  raises `ValueError` naming the stage and the setting.
- An object that is not one of `GencanPack`, `CbmcGrow`, or `LatticeGrow`
  passed to `Pipeline([...])` or `.with_stage(x)` raises `TypeError`, listing
  the three supported engines.

## Reproducibility

Packing is deterministic for a given `(targets, tolerance, precision,
seed)` tuple. Capture the builder chain and `max_loops` to reproduce
a result later.
