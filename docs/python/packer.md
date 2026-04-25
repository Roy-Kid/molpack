# Python Packer

`Molpack` is the run coordinator in the Python API. Targets describe molecules;
the packer describes how the run should behave.

## The Default Constructor Is Intentional

You start with an empty builder:

```python
from molpack import Molpack

packer = Molpack()
```

Then you turn the knobs you actually care about.

## The Knobs That Matter Most

Here is the usual starter chain:

```python
packer = (
    Molpack()
    .with_tolerance(2.0)
    .with_precision(0.01)
    .with_seed(42)
    .with_progress(True)
)
```

The rest of the builder surface exists for harder jobs:

- `with_inner_iterations(...)` for the inner optimizer budget
- `with_init_passes(...)` and `with_init_box_half_size(...)` for setup behavior
- `with_perturb(...)`, `with_perturb_fraction(...)`, and `with_random_perturb(...)` for stall recovery
- `with_parallel_eval(...)` for Rayon-backed pair evaluation
- `with_global_restraint(...)` for one rule shared by every target
- `with_handler(...)` for custom progress handling or early stop

Like `Target`, `Molpack` is immutable on the Python side. Each builder returns a
new object.

## Running a Job

The actual run call is compact:

```python
result = packer.pack(targets, max_loops=200)
```

The list of targets must be non-empty. `max_loops` is the outer-loop budget for
the run.

## Read the Result Like a Scientist, Not Like a Screenshot

The first four result fields to check are:

- `result.converged`
- `result.fdist`
- `result.frest`
- `result.natoms`

`result.positions` and `result.frame` matter, of course, but the diagnostics
tell you whether the geometry is actually acceptable.

```python
if not result.converged:
    print(f"still not there: fdist={result.fdist:.4f} frest={result.frest:.4f}")
```

## Handlers

Handlers are useful when you want progress updates or controlled early stops:

```python
class QuietStopper:
    def on_step(self, info):
        if info.loop_idx > 50 and info.fdist < 1e-3 and info.frest < 1e-3:
            return True
        return None

packer = Molpack().with_handler(QuietStopper())
```

Any subset of `on_start`, `on_step`, and `on_finish` is allowed.

## Errors You Will Actually See

Packing failures raise typed exceptions rooted at `PackError`, including:

- `NoTargetsError`
- `InvalidPBCBoxError`
- `ConflictingPeriodicBoxesError`
- `MaxIterationsError`
- `ConstraintsFailedError`

That means you can either catch one specific case or handle all packing failures
through the base class.

## Takeaway

Keep the split clean: targets describe the molecules, the packer describes the
run. If a setting changes the whole job, it probably belongs on `Molpack`.
