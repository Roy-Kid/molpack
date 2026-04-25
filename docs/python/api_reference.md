# Python API Reference

This page is the quick map of the Python binding. It is not meant to replace the
Rust API docs; it is here to make the Python surface easy to scan from one
place.

## Core Types

```python
from molpack import (
    Target,
    Molpack,
    PackResult,
    StepInfo,
)
```

- `Target` describes one molecule type plus its copy count and target-local geometry.
- `Molpack` owns run-wide settings and executes the job.
- `PackResult` exposes packed coordinates and convergence diagnostics.
- `StepInfo` is the snapshot passed to handler callbacks.

## Utility Types

```python
from molpack import Angle, Axis, CenteringMode
```

- `Angle` makes rotational units explicit.
- `Axis` identifies `X`, `Y`, or `Z`.
- `CenteringMode` controls how reference coordinates are centered before packing.

## Built-In Restraints

```python
from molpack import (
    InsideBoxRestraint,
    InsideSphereRestraint,
    OutsideSphereRestraint,
    AbovePlaneRestraint,
    BelowPlaneRestraint,
)
```

These cover the common Packmol-style geometries. `InsideBoxRestraint` also
carries periodic-boundary declarations.

## `Target` at a Glance

Construction:

```python
Target(frame, count: int)
```

Common builders:

- `.with_name(name)`
- `.with_restraint(r)`
- `.with_atom_restraint(indices, r)`
- `.with_centering(mode)`
- `.fixed_at(position)`
- `.with_orientation((ax, ay, az))`
- `.with_rotation_bound(axis, center, half_width)`
- `.with_perturb_budget(n)`

Common properties:

- `.name`
- `.natoms`
- `.count`
- `.elements`
- `.radii`
- `.is_fixed`

## `Molpack` at a Glance

Construction:

```python
Molpack()
```

Common builders:

- `.with_tolerance(t)`
- `.with_precision(p)`
- `.with_seed(seed)`
- `.with_progress(enabled)`
- `.with_parallel_eval(enabled)`
- `.with_global_restraint(r)`
- `.with_handler(handler)`
- `.with_inner_iterations(n)`
- `.with_init_passes(n)`
- `.with_init_box_half_size(h)`
- `.with_perturb(enabled)`
- `.with_perturb_fraction(f)`
- `.with_random_perturb(enabled)`

Execution:

```python
.pack(targets, max_loops=200) -> PackResult
```

## `PackResult` at a Glance

The fields you will use most often are:

- `.positions`
- `.frame`
- `.elements`
- `.natoms`
- `.converged`
- `.fdist`
- `.frest`

## Handler Protocol

A handler may implement any subset of:

```python
class Handler:
    def on_start(self, ntotat, ntotmol): ...
    def on_step(self, info): ...
    def on_finish(self): ...
```

Returning `True` from `on_step` requests an early stop.

## Exceptions

The main exception hierarchy is:

- `PackError`
- `ConstraintsFailedError`
- `MaxIterationsError`
- `NoTargetsError`
- `EmptyMoleculeError`
- `InvalidPBCBoxError`
- `ConflictingPeriodicBoxesError`

## Takeaway

The Python binding stays small on purpose. Once you know `Target`, `Molpack`,
the built-in restraints, and the result diagnostics, you already know most of
what matters day to day.
