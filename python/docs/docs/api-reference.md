# API Reference

Import surface:

```python
from molpack import (
    # Core
    Target, Molpack, PackResult, StepInfo,
    # Typed values
    Angle, Axis, CenteringMode,
    # Restraints
    InsideBoxRestraint, InsideSphereRestraint, OutsideSphereRestraint,
    AbovePlaneRestraint, BelowPlaneRestraint,
    # Protocols
    Handler, Restraint,
    # Errors
    PackError,
    ConstraintsFailedError,
    MaxIterationsError,
    NoTargetsError,
    EmptyMoleculeError,
    InvalidPBCBoxError,
    ConflictingPeriodicBoxesError,
)
```

---

## `Angle`

Angular quantity with explicit units at the call site.

```python
Angle.from_degrees(30.0).radians   # 0.5235...
Angle.from_radians(0.5).degrees    # 28.6...
Angle.ZERO                         # identity rotation
```

---

## `Axis`

Cartesian axis enum: `Axis.X`, `Axis.Y`, `Axis.Z`.

---

## `CenteringMode`

Centering policy for target reference coords:

- `CenteringMode.AUTO` — free targets centered, fixed kept in place (default).
- `CenteringMode.CENTER` — always center.
- `CenteringMode.OFF` — keep input coords unchanged.

---

## `Target`

Molecule-type specification. Immutable — builder methods return new
instances.

**Constructor**

```python
Target(frame, count: int)
```

- `frame` — any object with `frame["atoms"]` returning a block with
  columns `"x"`, `"y"`, `"z"`, and `"element"` (or `"symbol"` for
  `molrs` PDB frames). Accepts `molrs.Frame`, `molpy.Frame`, or a
  plain dict.
- `count` — number of copies to produce.

**Builders**

- `.with_name(name: str)` — display label.
- `.with_restraint(r)` — attach to all atoms (stackable).
- `.with_atom_restraint(indices: Sequence[int], r)` — 0-based indices.
- `.with_perturb_budget(n: int)` — per-target perturbation budget.
- `.with_centering(mode: CenteringMode)`.
- `.with_rotation_bound(axis: Axis, center: Angle, half_width: Angle)`.
- `.fixed_at(position: [x, y, z])` — pin the target.
- `.with_orientation((ax, ay, az))` — Euler tuple of `Angle`s; must
  follow `fixed_at`.

**Properties**

- `.name : str | None`
- `.natoms : int`
- `.count : int`
- `.elements : list[str]`
- `.radii : list[float]`
- `.is_fixed : bool`

---

## `Molpack`

Orchestrator for the three-phase GENCAN optimizer. Zero-arg
constructor — all tuning is via `with_*` builders.

**Constructor**

```python
Molpack()
```

**Builders**

- `.with_tolerance(t: float)` — minimum pairwise distance (Å; default 2.0).
- `.with_precision(p: float)` — convergence threshold (default 0.01).
- `.with_inner_iterations(n: int)` — GENCAN inner-loop cap (default 20).
- `.with_init_passes(n: int)` — init compaction passes (0 = auto).
- `.with_init_box_half_size(h: float)` — init placement bound (default 1000 Å).
- `.with_perturb_fraction(f: float)` — stall perturbation fraction (default 0.05).
- `.with_random_perturb(enabled: bool)`.
- `.with_perturb(enabled: bool)` — master switch (default True).
- `.with_seed(seed: int)` — deterministic RNG (default 0).
- `.with_parallel_eval(enabled: bool)` — rayon-backed pair eval.
- `.with_progress(enabled: bool)` — built-in progress reporter.
- `.with_handler(handler)` — attach a custom `Handler` (stackable).
- `.with_global_restraint(r)` — broadcast to every target (stackable).

**Run**

```python
.pack(targets: list[Target], max_loops: int = 200) -> PackResult
```

Raises a typed `PackError` subclass on failure.

---

## `PackResult`

Read-only output container.

**Properties**

- `.positions : ndarray (N, 3) float64`
- `.frame : dict` — Frame-compatible dict with an `atoms` block.
- `.elements : list[str]`
- `.natoms : int`
- `.converged : bool`
- `.fdist : float`
- `.frest : float`

---

## `StepInfo`

Read-only snapshot passed to `Handler.on_step`.

```python
info.loop_idx          # outer-loop iteration
info.max_loops
info.phase             # phase index
info.total_phases
info.molecule_type     # int | None
info.fdist
info.frest
info.improvement_pct
info.radscale
info.precision
info.relaxer_acceptance  # list[tuple[int, float]]
```

---

## Restraints

All restraint classes are immutable.

### `InsideBoxRestraint(min, max, periodic=(False, False, False))`

Axis-aligned box. `periodic` is a 3-tuple of booleans declaring per-axis
periodicity — see [Periodic boundaries](guide/periodic-boundaries.md).

### `InsideSphereRestraint(center, radius)`

Closed ball.

### `OutsideSphereRestraint(center, radius)`

Complement of closed ball.

### `AbovePlaneRestraint(normal, distance)`

Half-space $\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} \ge d\}$.

### `BelowPlaneRestraint(normal, distance)`

Half-space $\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} \le d\}$.

---

## Duck-type protocols

### `Restraint`

```python
class Restraint(Protocol):
    def f(self, x: tuple[float, float, float], scale: float, scale2: float) -> float: ...
    def fg(
        self, x: tuple[float, float, float], scale: float, scale2: float,
    ) -> tuple[float, tuple[float, float, float]]: ...
```

### `Handler`

```python
class Handler(Protocol):
    def on_start(self, ntotat: int, ntotmol: int) -> None: ...
    def on_step(self, info: StepInfo) -> bool | None: ...   # True → stop
    def on_finish(self) -> None: ...
```

All `Handler` methods are optional — missing ones are silently skipped.

---

## Exceptions

Typed hierarchy rooted at `PackError` (itself a `RuntimeError`
subclass). Catch the base to handle any packing failure uniformly.

- `PackError` — base.
- `ConstraintsFailedError` — solver could not satisfy restraints even
  without distance tolerances.
- `MaxIterationsError` — ran out of outer loops.
- `NoTargetsError` — empty target list.
- `EmptyMoleculeError` — a target has zero atoms.
- `InvalidPBCBoxError` — periodic box has a non-positive extent.
- `ConflictingPeriodicBoxesError` — two restraints declared
  incompatible periodic boxes.

`ValueError` / `TypeError` still surface on Python-side invariants
(bad atom indices, wrong restraint object, etc.).
