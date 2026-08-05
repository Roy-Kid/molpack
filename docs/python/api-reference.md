# API Reference

Import surface:

```python
from molpack import (
    # Core
    Target, Molpack, PackResult, StepInfo,
    # Typed values
    Angle, Axis, CenteringMode,
    # Geometric (per-atom) restraints
    InsideBoxRestraint, InsideSphereRestraint, OutsideSphereRestraint,
    AbovePlaneRestraint, BelowPlaneRestraint,
    # Collective (distribution-matching) restraints
    GaussianPlane, GaussianPoint,
    ExponentialPlane, ExponentialPoint,
    TabulatedPlane, TabulatedPoint,
    # In-loop relaxers
    TorsionMcRelaxer, LBFGSRelaxer,
    # Script loader (`.inp`)
    ScriptJob, load_script,
    # Parallel evaluation
    rayon_enabled, num_threads, init_thread_pool,
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

# Post-pack whole-system relaxation (optional molpy/LAMMPS backend)
from molpack.relaxer import LAMMPSRelaxer
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

- `frame` — a `molrs.Frame` or `molpy.Frame` with atom columns `"x"`,
  `"y"`, `"z"`, and `"element"` (or `"symbol"` for `molrs` PDB frames).
  Resolved zero-copy via its FFI capsule; a plain dict is not accepted.
- `count` — number of copies to produce.

**Builders**

- `.with_name(name: str)` — display label.
- `.with_restraint(r)` — attach a restraint to every atom (stackable).
  Accepts a geometric built-in, a collective (distribution-matching)
  restraint, or any duck-typed `f`/`fg` object — see [Restraints](#restraints).
- `.with_atom_restraint(indices: Sequence[int], r)` — 0-based indices.
- `.with_relaxer(relaxer)` — attach an in-loop geometry relaxer
  (`TorsionMcRelaxer` or `LBFGSRelaxer`); requires `count == 1`. See
  [In-loop relaxers](#in-loop-relaxers).
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
- `.with_periodic_box(min: [x,y,z], max: [x,y,z])` — declare a
  fully-periodic box directly on the packer (Packmol `pbc`). Alternative
  to a periodic `InsideBoxRestraint`; see
  [Periodic boundaries](guide/periodic-boundaries.md).
- `.with_perturb_fraction(f: float)` — stall perturbation fraction (default 0.05).
- `.with_random_perturb(enabled: bool)`.
- `.with_perturb(enabled: bool)` — master switch (default True).
- `.with_avoid_overlap(enabled: bool)` — reject initial random placements
  overlapping a fixed molecule (Packmol `avoid_overlap`; default True).
- `.with_seed(seed: int)` — deterministic RNG (default 0).
- `.with_parallel_eval(enabled: bool)` — rayon-backed pair eval. Raises
  `RuntimeError` if the wheel lacks the `rayon` feature (fail-fast).
- `.with_lammps_output(enabled: bool)` — enable LAMMPS-style screen output.
- `.with_log_level(level: str)` — `quiet`, `summary`, `progress`, or `verbose`.
- `.with_log_frequency(n: int)` — print every `n` outer steps.
- `.with_progress(enabled: bool)` — compatibility alias for progress output.
- `.with_handler(handler)` — attach a custom `Handler` (stackable).
- `.with_xyz_output(path: str, every: int = 1)` — record the packing
  trajectory to a multi-frame extended-XYZ file (a frame every `every`
  loops, loop 0 included).
- `.with_global_restraint(r)` — broadcast to every target (stackable).

**Run**

```python
.pack(targets: list[Target], max_loops: int = 200) -> molrs.Frame
.pack_with_report(targets: list[Target], max_loops: int = 200) -> PackResult
```

Raises a typed `PackError` subclass on failure.

---

## `PackResult`

Read-only output container returned by `pack_with_report()`.

**Properties**

- `.positions : ndarray (N, 3) float64`
- `.frame : molrs.Frame` — topology-complete frame (periodic box stamped if one was declared).
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

All restraint classes are immutable. Two families, both attached with
`target.with_restraint(r)` (or `Molpack.with_global_restraint(r)`):
**geometric** per-atom region restraints (below) and **collective**
distribution-matching restraints ([next section](#collective-distribution-matching-restraints)).

### Geometric (per-atom) restraints

Their `f`/`fg` see **one atom** at a time — a soft quadratic penalty that
is zero inside the region and rises outside.

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

### Collective (distribution-matching) restraints

Attached the same way (`target.with_restraint(r)`), but their `f`/`fg`
see **every copy** of the target at once and drive the species' spatial
distribution toward a target profile via a squared 1-D Wasserstein
(sorted-CDF) penalty. The reaction coordinate is either signed distance
to a plane ($\xi = \mathbf{x}\cdot\hat{\mathbf{n}} - \text{offset}$) or
radial distance to a point ($\xi = \lVert\mathbf{x} - \text{center}\rVert$).

| Class | Constructor | Target profile |
|-------|-------------|----------------|
| `GaussianPlane` | `(normal, offset, strength, mu, sigma)` | Gaussian $N(\mu, \sigma)$ slab |
| `GaussianPoint` | `(center, strength, mu, sigma)` | Gaussian shell (radius `mu`, thickness `sigma`) |
| `ExponentialPlane` | `(normal, offset, strength, lambda_)` | diffuse layer $\propto e^{-\xi/\lambda}$, $\xi \ge 0$ |
| `ExponentialPoint` | `(center, strength, lambda_)` | radial atmosphere $\propto e^{-\xi/\lambda}$ |
| `TabulatedPlane` | `(normal, offset, strength, xs, rho)` | arbitrary planar prior on grid `(xs, rho)` |
| `TabulatedPoint` | `(center, strength, xs, rho)` | arbitrary radial prior on grid `(xs, rho)` |

`strength` is the overall penalty multiplier $\lambda$. `sigma` /
`lambda_` must be `> 0`; tabulated `xs` must be strictly ascending
(≥ 2 points) with non-negative `rho` of positive total mass. Invalid
arguments raise `ValueError` at construction.

---

## In-loop relaxers

Relaxation-assisted packing: attach to a `Target` via
`target.with_relaxer(r)` to reshape a single molecule's reference
geometry *during* the pack loop. Both require the target's `count == 1`
(every copy shares the reference geometry the relaxer rewrites). Immutable.

### `TorsionMcRelaxer(frame)`

Monte-Carlo torsion-angle sampling — engine-free and force-field-free.
Rotatable bonds are detected from the frame's bond topology; proposed
rotations are accepted against the packer objective (Metropolis).

- `.with_temperature(t: float)` — Metropolis temperature (default 1.0).
- `.with_steps(n: int)` — MC steps proposed per packing iteration (default 10).
- `.with_max_delta(rad: float)` — max per-step rotation, radians (default π/6).
- `.with_self_avoidance(radius: float)` — quadratic overlap penalty on
  non-bonded intramolecular pairs closer than `2 * radius`; `0.0`
  disables (default).

```python
from molpack import Target, TorsionMcRelaxer

chain = TorsionMcRelaxer(frame).with_steps(20).with_self_avoidance(1.5)
target = Target(frame, count=1).with_relaxer(chain)
```

### `LBFGSRelaxer(forcefield)`  *(requires the `ff` feature)*

Force-field L-BFGS geometry minimization. Built from a
`molrs.ForceField` / `molpy.ForceField` (zero-copy FFI capsule); the
potential is compiled lazily against the molecule's frame when packing
starts.

- `.with_fmax(fmax: float)` — stop when the max per-atom force drops
  below `fmax` (kcal/mol/Å; default 0.05).
- `.with_max_steps(max_steps: int)` — L-BFGS iteration cap per
  relaxation call (default 500).

---

## Post-pack relaxation

### `LAMMPSRelaxer` *(from `molpack.relaxer`)*

Whole-system relaxation of a **finished** packed box via LAMMPS — a
different axis from the in-loop relaxers above. A thin façade over
`molpy.engine.LAMMPSEngine`, imported **lazily** so `import molpack`
never requires molpy.

```python
LAMMPSRelaxer(
    ff,                              # a typified molpy ForceField
    *,
    executable: str | None = None,   # None auto-detects lmp / lmp_serial / lmp_mpi
    launcher: list[str] | None = None,   # e.g. ["mpirun", "-np", "8"]
    pair_style: str = "lj/cut/coul/cut 10.0",
    atom_style: str = "full",
    units: str = "real",
    workdir: str | Path | None = None,
)
```

- `.minimize(target, **options) -> molrs.Frame` — energy minimisation.
- `.md(target, **options) -> molrs.Frame` — short MD settle.
- `.relax(target, **options)` — alias for `minimize` (also `__call__`).

`target` is a `PackResult` (its `.frame` is used) or a bare
`molrs.Frame` carrying a periodic box; the input is never mutated.
Raises `ImportError` if `molcrafts-molpy` is not installed.

---

## Script loader

### `load_script(path, *, read_frame=None) -> ScriptJob`

Parse and lower a Packmol-compatible `.inp` script. Template files are
read on the Python side (defaulting to `molrs.io.read_pdb` / `read_xyz` by
extension), so the wheel stays free of `molrs-io`. Pass `read_frame`
— a callable `(path, filetype) -> molrs.Frame` — to plug in another
loader (mdtraj, ASE, …).

### `ScriptJob`

Bundle returned by `load_script`. Access fields by attribute **or**
tuple-unpack it:

```python
job = load_script("mix.inp")
packer, targets, output, nloop = load_script("mix.inp")   # same object
```

- `.packer : Molpack` — pre-configured with the script's `tolerance`,
  `seed`, and any `pbc` box.
- `.targets : list[Target]`
- `.output : pathlib.Path` — resolved output path.
- `.nloop : int` — outer-loop cap (`nloop` keyword; default 400).

---

## Parallel evaluation

The parallel evaluator runs on rayon's process-global thread pool, built
**once** per process and not resizable afterwards.

- `rayon_enabled() -> bool` — was the wheel built with the `rayon` feature?
- `num_threads() -> int` — worker count the pool will use (1 on a serial build).
- `init_thread_pool(n: int)` — pin the pool size; must be called
  **before** the first pack. Raises `RuntimeError` without the `rayon`
  feature or if the pool was already initialized, and `ValueError` if
  `n == 0`. For a scaling sweep, set the count once per process (or via
  `RAYON_NUM_THREADS`) and launch one process per data point.

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
