# API Reference

Import surface:

```python
from molpack import (
    # Core
    Target, PackResult, StepInfo,
    # Engine entries — one per packing algorithm
    GenCanPack, CbmcGrow,
    # Typed values
    Angle, Axis, CenteringMode,
    # Geometric (per-atom) restraints
    InsideBoxRestraint, InsideSphereRestraint, OutsideSphereRestraint,
    AbovePlaneRestraint, BelowPlaneRestraint,
    # Collective (distribution-matching) restraints
    GaussianPlane, GaussianPoint,
    ExponentialPlane, ExponentialPoint,
    TabulatedPlane, TabulatedPoint,
    # Chain-growth priors
    TorsionPrior, AnglePrior,
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
- `.with_perturb_budget(n: int)` — per-target perturbation budget.
- `.with_centering(mode: CenteringMode)`.
- `.with_rotation_bound(axis: Axis, center: Angle, half_width: Angle)`.
- `.fixed_at(position: [x, y, z])` — pin the target.
- `.with_orientation((ax, ay, az))` — Euler tuple of `Angle`s; must
  follow `fixed_at`.
- `.with_mass(amu: float)` — per-copy total mass override for
  `with_density` when element symbols cannot provide one (CG
  beads, bare-coordinate targets).

**Static constructor**

- `Target.fixed_from(result: PackResult)` — wrap a previous run's whole
  output as one fixed obstacle, coordinates kept verbatim. The chaining
  primitive for staged packs: grow with `CbmcGrow`, then pack the next
  species around the frozen matrix with `GenCanPack`. See
  [Engine entries](#engine-entries).

**Properties**

- `.name : str | None`
- `.natoms : int`
- `.count : int`
- `.elements : list[str]`
- `.radii : list[float]`
- `.is_fixed : bool`

---

## Engine entries

One entry class per packing algorithm; you pick the algorithm by picking
the entry. Both are immutable builders — every `with_*` returns a new
instance — and both expose the same terminal verb:

```python
.run(targets: list[Target], max_loops: int) -> PackResult
```

`run()` consumes the entry (one engine, one run): a second call on the
same object raises `RuntimeError`. Raises a typed `PackError` subclass on
packing failure.

### Shared builders

Available on `GenCanPack` **and** `CbmcGrow`:

- `.with_tolerance(t: float)` — minimum pairwise distance (Å; default 2.0).
- `.with_precision(p: float)` — convergence threshold (default 0.01).
- `.with_seed(seed: int)` — deterministic RNG (default Packmol's 1234567).
- `.with_periodic_box(min: [x,y,z], max: [x,y,z])` — declare a
  fully-periodic box directly on the entry (Packmol `pbc`). Alternative
  to a periodic `InsideBoxRestraint`; see
  [Periodic boundaries](guide/periodic-boundaries.md).
- `.with_density(rho: float)` — size the box from a target mass density
  (g/cm³) instead of declaring it: `run()` resolves a cubic, fully
  periodic box holding the total mass of all targets at `rho`. Masses
  come from element symbols or `Target.with_mass`; an unresolvable mass
  raises `ValueError`. Mutually exclusive with `.with_periodic_box`.
- `.with_parallel_eval(enabled: bool)` — rayon-backed pair eval. Raises
  `RuntimeError` if the wheel lacks the `rayon` feature (fail-fast).
- `.with_progress(on: bool = True)` — LAMMPS-style screen output
  (off by default).
- `.with_handler(handler)` — attach a custom `Handler` (stackable).
- `.with_global_restraint(r)` — broadcast to every target (stackable).

### `GenCanPack`

Rigid-body placement driven by the three-phase GENCAN optimizer — the
Packmol algorithm. Zero-arg constructor.

```python
GenCanPack()
```

GENCAN-only builders, on top of the shared ones:

- `.seeded_from(result: PackResult)` — continue on a previous run's
  placement solution (the explicit push-off chain): the free copies start
  exactly where `result` left them, `initial()` is skipped, and the
  stall-perturbation moves stay off — molecules are pushed apart by
  rigid-body descent only. The cell travels with the seed; declaring a
  box, density, or cell on a seeded engine is a named error, and the run's
  free targets must match the seed's shape.
- `.with_inner_iterations(n: int)` — GENCAN inner-loop cap (default 20).
- `.with_init_passes(n: int)` — init compaction passes (0 = auto).
- `.with_init_box_half_size(h: float)` — init placement bound (default 1000 Å).
- `.with_perturb(fraction: float, random: bool = False, enabled: bool = True)`
  — stall-perturbation heuristic: fraction re-sampled per stall, random
  vs worst-first selection, and the master switch.
- `.with_avoid_overlap(on: bool = True)` — reject initial random placements
  overlapping a fixed molecule (Packmol `avoid_overlap`; default True).

### `CbmcGrow`

Configurational-bias chain growth for dense melts — see the
[Chain growth guide](guide/growth.md). The torsion prior is the one
mandatory constructor argument, because it decides the chain statistics
of the product.

```python
CbmcGrow(torsion_prior: TorsionPrior)
```

Growth-only builders, on top of the shared ones:

- `.with_trials(k: int)` — torsion candidates per growth step (CBMC
  `k`; clamped ≥ 1, default 12).
- `.with_selectivity(beta: float)` — Rosenbluth inverse temperature
  applied to the soft-shell crowding penalty (clamped ≥ 0, default 2.0).
- `.with_soft_shell(width: float)` — soft-shell width in Å beyond the
  hard core: allowed but charged (clamped ≥ 0, default 1.0).
- `.with_retract(steps: int)` — steps retracted on a dead end; repeated
  dead ends retract exponentially deeper (clamped ≥ 1, default 10).
- `.with_relax(every: int, window: int)` — regrow each chain's last
  `window` steps every `every` rounds; `every=0` disables (default
  `(25, 6)`). A regrown tail is only kept when its Rosenbluth weight
  does not degrade.
- `.with_soften_after(attempts: int)` — consecutive dead ends at one
  step before the hard core softens (clamped ≥ 1, default 50).
- `.with_min_hard_scale(scale: float)` — softening floor, clamped to
  `[0, 1]` (default 0.8, the classic push-off bound).
- `.with_exclusion_depth(depth: int)` — intramolecular exclusion depth
  in bonds: 3 (default) is the all-atom 1-2/1-3/1-4 convention; CG
  templates conventionally use 1 or 2.
- `.with_angle_prior(prior: AnglePrior)` — placement-angle prior
  (default `AnglePrior.template()`).
- `.with_serial(serial: bool = True)` — grow one chain to completion
  before starting the next.
- `.with_void_bias(void_bias: bool = True)` — seed chains in empty field
  cells (cavity seeding).

### `LatticeGrow`

Diamond-lattice growth for melt density and above (see the
[Chain growth guide](guide/growth.md)): the backbone grows as an on-lattice
self-avoiding walk with RIS weights, then decorates back to the template's
exact bonded geometry. Same mandatory torsion-prior constructor and shared
builders; one extra knob:

- `.with_occupancy_guard(on: bool = True)` — nearest-neighbour site
  exclusion (keeps non-bonded pairs at ≥ the 2nd-neighbour distance).

Linear sp³ heavy-atom backbones only in v1 — branched and non-tetrahedral
templates are named rejections.

### Chaining two entries

There is no mixed-algorithm pack inside one `run()`, and no hidden
fallback. Stage it in user code, in one of two shapes:

```python
# Push-off: continue the SAME free targets on the grown state.
grown = CbmcGrow(prior).with_density(0.9).run([chain], max_loops=60)
pushed = GenCanPack().seeded_from(grown).with_seed(7).run([chain], max_loops=60)

# Fixed matrix: freeze the first result, pack new species around it.
full = GenCanPack().run([Target.fixed_from(grown), solvent], max_loops=200)
```

---

## `PackResult`

Read-only output container returned by `run()`.

**Properties**

- `.positions : ndarray (N, 3) float64`
- `.frame : molrs.Frame` — topology-complete frame (periodic box stamped if one was declared).
- `.elements : list[str]`
- `.natoms : int`
- `.converged : bool`
- `.fdist : float`
- `.frest : float`
- `.softened : int` — how many times the growth solver had to relax its
  constructive hard-core guarantee (always 0 on the GENCAN path). See
  [`CbmcGrow`](#cbmcgrow).

---

## Chain-growth priors

Inputs to [`CbmcGrow`](#cbmcgrow): the torsion prior is its mandatory
constructor argument, the angle prior an optional builder. Both are
geometric data — see the [Chain growth guide](guide/growth.md) for when
and why each matters.

### `TorsionPrior`

Geometric data only — never a force field. Static constructors:

- `TorsionPrior.uniform()` — uniform on (−π, π]. Freely-rotating-chain
  statistics (C∞ = 2.0); negative control only, quantitatively wrong for
  melts.
- `TorsionPrior.template(kappa: float)` — von-Mises-like spread of
  concentration `kappa` around the template's own torsion values.
- `TorsionPrior.states(states: list[tuple[float, float]])` — RIS-style
  discrete `(angle_rad, weight)` states; weights are normalized at use.
- `TorsionPrior.three_state_from_c_inf(c_inf: float, theta_rad: float)`
  — trans/gauche± prior whose trans fraction is solved from a target
  characteristic ratio. PEO with tetrahedral backbone angles:
  `TorsionPrior.three_state_from_c_inf(5.5, 1.9106)`.

### `AnglePrior`

- `AnglePrior.template()` — copy bond angles verbatim from the template
  (the all-atom default).
- `AnglePrior.wlc(kappa: float)` — discrete worm-like chain with tilt
  `kappa` (CG persistence control).
- `AnglePrior.wlc_from_c_inf(c_inf: float)` — WLC tilt calibrated from
  a target characteristic ratio; assumes uniform torsions.
  Kremer–Grest melts: `AnglePrior.wlc_from_c_inf(1.76)`.

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
`target.with_restraint(r)` (or the entry's `with_global_restraint(r)`):
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

- `.packer : GenCanPack` — pre-configured with the script's `tolerance`,
  `seed`, and any `pbc` box. `.inp` scripts always lower to the
  rigid-body entry.
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
(bad atom indices, wrong restraint object, etc.). Growth and density
declarations are input contracts, so their failures also raise
`ValueError`: a target that cannot be grown (no bond graph, < 3 atoms,
fixed placement, no box), a `with_density` fighting an explicit box, or
a mass the element symbols cannot resolve.
