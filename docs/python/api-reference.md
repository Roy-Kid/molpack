# API Reference

Import surface:

```python
from molpack import (
    # Core
    Target, State, IntraResidual, StepInfo,
    # Engine entries — one per packing algorithm
    GenCanPack, CbmcGrow, LatticeGrow,
    # Multi-stage composition
    Pipeline,
    # Typed values
    Angle, Axis, CenteringMode,
    # (geometric restraints are molrs regions — see "molrs regions as restraints")
    # Collective (distribution-matching) restraints
    GaussianPlane, GaussianPoint,
    ExponentialPlane, ExponentialPoint,
    TabulatedPlane, TabulatedPoint,
    # Collective (pairwise separation) restraint
    SelfSeparation,
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

- `frame` — a `molrs.core.Frame` (`molpy.Frame` is the same class) with atom
  columns `"x"`, `"y"`, `"z"`, and `"element"`.
  Resolved zero-copy via its FFI capsule; a plain dict is not accepted.
- `count` — number of copies to produce.

**Builders**

- `.with_name(name: str)` — display label.
- `.with_restraint(r)` — attach a restraint to every atom (stackable).
  Accepts a geometric built-in, a collective (distribution-matching)
  restraint, or any duck-typed `f`/`fg` object — see [Restraints](#restraints).
- `.with_atom_restraint(indices: Sequence[int], r)` — 0-based indices.
- `.with_radius(radius: float)` — packing radius for every atom
  (Packmol `radius`). Raises `ValueError` if not positive.
- `.with_atom_radius(indices, radius)` — 0-based; all-atom explicit
  hydrogen keeps the default skip table and uses ~0.85 Å here.
- `.with_fscale(fscale: float)` / `.with_atom_fscale(indices, fscale)` —
  overlap-penalty weight (Packmol `fscale`; default 1.0).
- `.with_short_radius(short_radius: float)` /
  `.with_atom_short_radius(indices, short_radius)` — second, shorter
  penalty radius (Packmol `short_radius`).
- `.with_short_radius_scale(scale: float)` /
  `.with_atom_short_radius_scale(indices, scale)` — weight of the
  short-radius penalty.
- `.with_special_bonds(table: Sequence[float])` — intramolecular skip
  weights. Slot 0 is 1-2; the last slot is the 1-N tail. Default
  `[0, 0, 0, 1]` (depth 3). Empty / non-finite / outside `[0, 1]` raise
  `ValueError`. Fractional `0.5` is stored and refused at
  `CbmcGrow.run`. This is not a force-field `special_bonds` triple.
- `.with_hydrogens(indices: Sequence[int])` — the atoms lattice growth
  treats as hydrogens (0-based), placed off their backbone neighbour rather
  than on a lattice site. Default: element symbol `H`. Pass `[]` for a
  coarse-grained model. An out-of-range index raises `ValueError`.
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

- `Target.fixed_from(result: State)` — wrap a previous run's whole
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
- `.special_bonds : list[float]` — intramolecular skip table; default
  `[0.0, 0.0, 0.0, 1.0]`.
- `.is_fixed : bool`

---

## Engine entries

One entry class per packing algorithm; you pick the algorithm by picking
the entry. All three are immutable builders — every `with_*` returns a new
instance — and all three expose the same terminal verb:

```python
.run(targets: list[Target], max_loops: int) -> State
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
  fully-periodic box directly on the entry (Packmol `pbc`); see
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

- `.with_restart(result: State)` — continue on a previous run's
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
- `.with_soften_after(attempts: int)` — cumulative dead ends on one
  chain before that chain's hard core softens by one rung (×0.97); a
  successful placement does not reset the count (clamped ≥ 1, default 50).
  Softening is per chain: a rung on one chain leaves every other chain at
  full contact.
- `.with_min_hard_scale(scale: float)` — softening floor, clamped to
  `[0, 1]` (default 0.8, the classic push-off bound).
- `.with_angle_prior(prior: AnglePrior)` — placement-angle prior
  (default `AnglePrior.template()`).
- `.with_serial(serial: bool = True)` — grow one chain to completion
  before starting the next.
- `.with_void_bias(void_bias: bool = True)` — seed chains in empty field
  cells (cavity seeding).

### `LatticeGrow`

Diamond-lattice growth for melt density and above (see the
[Chain growth guide](guide/growth.md)): the backbone grows as an on-lattice
self-avoiding walk with RIS weights, and each backbone atom is then seated on
the site the walk chose. **A template supplies its topology, not its
geometry**: the backbone's bond lengths and angles are the lattice's — one
step throughout, sized from the template's own mean backbone bond and moved a
percent or two by how the cell divides — and its torsions are exactly the
trans/gauche± the prior drew. Only hydrogens and side atoms keep the
template's local geometry. Same mandatory torsion-prior constructor and shared
builders; one extra knob:

- `.with_occupancy_guard(on: bool = True)` — nearest-neighbour site
  exclusion (keeps non-bonded pairs at ≥ the 2nd-neighbour distance).

Trees (including branched) are accepted; a cycle raises named
`RingTemplate` `ValueError`; non-tetrahedral templates stay named
rejections.

### `Pipeline`

Composes stage objects — `GenCanPack`, `CbmcGrow`, and `LatticeGrow`
instances — into one multi-algorithm run: grow a chain, then hand it to
rigid-body descent, in one lifecycle and one `State` rather than two
separate `run()` calls. See [Composing stages](guide/packer.md#composing-stages)
for the full walkthrough.

```python
Pipeline(stages: Sequence[GenCanPack | CbmcGrow | LatticeGrow] | None = None)
```

**Builders**

- `.with_stage(stage: GenCanPack | CbmcGrow | LatticeGrow)` — append one
  more stage; returns a new `Pipeline`.
- The shared builders — same as [`GenCanPack`](#shared-builders):
  `.with_tolerance`, `.with_precision`, `.with_seed`,
  `.with_periodic_box`, `.with_density`, `.with_parallel_eval`,
  `.with_progress`, `.with_handler`, `.with_global_restraint`. Set these on
  the `Pipeline`, never on a stage object that goes into one — a stage
  carrying a non-default shared setting raises `ValueError` naming the
  stage and the setting.

**Running**

```python
.run(targets: list[Target], max_loops: int) -> State
```

Each stage's own `.with_handler(...)` callbacks are adopted into the
pipeline's handler set and fire for every stage in the run, not only the
one they were attached to. Raises `ValueError` for an empty pipeline, a
stage-ordering error, or a stage carrying a non-default shared setting
(each naming the offending stage); raises `TypeError`, listing the three
supported entries, for any object passed to `Pipeline([...])` or
`.with_stage(x)` that is not a `GenCanPack`, `CbmcGrow`, or `LatticeGrow`.

### Chaining two entries

An entry's own `run()` never mixes algorithms, and there is no hidden
fallback between them. Outside `Pipeline` (above), stage it in user code
instead, in one of two shapes:

```python
# Push-off: continue the SAME free targets on the grown state.
grown = CbmcGrow(prior).with_density(0.9).run([chain], max_loops=60)
pushed = GenCanPack().with_restart(grown).with_seed(7).run([chain], max_loops=60)

# Fixed matrix: freeze the first result, pack new species around it.
full = GenCanPack().run([Target.fixed_from(grown), solvent], max_loops=200)
```

---

## `State`

Frozen outcome of one `run()`. Diagnostics (`frame`, `fdist`, `frest`,
`converged`, `degraded`, `intra`) live on this object; pass the same
object to `GenCanPack.with_restart` or `Target.fixed_from` to continue.

**Properties**

- `.positions : ndarray (N, 3) float64`
- `.frame : molrs.core.Frame` — topology-complete frame (periodic box stamped if one was declared).
- `.elements : list[str]`
- `.natoms : int`
- `.converged : bool`
- `.fdist : float`
- `.frest : float`
- `.degraded : int` — how many molecules were placed below the growth
  solver's own guarantee, one count per demotion. `LatticeGrow` demotes a chain
  when the walk cannot keep the occupancy guard (no two non-bonded atoms closer
  than the 2nd lattice neighbour): the chain goes in with site self-avoidance
  only, or as a forced zigzag. `CbmcGrow` counts each hard-core softening rung
  the same way. `0` on the GENCAN path, which promises nothing constructively.
  A non-zero count is the honest reading of a crowded box — those molecules
  carry the close contacts `fdist` reports, and `converged` is false while it
  stands.
- `.intra : IntraResidual` — same-copy scored vs exempted minima (Å,
  minimum image). Forwards the assembled residual; does not recompute
  from positions.

### `IntraResidual`

Nested diagnostic on [`State`](#state). Empty class is `+∞`.

- `.scored : float` — minimum same-copy pair distance among pairs the
  target's skip table scores (Å).
- `.exempted : float` — minimum same-copy pair distance among pairs the
  table exempts (Å).

There are no `min_intra_*` aliases.

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
info.stage             # StageInfo — which packing algorithm emitted this step
```

### `StageInfo`

Read-only triple identifying the stage a step belongs to. Load-bearing
inside a multi-stage [`Pipeline`](#pipeline); present, with `index = 0` and
`total = 1`, on every single-entry run too.

- `.index : int` — 0-based position of this stage in the run; monotonic
  across a multi-stage `Pipeline`.
- `.total : int` — number of stages in the run.
- `.name : str` — the stage's own name: `"gencan"`, `"growth"`, or
  `"lattice"`.

---

## Restraints

All restraint classes are immutable. Two families, both attached with
`target.with_restraint(r)` (or the entry's `with_global_restraint(r)`):
**molrs regions** lifted to a per-atom penalty (below) and **collective**
distribution-matching restraints ([next section](#collective-distribution-matching-restraints)).

### molrs regions as restraints

Any molrs region object attaches as a restraint: `molrs.core.Sphere`, `Cuboid`,
`Parallelepiped`, `HalfSpace`, `Cylinder`, `Ellipsoid`, `Polyhedron`,
`SphereUnion`, or any `&` / `|` / `~` composition of them. The region
crosses the wheel boundary as a `molrs.RegionRef` capsule (both wheels on
one molrs minor line) and is lifted to `scale · max(0, distance)²` — a
soft quadratic penalty on the **atom centre** that is zero inside the
region and on its boundary. molpack defines no geometric restraint class;
the shapes, their constructors and their `contains` / `distance` /
`bounds` queries are documented with molrs. Solver split (the same for
every region, not inferred from the shape):

- `GenCanPack` — soft exterior penalty.
- `CbmcGrow` — hard reject on propose; `force_place` may sit outside.
- `LatticeGrow` — sites outside the mesh are blocked (Region ∩ lattice), and
  the backbone atoms **are** those sites, so the mask's guarantee is the
  molecule's. Hydrogens and side atoms hang off the backbone with the
  template's local geometry and can reach about a bond length past the
  surface; author the mesh with that clearance in it if the wall has to hold.

The region answers its own questions for an `(n, 3)` array: `contains`
returns `(n,)` bool; `distance` returns `(n,)` Å, negative inside. Together
they say what the packer was told to enforce:

```python
cavity = molrs.core.Polyhedron(molrs.io.read_stl("dendrite.stl"))
depth = cavity.distance(state.positions)
print(f"{(depth > 0).sum()} atoms outside, worst {depth.max():.2f} Å")
```

A restraint is satisfied to the run's `precision`, not exactly. `frest` is the
largest per-atom penalty `0.01 · d²`, so `frest < precision` means
`d < 10·√precision`: the default `precision=1e-2` calls a run converged with an
atom 1 Å outside, `1e-4` with 0.1 Å. Tighten `with_precision` when the wall is
the point — and remember the mesh, not the solver, is where clearance belongs.

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

### Collective (pairwise separation) restraint

| Class | Constructor | Meaning |
|-------|-------------|---------|
| `SelfSeparation` | `(d_min, strength=1.0)` | no two copies of this species closer than `d_min`, centre to centre |

The anti-clustering restraint. A distribution target says *where* copies
should be; this says how close two of them may come. The pair term keeps
atoms from overlapping and then stops caring, and it does not distinguish
two copies of one species from a copy of each of two — so without this a
species may pile its copies into one corner.

The penalty is silent above `d_min` and grows as $(d_{min} - D)^2$ below it,
where $D$ is the minimum-image distance between two copies' geometric
centroids. It is quadratic in the length by which the bound is missed — the
same shape as a geometric restraint's penalty — so `strength=1.0` weights a
shortfall like a region lift weights an equal overshoot, and the
convergence threshold `precision` means the same thing for both.

Centroids, not atoms: for a compact molecule the centroid stands in for the
whole, but two long chains can interdigitate with distant centroids. The
restraint states what it measures.

Cost is linear in `count`, not quadratic: the centres are binned into a cell
grid and only nearby pairs are examined, so the bound stays affordable on a
melt-sized species.

Feasibility is not checked, and does not need to be. A bound is either met or
not, so unlike the distribution restraints this one counts toward the restraint
verdict `frest` and gates convergence: if `count` copies cannot fit at `d_min`,
the run does not converge and says so — it is never silently relaxed. `d_min`
and `strength` must be `> 0`; otherwise `ValueError`.

```python
import molrs
from molpack import SelfSeparation, Target

ions = (
    Target(frame, count=27)
    .with_restraint(molrs.core.Cuboid([0, 0, 0], [40, 40, 40]))
    .with_restraint(SelfSeparation(10.0))
)
```

---

## Script loader

### `load_script(path, *, read_frame=None) -> ScriptJob`

Parse and lower a Packmol-compatible `.inp` script. Template files are
read on the Python side (defaulting to the `molrs.io` reader of the format the
script's `filetype` or the file name names: `read_pdb`, `read_xyz`, …), so the
wheel stays free of `molrs-io`. Pass `read_frame`
— a callable `(path, filetype) -> molrs.core.Frame` — to plug in another
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

`ValueError` / `TypeError` still surface on Python-side invariants
(bad atom indices, wrong restraint object, etc.). Growth and density
declarations are input contracts, so their failures also raise
`ValueError`: a target that cannot be grown (no bond graph, < 3 atoms,
fixed placement, no box), a `with_density` fighting an explicit box, or
a mass the element symbols cannot resolve.
