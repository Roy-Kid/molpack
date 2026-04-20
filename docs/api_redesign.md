# molpack API redesign — spec

Status: **draft, pending approval**
Author: Roy Kid
Date: 2026-04-17

This document nails down the target public-API surface for `molpack` and the
plan to get there. Scope is **Rust naming and ergonomics only** — no changes
to numerics, to Packmol parity in the internals, or to performance-relevant
dispatch. The Python binding is **out of scope** for this spec; a follow-up
pass will mirror the final Rust surface into `python/src/`.

The goal is a single consistent style so downstream users don't have to
learn three builder dialects to pack a solvated protein.

The plan is one mechanical rename pass applied in staged PRs, with the
outcome validated by the existing five canonical examples still converging
bit-for-bit against the Packmol reference.

---

## 1. Motivation

The current public surface grew incrementally during the Fortran port and
exhibits several styles simultaneously:

- `Molpack` uses **bare** setters (`.precision(f)`) plus **`add_*`** for
  collections.
- `Target` uses **`with_*`** for most setters and collections, **bare** for
  `fixed_at`/`constrain_rotation_*`, and **`without_*`** for one negation
  (`without_centering`).
- `EarlyStopHandler` and `TorsionMcRelaxer` use **`with_*`** consistently.

Calling code mixes all three in one expression:

```rust
let mut packer = Molpack::new()
    .precision(0.01)                // bare
    .tolerance(2.0)                 // bare
    .add_handler(ProgressHandler::new());  // add_*

let water = Target::new(water, 1000)
    .with_name("water")             // with_*
    .with_restraint(sphere)         // with_*
    .fixed_at([0.0, 0.0, 0.0]);     // bare, and silently mutates count=1
```

On top of that, several method names are Packmol Fortran jargon
(`maxit`, `nloop0`, `sidemax`, `movefrac`, `disable_movebad`, `discale`,
`pbc`) that a Rust user encountering this crate cold has no way to
interpret.

## 2. Scope

### In scope

- Rust public surface of the `molpack` crate: all items re-exported from
  `src/lib.rs`.
- All Rust examples, tests, benches, and rustdoc under the crate.
- Adding a `molpack::prelude` module (see §5.8).

### Out of scope

- **Python binding surface (`python/src/*.rs`).** Updated in a separate
  follow-up spec once the Rust side stabilizes.
- Numerical behavior, Packmol parity, GENCAN internals.
- `PackContext` internal fields (still `pub` for test access but not part
  of the supported user surface).
- `Molpack::validate_from_targets` and `ExampleCase` helpers — they are
  already well-named; they may pick up cosmetic signature touches but no
  structural changes.

### Non-goals

- No new features. No new restraints, no new handlers, no new solver
  options.
- No deprecation of numerically-meaningful knobs — every semantic stays,
  only the name changes.

## 3. Principles

Four rules drive every rename in §5. Everything here flows from these.

### P1 — one builder verb: `with_*`, and `with_*` means *optional*

Every builder method that returns `Self` is named `with_<noun>`, whether
it sets a scalar or appends to a collection. Collection methods clarify
the append semantic in one line of rustdoc.

**`with_*` is reserved for arguments that are genuinely optional** — the
library has a defensible default when the user doesn't call the method.
Anything that has no sensible default belongs in `new()` (or in a
`from_<source>` constructor), not behind `with_*`. The negative test:
if omitting the method would leave the builder in an unusable state, it
is not a `with_*`.

Concrete consequence for current code: `Molpack::pack(targets, max_loops, seed)`
splits — `seed` has a defensible default (0 → deterministic) and
becomes `.with_seed(u64)`; `max_loops` has no reasonable default (the
right value depends on system size and convergence difficulty) and
stays as a positional argument to the terminal `.pack(&targets, max_loops)`
call.

Rationale: `with_*` is the dominant Rust-ecosystem builder verb
(`reqwest`, `tokio::runtime`, `clap`, `axum`). Two verbs (`with_*` +
`add_*`) earn nothing except cognitive friction. Bare setters look like
mutators to Rust readers. Reserving `with_*` for "optional, has a
default" gives users a reliable reading: if it's `with_*`, you can skip
it; if it's in `new()` or a positional arg, you must supply it.

### P2 — no Packmol jargon in user-facing names

If a name appears only in `getinp.f90` or `packmol.f90`, it does not
appear in the Rust public API. Every public method/field reads as
English describing *what it controls*, not *what the Fortran call is*.

Rationale: molpack's audience is Rust/Python users doing molecular
packing, not Packmol power users reading a reference port. The
`docs/packmol_parity.md` page already maps Rust names to Fortran kind
numbers for people who need that bridge; the API itself should not
require it.

Exceptions (stay as-is):

- `fdist` / `frest` on `PackResult` — these appear in scientific papers
  and plots; renaming to `max_overlap` / `max_constraint_violation` is
  on the table but has a much larger surface and is deferred.
- `GENCAN` as a module name is a specific Birgin-Martínez solver, not
  generic jargon; it stays.

### P3 — no double negatives, no hidden side effects

- No method named `disable_X(bool)`; the shape `with_X(bool)` always
  reads "set X to the argument."
- No method that quietly mutates an unrelated field (`fixed_at` silently
  setting `count = 1`). Either the method asserts preconditions loudly,
  or the feature gets its own constructor.

### P4 — one method per concept, not three siblings

`constrain_rotation_x`, `constrain_rotation_y`, `constrain_rotation_z`
collapse into one `with_rotation_bound(axis, ...)`. Same for
`pbc` / `pbc_box`. Multiple entry points only where semantics genuinely
differ (`Target::new` from Frame vs `Target::from_coords` from arrays —
both stay).

## 4. Naming conventions

Summary reference; §5 applies these.

| Kind | Convention | Example |
|---|---|---|
| Builder that returns `Self` (optional) | `with_<noun>` | `Molpack::with_precision(f)` |
| Required construction argument | positional arg of `new` / `from_<source>` | `Target::new(frame, n)` |
| Terminal that consumes builder | verb | `Molpack::pack(&targets, max_loops)` |
| Handler callback | `on_<event>` past-tense | `on_initialized`, `on_step` |
| Predicate / query | `is_<adj>` / bare noun | `is_parallel_safe`, `acceptance_rate` |
| Angle arguments | `Angle` newtype | `Angle::from_degrees(30.0)` |
| Boolean switches | `with_<feature>(bool)` — `false` disables | `.with_perturb(false)` |
| Collection growth | `with_<singular>` | `.with_handler(h)` appends |
| Enum "no-op" variant | `Off` (not `None`) | `CenteringMode::Off` |

## 5. Target surface — rename tables

Every row is a concrete edit. Columns: **Current** · **Proposed** · **Kind
of change** (R = rename / M = merge / S = signature / D = delete /
N = new) · **Why**.

### 5.1 `Molpack`

**The type name stays `Molpack`.** "Packer" was considered but rejected:
it would suggest molpack contains multiple packer implementations, and
the single entry-point type is the right shape for the crate. The
`molpack::Molpack` path stutter is accepted as the cost of an
unambiguous, domain-specific type name.

Every method below follows P1 (`with_*` iff optional with a default):
all current builder methods have defaults, so all become `with_*`;
`max_loops` has no default and stays as a positional argument to the
terminal `pack()`.

| Current | Proposed | Kind | Why |
|---|---|---|---|
| `struct Molpack` | *(unchanged)* | — | |
| `Molpack::new()` | *(unchanged)* | — | |
| `.add_handler(h)` | `.with_handler(h)` | R | P1 |
| `.add_restraint(r)` | `.with_global_restraint(r)` | R | disambiguate from target-scope |
| `.precision(f)` | `.with_precision(f)` | R | P1 |
| `.tolerance(f)` | `.with_tolerance(f)` | R | P1 |
| `.maxit(n)` | `.with_inner_iterations(n)` | R | P2 |
| `.nloop0(n)` | `.with_init_passes(n)` | R | P2 |
| `.sidemax(f)` | `.with_init_box_half_size(f)` | R | P2 |
| `.movefrac(f)` | `.with_perturb_fraction(f)` | R | P2 — "movebad" is implementation |
| `.movebadrandom(b)` | `.with_random_perturb(b)` | R | P2 |
| `.disable_movebad(b)` | `.with_perturb(b)` (`false` disables) | R | P3 double-neg |
| `.parallel_pair_eval(b)` | `.with_parallel_eval(b)` | R | P1 |
| `.pbc(min, max)` | removed — periodicity lives on `InsideBoxRestraint` (§5.1.1) | D | §2 relaxed for per-axis PBC |
| `.pbc_box(lengths)` | removed — see above | D | |
| *(implicit seed=0 in `.pack`)* | `.with_seed(u64)` — default `0` | N+S | split seed off `.pack()` per P1 |
| `.pack(&targets, max_loops, seed)` | `.pack(&targets, max_loops)` | S | `seed` → builder; `max_loops` stays positional |

#### 5.1.1 PBC lives on `InsideBoxRestraint`, not on `Molpack`

**Decision: Path Y + A.** The `Molpack` builder gains no PBC method.
Instead, [`InsideBoxRestraint`] (§5.7) takes a `periodic: [bool; 3]`
argument in its constructor — a box that confines atoms softly and can
additionally declare "wrap axis `k` in the pair-kernel minimum image."
The packer derives the system periodic box by scanning every
restraint on every target for a `Restraint::periodic_box()` override
(new trait method); at most one such declaration is allowed per run.

```rust
// Non-periodic confinement (the 90% case)
let r = InsideBoxRestraint::new([0.;3], [40.;3], [false; 3]);

// Fully-periodic cell (packmol-style PBC)
let r = InsideBoxRestraint::new([0.;3], [40.;3], [true; 3]);

// Slab geometry — XY periodic, Z confined
let r = InsideBoxRestraint::new([0.;3], [40.;3], [true, true, false]);

Target::new(frame, 100).with_restraint(r);
// ...no PBC method on Molpack.
```

**Why:**

- Fewer APIs — one fewer method on `Molpack`; PBC declaration is
  colocated with the only object that has the box dimensions.
- User can't declare PBC with a box shape that disagrees with any
  restraint (the shape *is* the restraint).
- Per-axis periodicity is real physics (slab geometries, 2D-periodic
  membranes). This lifts spec §2's "no new features" for this one
  item — documented here.

**Mechanics:**

- `Restraint` trait gains
  `fn periodic_box(&self) -> Option<([F;3], [F;3], [bool;3])>` with
  default `None`.
- `InsideBoxRestraint::periodic_box()` returns `Some((min, max, periodic))`
  iff any axis is periodic.
- `Molpack::pack()` scans all restraints; duplicates with identical
  bounds+flags are allowed (broadcast via `with_global_restraint` and
  re-use across targets); mismatched declarations error as
  [`PackError::ConflictingPeriodicBoxes`].
- `PackContext` gains `pbc_periodic: [bool; 3]`; `cell::axis_cell` /
  `cell::setcell` / `cell::delta_vector` take per-axis flags; the
  pair-kernel minimum image only wraps axes marked periodic; the cell
  list clamps (instead of wrapping) on non-periodic axes.

**Cost accepted:**

- A user who wants PBC without an inside-box confinement has no
  expression path. This is not a supported use case: in practice
  periodic packing always pairs with a box-shaped confinement to keep
  atoms in the reference image.

### 5.2 `Target`

| Current | Proposed | Kind | Why |
|---|---|---|---|
| `Target::new(frame, count)` | `Target::new(frame, count)` | — | keep |
| `Target::from_coords(pos, radii, count)` | keep | — | |
| `.with_name(s)` | keep | — | |
| `.with_restraint(r)` | keep | — | molecule scope |
| `.with_restraint_for_atoms(&[i..], r)` | `.with_atom_restraint(&[i..], r)` | R+S | drop glue word; indices become 0-based |
| `.with_relaxer(r)` | keep | — | |
| `.with_maxmove(n)` | `.with_perturb_budget(n)` | R | P2 |
| `.with_center()` | `.with_centering(CenteringMode::Center)` | M | P4 |
| `.without_centering()` | `.with_centering(CenteringMode::Off)` | M | P3 double-neg + P4 |
| `.constrain_rotation_x(c_deg, w_deg)` | `.with_rotation_bound(Axis::X, center, half_width)` | M+S | P4; `center`/`half_width` are `Angle` |
| `.constrain_rotation_y(c_deg, w_deg)` | `.with_rotation_bound(Axis::Y, center, half_width)` | M+S | P4 |
| `.constrain_rotation_z(c_deg, w_deg)` | `.with_rotation_bound(Axis::Z, center, half_width)` | M+S | P4 |
| `.fixed_at(pos)` | `.fixed_at(pos)` | — | keep — but asserts `count == 1` loudly |
| `.fixed_at_with_euler(pos, eul_rad)` | `.fixed_at(pos).with_orientation([Angle; 3])` | M+S | P4; Euler entries are `Angle` |
| `struct FixedPlacement { position, euler }` | `struct Placement { position, orientation: [Angle; 3] }` | R | drop redundant `Fixed`; angle type explicit |

The 1-based → 0-based index flip on `with_atom_restraint` is the one
behavioral change in the rename set. The `docs/packmol_parity.md` page
already notes 0-based is Rust-native; Packmol input parsers live in a
separate layer anyway.

New types introduced by this stage:

```rust
pub enum Axis { X, Y, Z }

/// Angular quantity — constructors enforce the unit at the call site.
/// Internally stored as radians; implements `Copy`.
#[derive(Debug, Clone, Copy)]
pub struct Angle(F);

impl Angle {
    pub const fn from_radians(rad: F) -> Self { Self(rad) }
    pub fn from_degrees(deg: F) -> Self { Self(deg.to_radians()) }
    pub fn radians(self) -> F { self.0 }
    pub fn degrees(self) -> F { self.0.to_degrees() }
}
```

`Angle` replaces all bare-`F` rotation arguments in the public surface.
Examples update from `.constrain_rotation_x(30.0, 10.0)` to
`.with_rotation_bound(Axis::X, Angle::from_degrees(30.0), Angle::from_degrees(10.0))`
— verbose, but unambiguous about units and impossible to misuse.

### 5.3 `CenteringMode`

| Current | Proposed | Kind | Why |
|---|---|---|---|
| `CenteringMode::Auto` | keep | — | |
| `CenteringMode::Center` | keep | — | |
| `CenteringMode::None` | `CenteringMode::Off` | R | avoid collision with `Option::None` |

### 5.4 `Handler` + info structs

| Current | Proposed | Kind | Why |
|---|---|---|---|
| `Handler::on_start(n_at, n_mol)` | keep | — | |
| `Handler::on_initial(sys)` | `on_initialized(sys)` | R | past-tense consistency |
| `Handler::on_step(info, sys)` | keep | — | |
| `Handler::on_phase_start(info)` | keep | — | |
| `Handler::on_phase_end(info, report)` | keep | — | |
| `Handler::on_inner_iter(i, f, sys)` | keep | — | |
| `Handler::on_finish(sys)` | keep | — | |
| `Handler::should_stop()` | keep | — | |
| `StepInfo.improvement_pct` | keep | — | `_pct` carries the unit |
| `StepInfo.hook_acceptance` | `relaxer_acceptance` | R | drop Hook legacy (§5.5) |
| `EarlyStopHandler::with_warmup(n)` | keep | — | already P1 |
| `EarlyStopHandler::with_patience(n)` | keep | — | |
| `XYZHandler::new(path, every)` | keep | — | |

### 5.5 `Relaxer` + `TorsionMcRelaxer`

Drop the `Hook` legacy aliases entirely — no deprecation, since they
predate 1.0 and are already duplicates.

| Current | Proposed | Kind | Why |
|---|---|---|---|
| `pub use Relaxer as Hook` | deleted | D | legacy dup |
| `pub use RelaxerRunner as HookRunner` | deleted | D | |
| `pub use TorsionMcRelaxer as TorsionMcHook` | deleted | D | |
| `Relaxer::build(ref_coords)` | `Relaxer::spawn(ref_coords)` | R | "build" collides with common Builder idiom; "spawn" reads as "create a runner" |
| `pub fn compute_excluded_pairs(g)` | move to `TorsionMcRelaxer::excluded_pairs(g)` | R+S | namespace the helper |
| `pub fn self_avoidance_penalty(...)` | `TorsionMcRelaxer::self_avoidance_penalty(...)` | R+S | |
| `TorsionMcRelaxer::with_temperature(t)` | keep | — | |
| `TorsionMcRelaxer::with_steps(n)` | keep | — | |
| `TorsionMcRelaxer::with_max_delta(r)` | keep | — | |
| `TorsionMcRelaxer::with_self_avoidance(r)` | keep | — | |

### 5.6 `Region` + combinators

| Current | Proposed | Kind | Why |
|---|---|---|---|
| `trait Region` | keep | — | |
| `trait RegionExt` | keep | — | |
| `struct And<A, B>` | keep | — | |
| `struct Or<A, B>` | keep | — | |
| `struct Not<A>` | keep | — | |
| `struct FromRegion<R>(pub R)` | `struct RegionRestraint<R>(pub R)` | R | "From" is meaningless as a type name |
| `struct BBox` | `struct Aabb` | R | axis-aligned bounding box, standard abbrev |
| `struct InsideBoxRegion` | keep | — | |
| `struct InsideSphereRegion` | keep | — | |
| `struct OutsideSphereRegion` | keep | — | |

Also adds a convenience on `RegionExt`:

```rust
impl<R: Region + 'static> RegionExt for R {
    fn into_restraint(self) -> RegionRestraint<Self> { RegionRestraint(self) }
}
```

### 5.7 `Restraint` + 14 concrete impls

All 14 `*Restraint` names **stay as-is**. The suffix is the key
disambiguation signal against `*Region` counterparts; stripping it
would make `InsideSphere` ambiguous at the call site. Confirmed by
the examples: every current usage reads clearly.

The only touches in this bucket:

- `InsideBoxRestraint::new(min, max)` → `new(min, max, periodic: [bool; 3])`
  (see §5.1.1 — declares PBC when any axis flag is true).
- `InsideBoxRestraint::cube_from_origin(origin, side)` → same, plus
  `periodic: [bool; 3]`.
- `InsideBoxRestraint::from_simbox(sim)` → `from_simbox(sim, periodic: [bool; 3])`.
- `Restraint` trait gains
  `fn periodic_box(&self) -> Option<([F;3], [F;3], [bool;3])>` with
  default `None`. Only `InsideBoxRestraint` overrides it in the built-in set.

### 5.8 `molpack::prelude`

A new `prelude` module is added, containing the items a typical
packing script needs. The crate root keeps its current re-exports (for
direct access and rustdoc discoverability), but users can write:

```rust
use molpack::prelude::*;

let target = Target::new(frame, 100).with_restraint(sphere);
let result = Molpack::new().with_handler(progress).pack(&[target], 200)?;
```

Initial contents (exact list refined in Stage 1):

```rust
pub use crate::{
    Molpack, PackResult, PackError,
    Target, CenteringMode, Placement, Axis, Angle,
    Handler, NullHandler, ProgressHandler, EarlyStopHandler, XYZHandler,
    StepInfo, PhaseInfo, PhaseReport,
    Restraint,
    InsideBoxRestraint, InsideCubeRestraint, InsideSphereRestraint,
    InsideEllipsoidRestraint, InsideCylinderRestraint,
    OutsideBoxRestraint, OutsideCubeRestraint, OutsideSphereRestraint,
    OutsideEllipsoidRestraint, OutsideCylinderRestraint,
    AbovePlaneRestraint, BelowPlaneRestraint,
    AboveGaussianRestraint, BelowGaussianRestraint,
    Region, RegionExt, RegionRestraint, Aabb,
    InsideBoxRegion, InsideSphereRegion, OutsideSphereRegion,
    And, Or, Not,
    Relaxer, RelaxerRunner, TorsionMcRelaxer,
};
```

The 14 concrete restraints are all re-exported — they are the
primary user-facing vocabulary, and a prelude that makes users import
them one by one defeats the purpose.

## 6. Backward-compatibility policy

Two tiers:

### Tier A — keep a deprecated alias for one release (6 months)

Applied to **structure-level** renames that users reference by name:

- `FromRegion` → `RegionRestraint`: deprecated type alias.
- `FixedPlacement` → `Placement`: deprecated type alias.
- `BBox` → `Aabb`: deprecated type alias.

(No alias needed for `Molpack` — the type name stays.)

### Tier B — hard break

Applied to **method-level** renames on builders. Builder chains break at
compile time with a clear error, and the fix is a one-line sed. No
deprecated method stubs — they would fork the surface and defeat the
purpose of this cleanup.

Also hard-break:

- `Hook` / `HookRunner` / `TorsionMcHook` re-exports (never were
  canonical).
- `pub fn compute_excluded_pairs` / `pub fn self_avoidance_penalty` at
  crate root (move to impl).

This means the rename lands as **molpack 0.N+1** with a migration guide
(Stage 8).

## 7. Migration checklist

Ordered so each stage leaves the tree green (`cargo test` only — Python
binding is out of scope and will be addressed in a follow-up).

1. Spec approved (this document).
2. Audit rename call sites — produce line-count table.
3. **Stage 1** `Molpack` builder renames + `Angle` / `Periodic` new types + `prelude` module.
4. **Stage 2** `Target` builder renames + `CenteringMode` / `Placement` / `Axis`.
5. **Stage 3** Handler / Relaxer cleanup, drop `Hook` aliases.
6. **Stage 4** `Region` / `FromRegion` → `RegionRestraint` / `BBox` → `Aabb`.
7. **Stage 5** Rust examples + benches + tests rewritten.
8. **Stage 6** `docs/*.md` + rustdoc + lib.rs "surface at a glance" table refreshed.
9. **Stage 7** CHANGELOG entry + migration sed snippet.
10. **(Follow-up)** Python binding updated in a separate spec.

Each stage lands as its own commit or PR; all of them go through the
existing five-example parity check (mixture / interface / spherical /
bilayer / solvprotein) before merge. Python tests are allowed to
temporarily fail between the Rust merge and the follow-up Python pass;
the `python/` directory is left alone during Stages 1–7.

## 8. Open questions

Decided during review (no longer open):

- **Type name.** `Molpack` stays; not renamed to `Packer`.
- **`with_*` semantics.** `with_*` implies the argument is optional with
  a default; required args live in `new()` or positional terminal
  arguments. Applied consistently across §5.
- **Python scope.** Out of scope for this spec; handled in a follow-up.
- **`prelude` module.** Added (see §5.8).
- **`Angle` newtype.** Adopted; all rotation-bound arguments now take
  `Angle` (see §5.2).

Decided in Stage 1:

- **PBC shape.** Moved onto `InsideBoxRestraint` with a `periodic: [bool; 3]`
  argument (§5.1.1). §2's "no new features" is relaxed for per-axis
  periodicity, because it is load-bearing for slab geometries (2D
  membranes, surface adsorption) that are realistic workloads.

Still open — decide during Stage 3:

1. **Relaxer runner constructor name.** `Relaxer::build` vs
   `Relaxer::spawn` vs `Relaxer::into_runner`. "Build" collides with
   the generic builder idiom; "spawn" reads as async; "into_runner" is
   the most Rust-idiomatic conversion name.
