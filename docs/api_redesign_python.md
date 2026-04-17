# molpack Python API redesign

Companion to `docs/api_redesign.md` (Rust). The Rust redesign (commit `0c0a896`)
left the Python bindings in a broken state — 10 `cargo check` errors from
renamed/removed Rust methods. This spec defines the Python-side redesign,
aligns it with the Rust conventions where sensible, and calls out the places
where Python should diverge for ergonomics.

## 1. Principles

- **P1 — 1:1 Rust parity on type names and semantics.** Restraint classes
  carry the full `Restraint` suffix (`InsideBoxRestraint`, not `InsideBox`);
  builders mirror the Rust verbs (`with_perturb_budget`, etc.). Users
  reading one language can predict the other.
- **P2 — Frame is a first-class citizen.** `Target` accepts a Frame;
  `PackResult` exposes a Frame. molpack does not reimplement I/O — no
  `to_pdb`, `to_xyz`, `to_frame`. Users compose with `molrs` (or any
  Frame-compatible tool) to persist or transform results.
- **P3 — Python enums and Angle class, not raw strings/floats.** Units
  and states are explicit at the call site: `Angle.from_degrees(30)`,
  `Axis.Z`, `CenteringMode.AUTO`.
- **P4 — Hard break, no deprecated aliases.** Python 0.1.0 has no installed
  users beyond in-tree tests; migration cost is bounded by the test + example
  sweep. The CHANGELOG lists the renames, no shim code in the binding layer.
- **P5 — Indexing matches Rust.** `with_atom_restraint` takes **0-based**
  indices (was 1-based on Python side). Callers porting from Packmol `.inp`
  subtract 1 at the call site, just like the Rust side.

## 2. Scope

**In scope:**
- `python/src/*.rs` — pyclass wrappers
- `python/python/molpack/__init__.py` — re-exports
- `python/python/molpack/molpack.pyi` — type stubs
- `python/tests/*.py`, `python/examples/*.py` — all 14 test/example files
- `python/docs/docs/**/*.md` — 9 markdown files

**Out of scope (this spec):**
- New Python-only features not in Rust (e.g. `to_pdb()`, context-manager
  packing). Those go in a separate follow-up.
- `molcrafts-molpack` pypi publish — the build must compile green; publish is
  a subsequent deploy step.

## 3. Design choices (decide before implementing)

### 3.1 Restraint class names

**Decision (P1): 1:1 with Rust.** `InsideBoxRestraint`, `InsideSphereRestraint`,
`OutsideSphereRestraint`, `AbovePlaneRestraint`, `BelowPlaneRestraint`.

Rationale: predictability across languages beats terseness. Users reading
the Rust source find the same name in Python and vice versa. The 181
call-site rename is bounded by test/example files which we're already
rewriting.

### 3.1.1 Restraint constructor parameter order

Match Rust positional order exactly:

| Class | Signature |
|-------|-----------|
| `InsideBoxRestraint` | `(min, max, periodic=(False, False, False))` |
| `InsideSphereRestraint` | `(center, radius)` — flipped from current `(radius, center)` |
| `OutsideSphereRestraint` | `(center, radius)` — flipped |
| `AbovePlaneRestraint` | `(normal, distance)` |
| `BelowPlaneRestraint` | `(normal, distance)` |

### 3.2 Angle representation

| Option | Example | Trade-off |
|--------|---------|-----------|
| **A (recommend)** | `Angle.from_degrees(30)` / `Angle.from_radians(pi/6)` | Mirrors Rust `Angle`; units explicit at call site; type-checker catches mix-ups. |
| B | Raw float + suffix convention (`_deg`/`_rad`) | Regresses to the ambiguity the Rust redesign explicitly removed. |
| C | Module helpers `deg(30)` / `rad(pi/6)` returning a typed value | Short, but unusual — diverges from Rust naming without meaningful gain. |

**Recommendation: A.** Expose `molpack.Angle` as a Python class with
`from_degrees` / `from_radians` classmethods and `.degrees` / `.radians`
getters. Identical to Rust surface.

### 3.3 Axis and CenteringMode

| Option | Example | Trade-off |
|--------|---------|-----------|
| **A (recommend)** | `molpack.Axis.X` / `molpack.CenteringMode.AUTO` (Python `enum.Enum`) | Real enum — IDE introspection, `repr`, equality. |
| B | `Literal["x", "y", "z"]` strings | Minimalist, no import, but typos surface only at runtime. |

**Recommendation: A.** Enum subclass exported from the binding. No string
overload (keeps the surface tight).

### 3.4 `InsideBoxRestraint` periodic parameter

The Rust signature is `InsideBoxRestraint::new(min, max, [bool; 3])` with
`periodic` **required**. Python options:

| Option | Signature | Notes |
|--------|-----------|-------|
| **A (recommend)** | `InsideBoxRestraint(min, max, periodic=(False, False, False))` | Keyword-only default keeps the 90% non-periodic case clean; still accepts a tuple/list of 3 bools. |
| B | `InsideBoxRestraint(min, max, periodic)` required positional | Matches Rust exactly, but forces every non-periodic call site to type `(False, False, False)`. |
| C | Two constructors — `InsideBoxRestraint(...)` and `InsideBoxRestraint.periodic(...)` | Extra surface area; ugly when only one axis is periodic. |

**Recommendation: A.** Optional keyword, default to non-periodic, still
accepts per-axis booleans.

### 3.5 Molpack seed placement

Rust moved `seed` from `.pack(...)` into `.with_seed(n)`. Python should
match — consistency beats familiarity, and the only Python-side churn is
inside tests/examples (already being rewritten).

New signature:
```python
Molpack().with_seed(42).pack(targets, max_loops=200)
# instead of
Molpack().pack(targets, max_loops=200, seed=42)
```

## 4. Rename table

### 4.1 `Target` builder

| Old (Python) | New | Notes |
|--------------|-----|-------|
| `Target(name, frame, count)` | `Target(frame, count)` | `name` becomes optional via `.with_name()` |
| `.with_restraint_for_atoms(indices_1b, r)` | `.with_atom_restraint(indices_0b, r)` | Semantic change: flip to 0-based |
| `.with_maxmove(n)` | `.with_perturb_budget(n)` | — |
| `.with_center()` | `.with_centering(CenteringMode.CENTER)` | — |
| `.without_centering()` | `.with_centering(CenteringMode.OFF)` | — |
| `.constrain_rotation_x(c_deg, w_deg)` | `.with_rotation_bound(Axis.X, Angle.from_degrees(c), Angle.from_degrees(w))` | — |
| `.constrain_rotation_y(...)` | `.with_rotation_bound(Axis.Y, ...)` | — |
| `.constrain_rotation_z(...)` | `.with_rotation_bound(Axis.Z, ...)` | — |
| `.fixed_at_with_euler(pos, euler_rad)` | `.fixed_at(pos).with_orientation((a, b, c))` | Each `a/b/c` is an `Angle` |

Unchanged: `.with_name`, `.with_restraint`, `.fixed_at`, `.natoms`, `.count`,
`.elements`, `.is_fixed`, `__repr__`.

**Not added** (explicitly out of scope per P2 "Frame is a first-class
citizen"): `Target.from_coords(coords, radii, count)`. Users always go
through a Frame — if they have raw numpy arrays, they build a
Frame-compatible dict. This keeps one input path, one invariant.

### 4.2 `Molpack` builder

| Old | New | Notes |
|-----|-----|-------|
| `Molpack(tolerance=2.0, precision=0.01, maxit=20, nloop0=0, sidemax=1000.0, movefrac=0.05, movebadrandom=False, disable_movebad=False, progress=True)` | `Molpack()` with everything as `with_*` | `__init__` becomes zero-arg; all tuning goes through builders. |
| `.with_maxit(n)` | `.with_inner_iterations(n)` | Rust renamed this in Stage 1 |
| `.with_nloop0(n)` | `.with_warmup_iterations(n)` | Rust renamed this |
| `.with_sidemax(x)` | `.with_max_box_side(x)` | Rust renamed this |
| `.with_movefrac(x)` | `.with_perturb_fraction(x)` | Rust renamed this |
| `.with_movebadrandom(b)` | `.with_random_perturb(b)` | Rust renamed this |
| `.with_disable_movebad(b)` | `.without_perturb_phase(b)` | Rust renamed this |
| `.with_pbc(min, max)` | **DELETE** | PBC derives from `InsideBox(..., periodic=...)` |
| `.with_pbc_box(lengths)` | **DELETE** | Same |
| `.pack(targets, max_loops, seed=None)` | `.with_seed(n).pack(targets, max_loops=200)` | Seed moved to builder |
| — | `.with_global_restraint(r)` | **New** — mirror Rust |

Unchanged: `.with_tolerance`, `.with_precision`, `.with_progress`,
`.add_handler`, `pack()` base signature.

### 4.3 Restraint constructors

| Old | New |
|-----|-----|
| `InsideBox([0,0,0], [30,30,30])` | `InsideBoxRestraint([0,0,0], [30,30,30], periodic=(False, False, False))` |
| `InsideSphere(5.0, [0,0,0])` | `InsideSphereRestraint([0,0,0], 5.0)` — param order flipped |
| `OutsideSphere(5.0, [0,0,0])` | `OutsideSphereRestraint([0,0,0], 5.0)` — param order flipped |
| `AbovePlane([0,0,1], 10.0)` | `AbovePlaneRestraint([0,0,1], 10.0)` |
| `BelowPlane([0,0,1], 10.0)` | `BelowPlaneRestraint([0,0,1], 10.0)` |

### 4.4 New Python-exposed types

| Type | Source | Surface |
|------|--------|---------|
| `Angle` | `molpack.Angle` | `from_degrees(f) -> Angle`, `from_radians(f) -> Angle`, `.degrees`, `.radians`, `.ZERO` |
| `Axis` | `molpack.Axis` | `Enum` with `X`, `Y`, `Z` |
| `CenteringMode` | `molpack.CenteringMode` | `Enum` with `AUTO`, `CENTER`, `OFF` |

Not exposed (internal-only): `Placement`, `Aabb`, `RegionRestraint`,
`Relaxer`, `Handler` — users interact via duck-typed objects where relevant.

## 5. Stage plan

| Stage | Scope | Files |
|-------|-------|-------|
| 1 | Add `Angle`, `Axis`, `CenteringMode` pyclasses | `python/src/target.rs`, `python/src/lib.rs` |
| 2 | Rename restraint pyclasses to `*Restraint`; flip sphere param order; `InsideBoxRestraint` accepts `periodic` kwarg | `python/src/constraint.rs` |
| 3 | `Target` builder renames (includes 1→0 base index flip) + name optional | `python/src/target.rs` |
| 4 | `Molpack` zero-arg `__init__` + builder renames + drop PBC methods + `with_seed` | `python/src/packer.rs` |
| 5 | `Molpack::with_global_restraint`; typed `PackError` subclasses | `python/src/packer.rs`, `python/src/helpers.rs` |
| 6 | `PackResult.frame` property (molrs.Frame), cached `positions`; `StepInfo.relaxer_acceptance` | `python/src/packer.rs`, `python/src/handler.rs` |
| 7 | Rewrite `__init__.py` + `molpack.pyi` (add Protocols for Handler/Restraint, all new types, `Self` returns) | 2 files |
| 8 | Update all tests (`python/tests/*.py`) | 8 files |
| 9 | Update all examples (`python/examples/*.py`) | 6 files |
| 10 | Update `python/docs/docs/**/*.md` | 9 files |
| 11 | Verify: `cargo check`, `maturin develop`, `pytest`, `ruff`, `ty` | — |

## 6. Migration examples

**Before:**
```python
import molpack

target = (
    molpack.Target("water", frame, 50)
    .with_restraint(molpack.InsideBox([0,0,0], [30,30,30]))
    .with_restraint_for_atoms([1, 2], molpack.AbovePlane([0,0,1], 10.0))  # 1-based!
    .constrain_rotation_z(0.0, 30.0)
    .with_maxmove(20)
)

packer = molpack.Molpack(tolerance=2.0, progress=False).with_pbc_box([30,30,30])
result = packer.pack([target], max_loops=200, seed=42)
print(result.positions)  # user had to rebuild a frame by hand
```

**After:**
```python
import molpack
from molpack import Angle, Axis, CenteringMode

target = (
    molpack.Target(frame, 50)
    .with_name("water")
    .with_restraint(molpack.InsideBoxRestraint(
        [0,0,0], [30,30,30], periodic=(True, True, True),
    ))
    .with_atom_restraint([0, 1], molpack.AbovePlaneRestraint([0,0,1], 10.0))  # 0-based
    .with_rotation_bound(Axis.Z, Angle.from_degrees(0.0), Angle.from_degrees(30.0))
    .with_perturb_budget(20)
)

packer = (
    molpack.Molpack()
    .with_tolerance(2.0)
    .with_progress(False)
    .with_seed(42)
)
result = packer.pack([target], max_loops=200)

# Frame-first output — user picks the writer.
import molrs
molrs.write_pdb("out.pdb", result.frame)
```

## 7. Frame as the primary I/O boundary

`PackResult` exposes the packed result as a single **Frame** object
(`molrs.Frame`-compatible). molpack does **not** provide writers — users
compose with `molrs` / `molpy` / their own tooling.

Exposed on `PackResult`:

| Property | Type | Notes |
|----------|------|-------|
| `frame` | `molrs.Frame` | Canonical output. Contains `atoms` block with `x`, `y`, `z`, `element` columns at minimum. |
| `positions` | `np.ndarray[(N,3), float64]` | Convenience view; cached. |
| `natoms`, `converged`, `fdist`, `frest` | `int` / `bool` / `float` | Scalars. |

Explicitly **not** on `PackResult`:
- `to_pdb(path)`, `to_xyz(path)`, `to_frame()`, any writer. Out of scope.

## 8. Deferred / recorded for v2

- Python-idiomatic shortcuts if users request them after landing — e.g.
  `with_rotation_bound_degrees("z", 0, 30)` convenience. Not in v1; revisit
  based on real feedback.
- `TorsionMcRelaxer` Python class — Rust-side relaxer system is still
  maturing; expose in a follow-up once API stabilizes.
