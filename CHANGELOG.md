# Changelog

All notable changes to `molcrafts-molpack` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed — API redesign (breaking)

Consistency pass per `docs/api_redesign.md`. Every method on builder
surfaces now uses `with_*` for optional knobs; non-optional arguments
are positional. Packmol-Fortran-jargon names are gone from the user
surface.

**`Molpack` builder**

| Old | New |
|---|---|
| `.add_handler(h)` | `.with_handler(h)` |
| `.add_restraint(r)` | `.with_global_restraint(r)` |
| `.precision(f)` | `.with_precision(f)` |
| `.tolerance(f)` | `.with_tolerance(f)` |
| `.maxit(n)` | `.with_inner_iterations(n)` |
| `.nloop0(n)` | `.with_init_passes(n)` |
| `.sidemax(f)` | `.with_init_box_half_size(f)` |
| `.movefrac(f)` | `.with_perturb_fraction(f)` |
| `.movebadrandom(b)` | `.with_random_perturb(b)` |
| `.disable_movebad(b)` | `.with_perturb(b)` *(semantics inverted — `false` disables)* |
| `.parallel_pair_eval(b)` | `.with_parallel_eval(b)` |
| `.pbc(min, max)` / `.pbc_box(lengths)` | declared via `InsideBoxRestraint::new(min, max, [bool; 3])` — see below |
| `.pack(&targets, max_loops, Option<u64>)` | `.with_seed(u64).pack(&targets, max_loops)` |

**`Target` builder**

| Old | New |
|---|---|
| `.with_restraint_for_atoms(&[1-based], r)` | `.with_atom_restraint(&[0-based], r)` — **indices flipped** |
| `.with_maxmove(n)` | `.with_perturb_budget(n)` |
| `.with_center()` | `.with_centering(CenteringMode::Center)` |
| `.without_centering()` | `.with_centering(CenteringMode::Off)` |
| `.constrain_rotation_x/y/z(c_deg, w_deg)` | `.with_rotation_bound(Axis::X/Y/Z, Angle, Angle)` |
| `.fixed_at_with_euler(pos, eul)` | `.fixed_at(pos).with_orientation([Angle; 3])` |
| `CenteringMode::None` | `CenteringMode::Off` (avoids collision with `Option::None`) |
| `FixedPlacement` | `Placement` (*deprecated alias kept*) |

**Per-axis periodic boundaries (new feature)**

`InsideBoxRestraint` now takes a required `periodic: [bool; 3]`. The
packer derives the system PBC by scanning all target restraints and
uses per-axis wrap in the pair kernel. Slab geometries become
expressible: `[true, true, false]` gives an X/Y-periodic cell with Z
confined.

```rust
// Non-periodic (Stage 1 migration: add `[false; 3]` to every call)
InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3])

// Fully-periodic (replaces old Molpack::pbc_box)
InsideBoxRestraint::new([0.0; 3], [40.0; 3], [true; 3])
```

Two declarations with mismatched bounds or flags return
`PackError::ConflictingPeriodicBoxes`.

**Types introduced**

- `Axis` — Cartesian axis selector (`X` / `Y` / `Z`).
- `Angle` — angular newtype with `from_degrees` / `from_radians`.
- `Placement` — `{ position, orientation: [Angle; 3] }` (was `FixedPlacement`).
- `RegionRestraint<R>` — renamed from `FromRegion`; deprecated alias kept.
- `Aabb` — renamed from `BBox`; deprecated alias kept.
- `molpack::prelude` — bulk re-export module (`use molpack::prelude::*;`).

**Types / trait methods renamed**

- `Restraint::periodic_box()` — new default trait method (returns `None`).
- `Handler::on_initial` → `on_initialized`.
- `StepInfo.hook_acceptance` → `relaxer_acceptance`.
- `Relaxer::build(ref_coords)` → `Relaxer::spawn(ref_coords)`.

**Re-exports removed**

- `molpack::Hook`, `molpack::HookRunner`, `molpack::TorsionMcHook`
  (use `Relaxer`, `RelaxerRunner`, `TorsionMcRelaxer`).
- `molpack::compute_excluded_pairs`, `molpack::self_avoidance_penalty`
  at crate root — import from `molpack::relaxer::*` instead.

**Error variants**

- `PackError::ConflictingPeriodicBoxes { first, second }` — new.
- `PackError::InvalidPBCBox` — unchanged, now raised when any
  declared periodic axis has non-positive extent.

### Migration snippet

Most call sites can be migrated mechanically:

```sh
# In your Rust project:
sed -i '' \
  -e 's/\.add_handler(/.with_handler(/g' \
  -e 's/\.add_restraint(/.with_global_restraint(/g' \
  -e 's/\.tolerance(/.with_tolerance(/g' \
  -e 's/\.precision(/.with_precision(/g' \
  -e 's/\.maxit(/.with_inner_iterations(/g' \
  -e 's/\.movefrac(/.with_perturb_fraction(/g' \
  -e 's/\.with_restraint_for_atoms(/.with_atom_restraint(/g' \
  -e 's/\.with_center()/.with_centering(CenteringMode::Center)/g' \
  -e 's/\.without_centering()/.with_centering(CenteringMode::Off)/g' \
  -e 's/CenteringMode::None/CenteringMode::Off/g' \
  -e 's/FromRegion(/RegionRestraint(/g' \
  src/**/*.rs
```

Three items must be fixed by hand:

1. `InsideBoxRestraint::new(min, max)` → add `[false; 3]` third argument (or `[true; 3]` for what was formerly `Molpack::pbc_box`).
2. `.with_atom_restraint(&[...])` indices are now 0-based — subtract 1 from each entry when porting from 1-based Packmol `.inp` inputs.
3. `Molpack::pack(&t, n, Some(s))` → `.with_seed(s).pack(&t, n)`.

### Added
- Initial implementation of the three-phase packing algorithm (per-type init → geometric pre-fit → main loop with `movebad` heuristic)
- GENCAN quasi-Newton optimizer with bound constraints (`gencan` module)
- 14 concrete `*Restraint` types covering Packmol kinds 2–15
- `Region` trait with `And` / `Or` / `Not` combinators and `RegionRestraint` lift
- `Handler` trait with four built-ins: `NullHandler`, `ProgressHandler`, `EarlyStopHandler`, `XYZHandler`
- `Relaxer` trait with `TorsionMcRelaxer` built-in
- Per-axis periodic boundary conditions (slab geometries) via `InsideBoxRestraint`
- Optional `rayon` feature for parallel objective evaluation
- Optional `cli` feature that gates the `molpack` binary and its `clap` dependency
- Python bindings via PyO3 + maturin (`molcrafts-molpack` on PyPI)
- Five canonical Packmol workloads as examples (mixture, bilayer, spherical, interface, solvprotein)
- Criterion benchmark suite (five benches covering end-to-end and hot-path microbenchmarks)
- Rustdoc chapters: `getting_started`, `concepts`, `architecture`, `extending`

[Unreleased]: https://github.com/MolCrafts/molpack/compare/HEAD...HEAD
