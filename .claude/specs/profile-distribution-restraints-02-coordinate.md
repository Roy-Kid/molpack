---
title: Profile-restraint reaction coordinate ξ(x) with analytic gradient
status: done
created: 2026-06-12
---

# Profile-restraint reaction coordinate ξ(x) with analytic gradient

## Summary
Add a leaf geometry module `src/restraint/profile/coordinate.rs` that defines the
reaction coordinate ξ(x) shared by the whole profile-restraint family. A
`Coordinate` maps a Cartesian position x to a scalar ξ (a length in Å) and supplies
the analytic gradient ∇ξ, so that a later composing restraint can apply the chain
rule ∇_xU = (dU/dξ)·∇ξ. Four coordinate flavours are provided — planar, radial,
cylindrical, and region-distance — each minimum-image aware and each guarded so the
gradient stays finite even at its geometric singularity. This sub-spec is pure
geometry: no penalty/distribution math, no `Restraint` impl, no script parsing, and
no Python — those land in later sub-specs of the chain.

> **Chain bootstrap:** this is the first sub-spec to add a file under
> `src/restraint/profile/`, so it also creates `src/restraint/profile/mod.rs` and
> adds `pub mod profile;` to `src/restraint/mod.rs` (created by `-01-split`). Later
> siblings (`-03`, `-04`, `-05`) only append `pub mod <name>;` to the existing
> `profile/mod.rs`.

## Domain basis
Four reaction coordinates, each [length, Å], with analytic gradient:
- planar:  ξ = n̂·(x−x₀), ‖n̂‖=1 ⇒ ∇ξ = n̂ (constant; NO singularity).
- radial:  ξ = ‖x−c‖ ⇒ ∇ξ = (x−c)/‖x−c‖ = r̂ (undefined at centre x=c; clamp r ≥ r_guard).
- cylindrical: ρ⃗ = (x−a) − ((x−a)·û)û, ‖û‖=1, ξ = ‖ρ⃗‖ ⇒ ∇ξ = ρ⃗/‖ρ⃗‖ (⟂ axis; undefined on axis; clamp s ≥ r_guard). Derivation: projection (I−ûûᵀ) is symmetric idempotent ⇒ ∇ξ = ρ⃗/‖ρ⃗‖.
- region-distance: ξ = signed distance d(x) to a region surface ⇒ ∇ξ = ∇d = outward unit normal (‖∇d‖=1 a.e.; kinks at medial axis / sharp region features — gencan sees a finite kink there, acceptable).

Inner clamp (radial & cylindrical only): for ξ < r_guard use ξ_eff = r_guard and
freeze ∇ξ to a finite unit vector (or zero the radial force inside r_guard); choose
r_guard ~ one bin width or the packing tolerance radius. Justification: r̂ is
undefined at the centre and an unguarded 1/ξ force diverges; r_guard keeps ∇ξ finite.

Periodicity (§6.5): the reference plane/centre/axis is FIXED in the box;
planar/radial/cylindrical deltas must be wrapped under minimum image per active
periodic axis (slabs straddle the boundary). Region-distance uses the region's own
convention.

References: Packmol soft-objective framework (Martínez et al. 2009,
DOI 10.1002/jcc.21224); coordinate gradients are elementary calculus (no external
cite needed).

## Design
A `Coordinate` enum (immutable, `Debug`, `Clone` where the payload allows, `Send +
Sync`) with variants:
- `Planar { normal: [F; 3], point: [F; 3] }` — `normal` normalized to unit length on
  construction.
- `Radial { center: [F; 3], r_guard: F }`.
- `Cylindrical { axis_origin: [F; 3], axis_dir: [F; 3], r_guard: F }` — `axis_dir`
  normalized on construction; degenerate (near-zero-length) axis rejected at
  construction.
- `RegionDistance(Arc<dyn Region>)` — defers ξ to `Region::signed_distance` and ∇ξ to
  `Region::signed_distance_grad`.

Two methods form the public surface:
- `fn xi(&self, x: &[F; 3], pbc: &PbcWrap) -> F`
- `fn grad_xi(&self, x: &[F; 3], pbc: &PbcWrap) -> [F; 3]`

`PbcWrap` is the minimum-image carrier passed in by the caller; this module does NOT
reinvent wrapping. The planar/radial/cylindrical deltas (`x − point`, `x − center`,
`x − axis_origin`) are wrapped via the existing `delta_vector` minimum-image helper
(`src/cell.rs:108`) before any projection/norm; `RegionDistance` ignores `pbc` and
relies on the region's own convention. Constructors are fallible where a degenerate
parameter is possible (`Cylindrical::new`, `Planar::new` reject a zero-length
direction/normal) and return `Result`; `Radial::new` and `RegionDistance::new` are
infallible. Lifecycle: a `Coordinate` is constructed once from validated parameters,
holds no mutable state, and is shared read-only across packing iterations (`Send +
Sync`). The inner clamp lives entirely inside `Radial`/`Cylindrical` `grad_xi`: when
the wrapped radius/cylindrical-distance `s < r_guard`, ξ is reported as `r_guard` and
∇ξ is set to a finite unit vector (a fixed fallback direction when exactly on the
singularity) so the returned gradient norm is bounded and never NaN/Inf.

## Files to create or modify
- `src/restraint/profile/coordinate.rs` (new)
- `src/restraint/profile/mod.rs` (new) — created here (first profile file); add `pub mod coordinate;` and re-export `Coordinate`
- `src/restraint/mod.rs` — add `pub mod profile;` (one-line edit; file created by `-01-split`)
- `tests/coordinate.rs` (new)

## Tasks
- [x] Write failing analytic-vs-FD gradient tests for planar/radial/cylindrical/region-distance coordinates (tests/coordinate.rs)
- [x] Write failing r_guard finite-gradient tests asserting bounded, non-NaN ∇ξ at the centre and on the axis (tests/coordinate.rs)
- [x] Write failing minimum-image wrap test comparing ξ against a hand-computed wrapped delta (tests/coordinate.rs)
- [x] Implement Coordinate enum, fallible constructors (normal/axis normalization, degenerate-axis rejection), xi, and grad_xi in src/restraint/profile/coordinate.rs
- [x] Create src/restraint/profile/mod.rs, wire `pub mod coordinate;`, re-export Coordinate, and add `pub mod profile;` to src/restraint/mod.rs
- [x] Add docstrings with units (Å) and the ∇ξ derivation note on each variant
- [x] Run full check + test suite

## Testing strategy
- Happy path — for each of planar, radial, cylindrical, region-distance: evaluate `xi`
  at a representative off-singularity point and confirm it matches the closed-form
  value; confirm `grad_xi` matches a 3-point central finite difference of `xi`
  componentwise within a stated tolerance (e.g. ≤ 1e-6 absolute).
- Edge cases — a query AT the centre (radial) and ON the axis (cylindrical) returns
  finite ξ and a finite, bounded-norm ∇ξ with no NaN/Inf; a query just inside r_guard
  reports ξ = r_guard; a degenerate (zero-length) axis/normal is rejected at
  construction.
- Periodicity — with an active periodic box, ξ for planar/radial/cylindrical uses the
  minimum-image delta: a reference point near a box face plus a query across the
  boundary yields the wrapped distance, verified against a hand-computed value; with
  no periodic axis active the result equals the raw distance.
- Domain validation — planar ∇ξ is the constant unit normal; radial/cylindrical ∇ξ has
  unit norm off-singularity; region-distance ∇ξ has unit norm a.e. (away from medial
  axis), matching ‖∇d‖ = 1.

## Out of scope
- Distribution/penalty/profile math (dU/dξ and the histogram penalty) — `-03`/`-04`.
- The composed `Restraint` impl that applies the chain rule — `-05-compose`.
- `.inp` script parsing / lowering for profile restraints.
- Python bindings for the coordinate or restraint.
