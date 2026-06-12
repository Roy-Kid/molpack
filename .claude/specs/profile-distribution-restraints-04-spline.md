---
title: Monotone cubic spline + Boltzmann inversion for tabulated profile restraints
status: approved
created: 2026-06-12
---

# Monotone cubic spline + Boltzmann inversion for tabulated profile restraints

## Summary
This sub-spec adds `TabulatedProfile`, the general C¹ distribution primitive underlying the profile-distribution restraint family. It interpolates a user-supplied target density ρ*(ξ) tabulated on a grid with a monotone (clamped / Fritsch–Carlson) cubic spline, exposes the analytic spline derivative s'(ξ), and performs numerical Boltzmann inversion to a smooth penalty U(ξ) and its force −dU/dξ. The three analytic families (Gaussian / erf / exponential) become special cases: a tabulated sampling of any of them reproduces the analytic U within interpolation tolerance. The helper handles the shell-volume Jacobian conversion (count histogram → volumetric density) and the density floor shared with `-03-distribution`, producing a penalty and gradient that are C¹ everywhere so gencan converges.

## Domain basis
Tabulated Boltzmann inversion:

- U(ξ) = −kT · ln( max(ρ*(ξ), ρ_min) / ρ₀ )
- dU/dξ = −kT · s'(ξ) / max(ρ*(ξ), ρ_min)

where s(ξ) is a C¹ interpolant (monotone / clamped cubic spline) of the tabulated ρ*, and s'(ξ) is its analytic derivative.

§6.4 Smoothness: the penalty AND its gradient must be C¹ for gencan to converge. A piecewise-linear target gives discontinuous forces and must NOT be used directly — interpolate with a smooth scheme. A monotone (Fritsch–Carlson / clamped) cubic avoids the overshoot a naive natural cubic would introduce on steep tabulated data.

§6.2 Jacobian: a tabulated input declared as a count histogram n(ξ) must be converted to ρ*(ξ) = n(ξ) / (dV/dξ) BEFORE inversion (planar dV/dξ = A const; radial = 4π ξ²; cylindrical = 2π ξ L). Default-safe input kind = histogram.

§6.3 Floor: replacing ρ* by max(ρ*, ρ_min) is equivalent to capping U at U_max = −kT · ln(ρ_min/ρ₀); zero / empty bins stay finite. ρ₀ shifts U by a constant only and has no force effect.

Special-case identity: the three analytic families (Gaussian / erf / exponential) are special cases of the tabulated one — a tabulated sampling of an analytic ρ* must reproduce the analytic U within interpolation tolerance. For a Gaussian ρ* ∝ exp(−(ξ−μ)²/2σ²), inversion yields the harmonic U(ξ) = (kT/2σ²)(ξ−μ)² + const.

References:
- [1] Boltzmann inversion, E = −kT ln ρ + const (foundational, classical statistical mechanics).
- [2] Iterative Boltzmann inversion (exact-match upgrade path, OUT OF SCOPE here, cited only): D. Reith, M. Pütz, F. Müller-Plathe, *Deriving effective mesoscale potentials from atomistic simulations*, J. Comput. Chem. 24:1624 (2003), DOI 10.1002/jcc.10307.

## Design
A single new module `src/restraint/profile/spline.rs` exposing one public type `TabulatedProfile`, immutable, `Debug`, `Send + Sync`.

Construction (`TabulatedProfile::new`) takes sorted grid nodes ξ_i, raw tabulated values, an `input_kind` (histogram vs. density), the geometry needed for the shell Jacobian, and the density floor ρ_min. At construction it:

1. Validates the grid: ≥ 2 nodes, strictly increasing ξ_i, all values finite and non-negative.
2. If `input_kind` is histogram, converts each value n(ξ_i) → ρ*(ξ_i) = n(ξ_i) / (dV/dξ)(ξ_i) using the shell Jacobian (planar / radial / cylindrical). Stores the resulting volumetric ρ* — so downstream evaluation always sees a density.
3. Fits the monotone clamped cubic: computes secant slopes, initial node derivatives, then applies the Fritsch–Carlson limiter so no segment overshoots. Stores per-segment Hermite coefficients (one cubic per `[ξ_i, ξ_{i+1}]`).

Evaluation methods (all take ξ; floor applied at evaluation, not stored into ρ*):

- `value(ξ) → ρ*` — spline value s(ξ); flat extrapolation (clamp ξ to end nodes) outside the grid, continuous in value and derivative.
- `deriv(ξ) → ρ*'` — analytic spline derivative s'(ξ); zero in the extrapolated flat regions.
- `u(ξ, kt) → F` — U(ξ) = −kT · ln( max(s(ξ), ρ_min) / ρ₀ ).
- `du_dxi(ξ, kt) → F` — dU/dξ = −kT · s'(ξ) / max(s(ξ), ρ_min); zero force where the floor is active (s(ξ) ≤ ρ_min).

The shell-Jacobian helper and `input_kind` / `density_floor` concepts are shared with `-03-distribution`. If `-03` already defines `ShellJacobian`, `InputKind`, and the floor type, reuse them verbatim (`use` from the sibling module); otherwise mirror minimal definitions in this module so 04 stays self-contained and testable on its own. No analytic restraint composition lives here — `TabulatedProfile` is a pure numerical primitive consumed by `-05`.

## Files to create or modify
- `src/restraint/profile/spline.rs` (new) — `TabulatedProfile`: validation, Jacobian conversion, monotone clamped cubic fit, analytic derivative, Boltzmann inversion `u` / `du_dxi`, flat extrapolation, floor at evaluation.
- `src/restraint/profile/mod.rs` — add `pub mod spline;` and re-export `TabulatedProfile` (module home created by `-02-coordinate`; this sub-spec only appends the declaration / re-export).

## Tasks
- [ ] Write failing tests for grid validation (`src/restraint/profile/spline.rs`): rejects < 2 nodes, non-increasing ξ, non-finite or negative values
- [ ] Write failing tests for node interpolation and C¹ continuity across interior segment boundaries (`src/restraint/profile/spline.rs`)
- [ ] Write failing tests comparing analytic `deriv` against central finite difference of `value` over the grid interior (`src/restraint/profile/spline.rs`)
- [ ] Write failing tests for monotone (overshoot-free) behavior on steep monotone input and flat extrapolation outside the grid (`src/restraint/profile/spline.rs`)
- [ ] Write failing tests for Boltzmann inversion: tabulated Gaussian reproduces harmonic U within tolerance, zero bin yields finite `u` / `du_dxi` via floor, radial histogram reproduces intended density only with the ∝ξ² Jacobian (`src/restraint/profile/spline.rs`)
- [ ] Implement `TabulatedProfile::new` with validation and histogram→density shell-Jacobian conversion in `src/restraint/profile/spline.rs`
- [ ] Implement the monotone clamped cubic fit (Fritsch–Carlson limiter) and per-segment coefficient storage in `src/restraint/profile/spline.rs`
- [ ] Implement `value`, `deriv` with flat extrapolation, and `u` / `du_dxi` with the density floor in `src/restraint/profile/spline.rs`
- [ ] Add module-level and method docstrings per project doc style with units (ξ length, ρ* number density, U energy, kt energy) and wire `pub mod spline;` + re-export in `src/restraint/profile/mod.rs`
- [ ] Run full check + test suite (`cargo fmt`, `cargo clippy -- -D warnings`, `cargo test -p molcrafts-molpack --lib --tests`)

## Testing strategy
Happy path:
- Spline passes through every tabulated node `(ξ_i, ρ*_i)` to within f64 round-off.
- `value` and `deriv` agree across each interior node from the left and right segments (C¹).

Edge cases:
- Grid validation rejects < 2 nodes, non-strictly-increasing ξ, NaN/∞ values, negative densities.
- Query ξ below the first node and above the last node returns the end-node value with zero derivative (flat extrapolation), continuous at the boundary.
- A zero-valued bin yields a finite `u` (capped at U_max) and finite `du_dxi`.

Domain validation:
- Analytic `deriv` matches a central finite difference of `value` within a stated tolerance across the grid interior.
- A monotone-decreasing tabulated input produces a monotone spline (no value exceeds the local node bracket — no overshoot).
- A tabulated sampling of a Gaussian density ρ* ∝ exp(−(ξ−μ)²/2σ²) reproduces the analytic harmonic U(ξ) = (kT/2σ²)(ξ−μ)² within a stated interpolation tolerance (up to the additive ρ₀ constant), confirming the analytic families are special cases.
- A radial count histogram reproduces the intended uniform density only when the 4π ξ² shell Jacobian is applied; omitting it produces a visibly ξ-dependent (wrong) density.

## Out of scope
- Coordinate geometry / ξ extraction (`-02`).
- The composed `Restraint` impl wiring `TabulatedProfile` into `f` / `fg` (`-05`).
- Iterative Boltzmann inversion (IBI) outer loop (ref [2], cited only).
- `.inp` script parsing / lowering and the Python wheel surface.
