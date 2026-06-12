---
title: ProfileRestraint — compose coordinate × distribution into a Restraint
status: done
created: 2026-06-12
---

# ProfileRestraint — compose coordinate × distribution into a Restraint

## Summary

Introduce `ProfileRestraint`, the single composed restraint that turns the profile family into something the objective can actually call. It pairs a reaction coordinate ξ(x) (from sub-spec -02) with a 1-D distribution penalty U(ξ) (analytic from -03 or tabulated from -04) and implements the existing `Restraint` trait. Once landed, a biased site assigned to this restraint via the per-atom CSR sees a penalty U(ξ(x)) pulling it toward the target profile, scaled by the global soft-start weight so it never fights the overlap term. This is one type parameterised by (coordinate, distribution), not a matrix of named restraints.

> **Method-name contract (reconciled with siblings):** the coordinate exposes
> `xi(x, pbc)` and `grad_xi(x, pbc)` (sub-spec -02); the distribution exposes
> `u(ξ, kt)` and `du_dxi(ξ, kt)` (sub-spec -03, and the same surface on
> `TabulatedProfile` from -04). This spec composes those exact names.

## Domain basis

Composition by chain rule: a profile restraint's penalty on a biased site at position x is U(ξ(x)); its gradient is ∇ₓU = (dU/dξ)·∇ξ(x), where ξ and ∇ξ come from the coordinate (-02) and U, dU/dξ from the distribution (-03/-04). Penalty by Boltzmann inversion U(ξ) = −kT·ln(ρ*(ξ)/ρ₀); ρ₀ shifts U by a constant only (no force effect), so it is absorbed into the distribution and never appears in the gradient.

§6.8 schedule coupling: the profile weight rides the existing radius soft-start so an early strong bias does not fight the overlap term; the restraint must be active during the per-type pre-compaction phase AND the joint phase. In molpack this is automatic — the global two-scale weights `scale`/`scale2` are set per phase and ramped with the radius. `ProfileRestraint` multiplies its penalty and gradient through `scale` (linear, matching how built-in linear-energy restraints are weighted) and holds NO mutable ramp state.

§6.1 per-site semantics: the bias is applied to a SELECTED SUBSET of atoms (a representative site, e.g. lipid head P), not every atom — else it over-constrains the molecule. Selection is handled OUTSIDE the restraint by the existing per-atom CSR assignment (`atoms … end atoms`); `ProfileRestraint` itself only evaluates the position it is handed.

References: [3] Martínez et al., PACKMOL, J Comput Chem 30:2157 (2009), DOI 10.1002/jcc.21224.

## Design

`ProfileRestraint` is an immutable value type composing the two upstream families:

- `pub struct ProfileRestraint { coordinate: Coordinate, distribution: Distribution, kt: F }`, where `Coordinate` is the coordinate enum from -02 (planar / radial / cylindrical / region-distance) and `Distribution` is the shared distribution surface from -03 and -04 (analytic Gaussian / erf / tanh / exponential via the `ProfilePenalty` façade, and tabulated via `TabulatedProfile`). Both are immutable value types, so `ProfileRestraint` is `Send + Sync` with no interior mutability. Derives `Debug`; a constructor `ProfileRestraint::new(coordinate, distribution, kt)` stores the three fields verbatim. (If the analytic `ProfilePenalty` and `TabulatedProfile` are distinct types, wrap them behind a small `Distribution` enum here so `u`/`du_dxi` dispatch is uniform.)
- `f(x, scale, scale2)` computes ξ = `coordinate.xi(x, pbc)`, U = `distribution.u(ξ, kt)`, and returns `scale * U`. The coordinate's inner clamp (-02) keeps ξ finite near a centre/axis; the distribution floor (-03/-04) keeps U finite in a zero-density bin. The penalty is finite by construction.
- `fg(x, scale, scale2, g)` computes ξ and ∇ξ = `coordinate.grad_xi(x, pbc)`, dU/dξ = `distribution.du_dxi(ξ, kt)`, then accumulates `g[k] += scale * (dU/dξ) * ∇ξ[k]` for k in 0..3 — additive into `g` with `+=`, never overwriting, mirroring the guarded-gradient shape of `InsideSphereRestraint` (src/restraint/geometric.rs, ex-`restraint.rs:288–298`) and `RegionRestraint` (src/region.rs:237–247). It returns `self.f(x, scale, scale2)` so the value is computed once per contract.
- `periodic_box()` is overridden to return the coordinate's box when the coordinate is anchored in a periodic box (so the pair-kernel minimum-image wrap stays consistent); otherwise it defers to the default `None`.
- Weight choice: profile penalty is a linear-energy term (U has energy units, not squared-distance), so it is multiplied by `scale` — not `scale2`, which the quadratic overlap/region terms use. This is documented inline.

The objective integration is unchanged: the CSR loop in `accumulate_constraint_value` (src/objective.rs:520–535) and `accumulate_constraint_gradient` (src/objective.rs:537–548) already iterate every restraint assigned to an atom and pass `sys.scale` / `sys.scale2`, so `ProfileRestraint` is picked up automatically with no edit there.

## Files to create or modify

- `src/restraint/profile/mod.rs` — add the `ProfileRestraint` struct, constructor, and `impl Restraint`, plus a `Distribution` dispatch enum if -03/-04 expose distinct types; add `#[cfg(test)] mod tests;`. (File created by `-02-coordinate`; this sub-spec appends the composed type and its re-export.)
- `src/restraint/profile/tests/compose.rs` (new) — unit + property tests for the composed restraint (value, full gradient matrix, guard/floor, scale linearity, periodic box).

## Tasks

- [x] Write failing test: planar-Gaussian `f` returns scale·(kT/2σ²)(ξ−μ)² and `fg` returns the same value while accumulating g additively into a pre-seeded gradient (src/restraint/profile/tests/compose.rs)
- [x] Write failing test: full coordinate×distribution analytic-vs-FD gradient matrix — planar/radial/cylindrical/region-distance × Gaussian/erf/tanh/exponential/tabulated (src/restraint/profile/tests/compose.rs)
- [x] Write failing test: finite `f` and `g` for a site on the radial centre / cylindrical axis (inside r_guard) and a site in a zero-density bin (src/restraint/profile/tests/compose.rs)
- [x] Write failing test: doubling `scale` doubles `f` and `g`, and `periodic_box()` mirrors the coordinate's box when anchored (src/restraint/profile/tests/compose.rs)
- [x] Implement `ProfileRestraint` struct, `new` constructor, `Distribution` dispatch enum, and `Debug` derive in src/restraint/profile/mod.rs
- [x] Implement `Restraint::f` (ξ = coordinate.xi, U = distribution.u, return scale·U) in src/restraint/profile/mod.rs
- [x] Implement `Restraint::fg` with chain-rule gradient accumulated `g += scale·(dU/dξ)·∇ξ`, returning `f`, guarded per the -02 clamp and -03/-04 floor in src/restraint/profile/mod.rs
- [x] Implement `Restraint::periodic_box` override deferring to the coordinate's anchored box in src/restraint/profile/mod.rs
- [x] Add docstring documenting weight choice (`scale`, linear) and units (kT) on `ProfileRestraint` in src/restraint/profile/mod.rs
- [x] Run full check + test suite (`cargo fmt`, `cargo clippy -- -D warnings`, `cargo test -p molcrafts-molpack --lib --tests`)

## Testing strategy

- Happy path: planar-Gaussian `f` equals the closed-form harmonic value scale·(kT/2σ²)(ξ−μ)² at several sample ξ; `fg` returns the identical value and accumulates g += scale·(kT/σ²)(ξ−μ)·n̂ into a non-zero pre-seeded gradient, proving additivity (no overwrite).
- Domain validation (gradient contract): for every coordinate×distribution pair in the §6.9 matrix (4 coordinates × 5 distributions), the analytic gradient from `fg` agrees componentwise with a central finite difference of `f` within a stated tolerance (e.g. ≤1e-6 relative). This is the binding gradient check for the whole family.
- Edge cases: a site exactly on the radial centre / cylindrical axis (inside r_guard) yields finite `f` and finite `g` (no Inf/NaN), exercising the -02 inner clamp; a site whose ξ lands in a zero-density tabulated bin yields finite `f`/`g`, exercising the -03/-04 floor.
- Invariants: doubling `scale` exactly doubles both `f` and the accumulated `g` (linear weighting, no ramp state); `periodic_box()` returns the coordinate's box iff the coordinate is anchored in one. `Send + Sync + Debug` are compile-time-asserted via a static bound in the test module.

## Out of scope

- Script parsing and lowering of profile restraints (`.inp` `atoms … end atoms` wiring) — deferred to sub-spec -06.
- Python / PyO3 exposure of `ProfileRestraint`.
- Iterative Boltzmann Inversion (IBI) refinement of ρ*.
- Composition-ratio coupling between multiple profile restraints.
