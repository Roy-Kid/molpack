---
title: Profile distribution penalties — closed-form Boltzmann inversion with Jacobian and floor
status: approved
created: 2026-06-12
---

# Profile distribution penalties — closed-form Boltzmann inversion with Jacobian and floor

## Summary
This sub-spec adds the analytic target-distribution penalties that turn a desired 1-D profile of a scalar reaction coordinate ξ into a soft energy `U(ξ)` and force `−dU/dξ` by Boltzmann inversion, `U(ξ) = −kT·ln(ρ*(ξ)/ρ₀)`. It supplies four closed-form distribution shapes — Gaussian, error-function step, hyperbolic-tangent step, and exponential decay — each exposing `u(ξ, kt)` and `du_dxi(ξ, kt)` from exact derivatives (no numerical differentiation), plus the two correctness safeguards from the domain basis: a density floor that caps the penalty at `U_max = −kT·ln(ρ_min/ρ₀)` so a zero-probability bin cannot become an infinite wall, and a shell-volume Jacobian / input-kind flag that divides a count histogram by `dV/dξ` before inversion so radial and cylindrical profiles are not corrupted by a spurious centre-weighting force. The module is a pure leaf: it knows only the scalar ξ — no coordinate geometry, no tabulated profiles, and no `Restraint` implementation (those are sibling sub-specs).

## Domain basis

Penalty by Boltzmann inversion of a target volumetric density ρ*(ξ):

- General: `U(ξ) = −kT·ln(ρ*(ξ)/ρ₀)`; `dU/dξ = −kT·ρ*'(ξ)/ρ*(ξ)`.
- `kT` is a FIXED energy-scale factor (not a simulated temperature) sharing the overlap-penalty's units; it sets bias steepness. `ρ₀` is an arbitrary reference entering `U` only as an additive constant ⇒ it does NOT affect `dU/dξ` or the force (it sets the energy zero only).

Closed forms (implement these exactly; do not finite-difference):

- Gaussian `ρ* ∝ exp(−(ξ−μ)²/2σ²)`: `U = (kT/2σ²)(ξ−μ)² + const`; `dU/dξ = (kT/σ²)(ξ−μ)`. EXACTLY a harmonic well, spring `k = kT/σ²`, centred at μ (harmonic tethering is the stiff/small-σ limit).
- Erf step `ρ* ∝ ½[1+erf((ξ−ξ₀)/(√2 w))]` (DEFAULT interface shape — capillary-wave-averaged physical liquid interface): `U = −kT·ln{½[1+erf((ξ−ξ₀)/(√2 w))]} + const`; `dU/dξ = −kT·(√(2/π)/w·exp(−(ξ−ξ₀)²/2w²)) / [1+erf((ξ−ξ₀)/(√2 w))]`. Falling side: flip the sign of the erf argument.
- Tanh step `ρ* ∝ ½[1+tanh((ξ−ξ₀)/w)]` (mean-field variant): `dU/dξ = −kT·(1/w)·sech²((ξ−ξ₀)/w) / [1+tanh((ξ−ξ₀)/w)]`. Expose both erf and tanh.
- Exponential `ρ* ∝ exp(−ξ/λ)`, ξ≥0: `U = (kT/λ)ξ + const`; `dU/dξ = kT/λ = CONSTANT` force toward small ξ. Valid only ξ≥0 (the inner clamp in -02 guarantees ξ never goes negative; this module assumes a non-negative ξ for the exponential and does not re-clamp).

Shell-volume Jacobian (§6.2 — load-bearing correctness point): Boltzmann inversion must be applied to a VOLUMETRIC density ρ*(ξ) [count/volume]. A count histogram n(ξ) [count/ξ] must be converted `ρ*(ξ) = n(ξ)/(dV/dξ)` BEFORE inversion. `dV/dξ`: planar = A (const); radial = 4π ξ²; cylindrical = 2π ξ L. If omitted in radial/cylindrical, `U` gains a spurious `+kT·ln(4π ξ²)` [radial] / `+kT·ln(2π ξ L)` [cyl] whose derivative adds a spurious `+2kT/ξ` [radial] / `+kT/ξ` [cyl] outward force diverging as ξ→0, over-weighting the centre/near-axis. RULE: the distribution input MUST declare density-vs-histogram; default-safe = treat input as a histogram and divide by the geometry's shell factor. (The geometry kind needed for `dV/dξ` is supplied by the caller / -05 compose; this module takes the shell factor as a `ShellJacobian` parameter.)

Divergence control (§6.3): replace ρ* by `max(ρ*, ρ_min)` ⇔ cap `U_max = −kT·ln(ρ_min/ρ₀)`. A zero-probability bin is unobservable, not infinitely forbidden; a +∞ wall has no finite gradient and stalls gencan. Pick `ρ_min ~ 1e-3–1e-2` of the peak.

References:
- [1] Boltzmann distribution ρ ∝ exp(−E/kT) ⇒ E = −kT·ln ρ + const.
- [3] Martínez et al., PACKMOL, J Comput Chem 30:2157 (2009), DOI 10.1002/jcc.21224.
- [4] erf vs tanh interfacial profiles (capillary-wave vs mean-field), arXiv:cond-mat/0605104.
- [5] Gouy–Chapman / linearized Poisson–Boltzmann exponential decay (Debye length λ_D).

## Design

New leaf module `src/restraint/profile/distribution.rs`, registered in `src/restraint/profile/mod.rs` (created by `-02-coordinate`; this sub-spec only appends `pub mod distribution;` + re-exports). All types are immutable value types deriving `Debug, Clone, Copy`, usable from `Send + Sync` contexts (no interior mutability, no globals).

Entities:

- `Distribution` — an enum of the four closed-form shapes:
  - `Gaussian { mu: F, sigma: F }`
  - `ErfStep { xi0: F, w: F, rising: bool }`
  - `TanhStep { xi0: F, w: F, rising: bool }`
  - `Exponential { lambda: F }`
  Methods `fn u(&self, xi: F, kt: F) -> F` and `fn du_dxi(&self, xi: F, kt: F) -> F` return the raw closed-form value/derivative (additive `ρ₀`/`const` reference dropped to zero — `ρ₀` does not affect the force, and the energy zero is immaterial to gencan). `rising = false` flips the erf/tanh argument sign for a falling interface.
- `ShellJacobian` — an enum selecting `dV/dξ` for histogram→density conversion:
  - `Planar` (`dV/dξ = const`, factor 1 — no centre-weighting term)
  - `Radial` (`dV/dξ ∝ ξ²`)
  - `Cylindrical { length: F }` (`dV/dξ ∝ ξ`)
  Method `fn log_shell_correction(&self, xi: F) -> F` returns the additive `ln(dV/dξ)` term whose ξ-derivative is the spurious force, and `fn d_log_shell(&self, xi: F) -> F` returns that derivative (`0` planar, `2/ξ` radial, `1/ξ` cyl). These are used only when the input is a histogram.
- `InputKind` — `VolumetricDensity | CountHistogram`. The default-safe constructor / `Default` impl yields `CountHistogram` (matches the domain RULE).
- `DensityFloor` — holds `rho_min` and the precomputed cap `u_max = −kt·ln(rho_min/rho0)`; exposes `fn cap(&self, u: F) -> F` clamping `u` to `u_max` and a companion predicate so the derivative can be zeroed when the cap is active (a capped bin contributes no force). Stores the `kt`/`rho0` it was built against, or takes them at construction.
- A small façade `ProfilePenalty { dist: Distribution, jacobian: ShellJacobian, input_kind: InputKind, floor: DensityFloor }` (immutable, `Debug, Clone, Copy` — `Cylindrical` carries an `F`, still `Copy`) exposing the composed `fn u(&self, xi: F, kt: F) -> F` and `fn du_dxi(&self, xi: F, kt: F) -> F`: it inverts the closed-form `ρ*`, divides by the shell factor when `input_kind == CountHistogram` (subtracting `log_shell_correction` from `U` / `d_log_shell` from `dU/dξ`), then applies the floor cap (and zeroes the derivative where capped). For `VolumetricDensity` the shell correction is skipped — the input is already a density.

Lifecycle/ownership: every type is a plain `Copy` value built by the caller (-05 compose supplies the geometry-derived `ShellJacobian` and the `kt`); this module owns no buffers and performs no allocation. Functions are referentially transparent in (ξ, kt). No coordinate vectors enter this module — ξ is already reduced to a scalar by -02.

## Files to create or modify
- `src/restraint/profile/distribution.rs` (new) — `Distribution`, `ShellJacobian`, `InputKind`, `DensityFloor`, `ProfilePenalty`, their `u`/`du_dxi` closed forms, and the in-file `#[cfg(test)]` module.
- `src/restraint/profile/mod.rs` — add `pub mod distribution;` and re-export the public types (file created by `-02-coordinate`; this sub-spec only appends the declaration / re-export line).

## Tasks
- [ ] Write failing tests for closed-form `du_dxi` vs central finite difference, all four distributions (src/restraint/profile/distribution.rs)
- [ ] Write failing tests for Gaussian≡harmonic (U, dU/dξ vs (kT/2σ²)(ξ−μ)², (kT/σ²)(ξ−μ)) and Exponential≡constant-force (dU/dξ == kT/λ) (src/restraint/profile/distribution.rs)
- [ ] Write failing tests for histogram-vs-density Jacobian equivalence on a radial profile and for the spurious-force case when the shell factor is omitted (src/restraint/profile/distribution.rs)
- [ ] Write failing tests for the density floor: zero-density input yields finite U == U_max and finite (non-NaN/Inf) du_dxi (src/restraint/profile/distribution.rs)
- [ ] Implement `Distribution` enum with exact `u`/`du_dxi` closed forms for Gaussian, ErfStep, TanhStep, Exponential (src/restraint/profile/distribution.rs)
- [ ] Implement `ShellJacobian` (`log_shell_correction`, `d_log_shell`) and `InputKind` with `CountHistogram` default-safe (src/restraint/profile/distribution.rs)
- [ ] Implement `DensityFloor` cap (`U_max = −kT·ln(ρ_min/ρ₀)`) and the `ProfilePenalty` façade composing inversion, Jacobian division, and floor (src/restraint/profile/distribution.rs)
- [ ] Add module-level docstrings with units and a `pub mod distribution;` + re-export in src/restraint/profile/mod.rs
- [ ] Run full check + test suite (cargo fmt, cargo clippy -- -D warnings, cargo test -p molcrafts-molpack --lib --tests)

## Testing strategy
- Happy path: each `Distribution.u`/`du_dxi` returns the documented closed form at representative interior ξ (Gaussian near and far from μ; erf/tanh below, at, and above ξ₀ on both rising and falling sides; exponential at several ξ≥0).
- Domain validation — analytic-vs-numeric: central finite difference of `u(ξ, kt)` matches `du_dxi(ξ, kt)` within a stated tolerance (e.g. relative ≤ 1e-6 at step h≈1e-5, away from the floor-capped region) for all four shapes (ac-001).
- Domain validation — closed-form identities: Gaussian `u`/`du_dxi` equal `(kT/2σ²)(ξ−μ)²` and `(kT/σ²)(ξ−μ)` to floating tolerance; exponential `du_dxi` equals the constant `kT/λ` independent of ξ (ac-002).
- Domain validation — Jacobian (two-sided §6.2): build a count histogram n(ξ) and the volumetric density ρ*(ξ)=n(ξ)/(4π ξ²) describing the SAME radial profile; assert `ProfilePenalty` with `CountHistogram + Radial` and with `VolumetricDensity` give equal `U(ξ)` (to floating tolerance); assert that omitting the shell factor (treating the histogram as a density) gives a demonstrably different `U` whose `dU/dξ` differs by the predicted `+2kT/ξ` (ac-003).
- Edge cases / divergence control: a zero-density region (ρ* below ρ_min, or an erf/tanh tail driven to numerical zero) yields finite `U` capped at `U_max` and finite `du_dxi` — no `Inf`/`NaN`; rising/falling flag and `Cylindrical { length }` factor exercised; exponential evaluated only at ξ≥0 (ac-004).
- Smoke: module compiles, types are `Debug` + `Copy`; covered implicitly by the full check + test suite.

## Out of scope
- Coordinate geometry — projecting a 3-D position onto the scalar ξ and the inner clamp (sub-spec -02).
- Tabulated / spline target profiles (sub-spec -04).
- The `Restraint` trait implementation that wires ξ-projection, this penalty, and gradient accumulation into the packer (sub-spec -05 compose).
- Script (`.inp`) parsing of the new keyword surface and the Python binding (later sub-specs in the chain).
- Choosing ρ_min adaptively from the observed peak — the caller passes ρ_min; auto-derivation is deferred.
