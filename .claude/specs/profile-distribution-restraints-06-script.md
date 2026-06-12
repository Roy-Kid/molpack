---
title: Wire profile-distribution restraints into the .inp script layer
status: done
created: 2026-06-12
---

# Wire profile-distribution restraints into the .inp script layer

## Summary

This capstone sub-spec exposes the `ProfileRestraint` family (built in sub-spec `-05-compose`) through molpack's Packmol-compatible `.inp` script grammar, so users can bias selected sites toward a target spatial distribution without touching Rust. A new lowercase `profile` keyword — named after the distribution it imposes, with geometry as a parameter — parses inside a `structure` block or an `atoms … end atoms` sub-block, lowers to a concrete `ProfileRestraint` on the corresponding `Target`, and rides the existing scale/scale2 soft-start automatically. The Python wheel inherits the feature for free through `Script::lower`, with no new PyO3 surface. End-to-end packing tests land the four user-facing acceptance criteria, proving the distribution actually shapes the packed configuration.

## Domain basis

The profile-restraint family biases selected sites toward a target distribution ρ\*(ξ) by Boltzmann inversion of a one-dimensional collective coordinate ξ:

```
U(ξ) = −kT · ln( ρ*(ξ) / ρ₀ )
```

A restraint composes a **coordinate** ξ (planar signed distance, radial distance, cylindrical radius, or region-distance) with a **distribution** ρ\* (Gaussian, erf, tanh, exponential, or tabulated). Restraints are named after the distribution; the geometry is a parameter (design discussion §2, §5).

Load-bearing rules this keyword must honor:

- **§6.1 per-site subset** — the bias targets a representative-site subset of atoms, exposed through the *existing* `atoms <idx…> … end atoms` sub-block; no new selection syntax is introduced.
- **§6.2 density-vs-count** — for a tabulated input the user must declare whether values are a volumetric density or a count histogram. Default-safe is `histogram`: the shell Jacobian (∝ ξ² for radial, ∝ ξ for cylindrical) is divided out before inversion so a uniform-count target does not produce a spurious center-seeking bias.
- **§6.7 superposition** — a two-component asymmetric interface is expressed as two *independent* profile restraints (opposite erf sides, or Gaussians at +d and −d) that superpose; no special composite keyword is needed.
- **§6.8 soft-start** — the restraint weight rides the existing scale/scale2 soft-start automatically once lowered onto a `Target`; this sub-spec adds no new weighting knob.

Reference: [3] Martínez, Andrade, Birgin, Martínez, *PACKMOL: A package for building initial configurations for molecular dynamics simulations*, J. Comput. Chem. **30**, 2157 (2009), DOI 10.1002/jcc.21224.

## Design

Extend the script layer along its single, already-established spec→restraint pipeline:

- **`RestraintSpec::Profile` variant** (`src/script/parser.rs`, enum at line 87). Carries: distribution kind (an enum mirroring `-05`'s distribution set — gaussian / erf / tanh / exponential / tabulated), geometry parameters (coordinate kind plus its numeric args — plane normal+origin, radial center, cylinder center+axis, region-distance handle), the distribution's own numeric params (e.g. `mu`, `sigma` for gaussian), an `input_kind` flag (`Density | Histogram`, default `Histogram`) used only by the tabulated path, and an optional table source (file path or inline nodes). The variant is `Debug + Clone` like its siblings.

- **`parse_profile()`** (`src/script/parser.rs`, mirroring `parse_inside` at line 401). Tokenises a Packmol-style whitespace line `profile <distribution> <geometry> <geometry-args> <distribution-params> [density|histogram]`, validates token counts and keyword spelling via the existing `parse_f64` / `parse_vec3` / `parse_err` helpers, and returns `RestraintSpec::Profile`. Dispatched from a new `"profile"` arm added to **both** the `State::InStructure` match (around parser.rs:281–319) and the `State::InAtoms` match (around parser.rs:338) so the keyword works whole-molecule and per-atom-subset identically to `inside`/`outside`/`over`/`below`. If parser.rs exceeds its 800-LOC budget after these additions, the profile-parsing helpers are extracted into a sibling `src/script/parser_profile.rs` re-exported from `parser.rs` (high-cohesion split, no public-surface change).

- **Lowering arm** in `restraint_from_spec` (`src/script/build.rs`, line 157). A new `RestraintSpec::Profile { … }` arm constructs the `-05` `ProfileRestraint` from the spec's coordinate + distribution + input_kind, boxed as `Box<dyn Restraint>`. Because `restraint_from_spec` is the single source of truth shared by `apply_mol_restraint` (build.rs:206) and `apply_atom_group` (build.rs:210), both the whole-molecule and `atoms`-subset paths get profile support from this one arm. The `ProfileRestraint` type and its constructor are imported into build.rs's `use crate::{…}` block.

- **Naming** — no identifier contains "packmol"; the keyword is the distribution name, geometry is a parameter (`profile gaussian plane …`). The grammar is Packmol-style: lowercase keyword, whitespace-delimited tokens, optional trailing `density|histogram` flag.

- **Python parity** — the wheel is built without `io`; the keyword flows through `Script::lower` → `StructurePlan::apply`, so Python callers reach it with zero new PyO3 type. A Python keyword-parity test is added only if the wheel test layer already asserts keyword coverage.

A measured-profile test helper (bin packed site positions along ξ, normalize, compare to the target curve within tolerance) is added to the e2e test file; it is the instrument behind acceptance criteria ac-003 and ac-004.

## Files to create or modify

- `src/script/parser.rs` — add `RestraintSpec::Profile` variant, `parse_profile()`, and `"profile"` dispatch arms in `State::InStructure` and `State::InAtoms`; add parser-level round-trip unit tests.
- `src/script/build.rs` — add the `RestraintSpec::Profile` arm in `restraint_from_spec`; import `ProfileRestraint`.
- `src/script/parser_profile.rs` (new) — *only if* parser.rs crosses its LOC budget; houses `parse_profile()` and its sub-helpers, re-exported from `parser.rs`.
- `tests/profile_pack_e2e.rs` (new) — end-to-end packing tests for the four user-facing criteria plus the measured-profile helper, mirroring `tests/restraint_pack_e2e.rs`.
- `python/tests/test_script_keywords.py` — add a `profile`-keyword parity assertion *only if* this file already asserts keyword coverage; otherwise untouched.

## Tasks

- [x] Write failing parser round-trip tests for the `profile` keyword (`src/script/parser.rs` `#[cfg(test)]`): `profile gaussian plane <nx ny nz> <x0 y0 z0> mu <μ> sigma <σ>` parses to `RestraintSpec::Profile` with correct distribution/geometry/params and default `Histogram`; an explicit `density` flag flips `input_kind`; the same line inside `atoms <idx…> … end atoms` lands on the atom group's restraints
- [x] Write failing lowering test (`src/script/build.rs` `#[cfg(test)]` or `tests/profile_pack_e2e.rs`) asserting a parsed profile spec lowers to a `ProfileRestraint` whose coordinate + distribution match the keyword args, checked via its `f`/`fg` at known points
- [x] Write failing end-to-end pack tests in `tests/profile_pack_e2e.rs` (new) with a measured-profile helper: single-component planar Gaussian → ρ(z) within tolerance; radial count-histogram → in tolerance only with the ξ² Jacobian and out of tolerance without it; two-component opposite-erf → asymmetric per-component layout; zero-density-region target → completes with finite forces, no NaN
- [x] Add the `RestraintSpec::Profile` variant to the enum in `src/script/parser.rs` (carrying distribution kind, geometry params, `input_kind`, optional table source)
- [x] Implement `parse_profile()` in `src/script/parser.rs` (extracting to `src/script/parser_profile.rs` if the file crosses 800 LOC) and wire `"profile"` dispatch arms into `State::InStructure` and `State::InAtoms`
- [x] Implement the `RestraintSpec::Profile` arm in `restraint_from_spec` (`src/script/build.rs`), importing `ProfileRestraint` and building it from the lowered coordinate + distribution + input_kind
- [ ] Add a `profile`-keyword parity assertion in `python/tests/test_script_keywords.py` only if that file already asserts keyword coverage *(N/A — `python/tests/test_script_keywords.py` does not exist; the keyword flows through `Script::lower` with no new PyO3 surface)*
- [x] Run full check + test suite (`cargo fmt`, `cargo clippy -- -D warnings`, `cargo test -p molcrafts-molpack --lib --tests`)

## Testing strategy

- **Happy path (parser)** — `profile gaussian plane 0 0 1 0 0 0 mu 10 sigma 2` round-trips to `RestraintSpec::Profile` with the expected fields; default `input_kind == Histogram`.
- **Happy path (lowering)** — the parsed spec lowers to a `ProfileRestraint` whose `f`/`fg` at a few known ξ values match the analytic Gaussian/erf form within `1e-9`.
- **Edge cases** — missing/extra tokens raise a `ScriptError` with line number; an unknown distribution or geometry sub-keyword is rejected; the `density|histogram` flag is case-insensitive and optional; an `atoms`-scoped `profile` applies to exactly the listed (1-based → 0-based) indices and no others.
- **Domain validation (ac-003)** — single-component planar Gaussian packs to a measured ρ(z) within stated tolerance; radial count-histogram matches its target *only* with the shell-volume (ξ²) Jacobian and is demonstrably out of tolerance without it (two-sided, §6.2).
- **Domain validation (ac-004)** — two-component opposite-erf produces asymmetric per-component measured profiles (each component biased to its own side, §6.7); a zero-density-region target packs to completion with finite forces and no NaN in positions or gradients.
- **Suite** — `cargo test -p molcrafts-molpack --lib --tests` stays green; new e2e tests stay in the fast tier (tiny atom counts, capped `max_loops`, seeded) at ≤ 200 ms each like `restraint_pack_e2e.rs`.

## Out of scope

- New site-selection syntax — the `atoms … end atoms` sub-block is reused verbatim.
- Per-molecule / center-of-mass biasing (only per-site atom subsets here).
- Composition-ratio coupling between components and Iterative Boltzmann Inversion (IBI) refinement.
- Runtime / MD biasing — this restraint shapes the *initial* packed configuration only.
- New PyO3 types — Python inherits the keyword through `Script::lower` with no binding change.
