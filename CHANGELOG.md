# Changelog

All notable changes to `molcrafts-molpack` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `ForceFieldRelaxer` (feature `ff`): a `Relaxer` that relaxes a flexible
  molecule's internal geometry during packing by energy-minimizing it under a
  caller-supplied molrs force-field `Potential` (via molrs's L-BFGS optimizer),
  accepting the relaxed conformer only when it does not worsen the packer
  objective. The `ff` feature enables `molrs/ff` and re-exports the `Potential`
  trait at the crate root.

### Changed
- Migrated the molrs dependency from the separate `molcrafts-molrs-core` +
  `molcrafts-molrs-io` crates to the unified `molcrafts-molrs` crate (v0.1.1):
  `core` is always-on (re-exported at the crate root), while `io` and `ff` are
  feature-gated modules. The molpack `io` feature now forwards to `molrs/io`, and
  `molrs_io::` import paths become `molrs::io::`. Behavior-preserving — the
  Packmol regression suite is unchanged.

### Removed
- The `validate_packed` example harness (external-file accuracy scorer for the
  molpack-jcc figures); dropped from `examples/` and `Cargo.toml`.

## [0.1.0] - 2026-06-13

Inaugural release of `molcrafts-molpack`.

### Added
- Three-phase packing algorithm (per-type init → geometric pre-fit → main loop with `movebad` heuristic)
- GENCAN quasi-Newton optimizer with bound constraints
- 14 concrete `*Restraint` types covering Packmol kinds 2–15, each with an analytic `f` energy and a matching `fg` gradient (verified by finite-difference unit tests)
- `Region` trait with `And` / `Or` / `Not` combinators and `RegionRestraint` lift
- `Handler` trait with built-ins: `NullHandler`, `ProgressHandler`, `EarlyStopHandler`, `XYZHandler`
- `Relaxer` trait with `TorsionMcRelaxer` built-in
- Per-axis periodic boundary conditions via `InsideBoxRestraint::new(min, max, [bool; 3])`
- `pbc` keyword in `.inp` scripts and `Molpack::with_periodic_box`
- `.inp` parser now covers 12 restraint forms — `inside`/`outside` for `box`, `cube`, `sphere`, `ellipsoid`, and `cylinder`, plus `over`/`below plane` and `over`/`below xygauss` (with positive-axis / non-zero-axis input validation)
- Public extension surface re-exported at the crate root: `Objective`, `Constraints`, `EvalMode`, `EvalOutput`
- Optional `rayon` feature for parallel objective evaluation
- Optional `cli` feature gating the `molpack` binary
- Python bindings via PyO3 + maturin (`molcrafts-molpack` on PyPI)
- `PackResult.frame` returns a topology-complete `molrs.Frame`: each target's source-frame topology (bonds/angles/dihedrals/impropers) is replayed onto the packed coordinates with indices offset per copy, `id`/`mol_id` regenerated, and `frame.box` stamped from the periodic box — adopt into molpy with `molpy.Frame.from_dict(...)`. Falls back to a coordinates-only `atoms` block for `.inp` script packing.
- Profile-distribution restraint family (`ProfileRestraint`): composes a reaction
  coordinate (planar / radial / cylindrical) with a target distribution
  (`gaussian`, `erf`, `tanh`, `exponential`, or a `tabulated` density/histogram via
  monotone-cubic interpolation + Boltzmann inversion), each with an analytic
  Jacobian. Exposed as the `.inp` `profile` keyword and lowered through
  `Script::lower` + `StructurePlan::apply`.
- Python restraint extensibility: user code can subclass the restraint base in
  pure Python and have it driven through the Rust optimizer (parallel-unsafe
  restraints are serialized automatically); see
  `python/examples/pack_profile_monolayer.py` and the
  `profile_monolayer_extensibility.md` write-up.
- Five canonical Packmol workloads as runnable examples (mixture, bilayer, spherical, interface, solvprotein), plus two measurement harnesses: `mt_scaling` (parallel speed-up-vs-size sweep) and `validate_packed` (scores an external packed file against its `.inp` through molpack's validator)
- Criterion end-to-end benchmark (`cargo bench --bench pack_end_to_end --features io`): drives all five canonical workloads through `Molpack::pack` as a catastrophic-regression alarm
- Rustdoc chapters: `getting_started`, `concepts`, `architecture`, `extending`

### Fixed
- `OutsideEllipsoidRestraint` (Packmol kind 9): `f` now multiplies the squared
  penalty by `scale2`, matching `fg`'s gradient and the documented "quadratic
  penalty group" convention. The previous transcription left `f` and `fg`
  100× out of phase at the default `scale2 = 0.01`, causing the optimizer to
  see a gradient much flatter than the function value reported.
- `movebad` worst-atom selection uses `f64::total_cmp` instead of `partial_cmp().unwrap()`, removing a NaN-triggered panic path.
- `initial` cell-grid sizing clamps total cell count (atom-count aware, hard ceiling) so a degenerate fallback box can no longer trigger an out-of-memory allocation.

### Changed
- Parallel objective gradient now uses a half-stencil scheme with parallel projection (under the `rayon` feature), halving redundant pair evaluations on the hot path.
- Restraint code split from a single `src/restraint.rs` into the `src/restraint/` module (`geometric/`, `profile/`) as the type count grew.
- Relaxer Monte-Carlo loop reuses a single trial buffer (`mem::swap` on accept) instead of cloning per step.
- `self_avoidance_penalty` evaluates the cheap distance gate before the exclusion-set hash probe.
- CI runs the test suites plus the Packmol regression (on non-PR events); formatting and lint (`fmt`, `clippy`, `ruff`, `ty`) are enforced by pre-commit hooks.
- `molrs-core` / `molrs-io` pinned to `0.1.0`.

[Unreleased]: https://github.com/MolCrafts/molpack/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/MolCrafts/molpack/releases/tag/v0.1.0
