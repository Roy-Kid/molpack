# Changelog

All notable changes to `molcrafts-molpack` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-06-10

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
- Five canonical Packmol workloads as runnable examples (mixture, bilayer, spherical, interface, solvprotein)
- Criterion benchmark suite (e2e + hot-path microbenchmarks)
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
- Relaxer Monte-Carlo loop reuses a single trial buffer (`mem::swap` on accept) instead of cloning per step.
- `self_avoidance_penalty` evaluates the cheap distance gate before the exclusion-set hash probe.
- CI gate is tests-only; formatting and lint (`fmt`, `clippy`, `ruff`, `ty`) are enforced by pre-commit hooks.
- `molrs-core` / `molrs-io` pin bumped to `0.0.18` (`MolGraph` → `Atomistic` API).

[Unreleased]: https://github.com/MolCrafts/molpack/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/MolCrafts/molpack/releases/tag/v0.1.0
