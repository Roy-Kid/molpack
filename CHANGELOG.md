# Changelog

All notable changes to `molcrafts-molpack` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- `OutsideEllipsoidRestraint` (Packmol kind 9): `f` now multiplies the squared
  penalty by `scale2`, matching `fg`'s gradient and the documented "quadratic
  penalty group" convention. The previous transcription left `f` and `fg`
  100× out of phase at the default `scale2 = 0.01`, causing the optimizer to
  see a gradient much flatter than the function value reported.

### Added
- Three-phase packing algorithm (per-type init → geometric pre-fit → main loop with `movebad` heuristic)
- GENCAN quasi-Newton optimizer with bound constraints
- 14 concrete `*Restraint` types covering Packmol kinds 2–15
- `Region` trait with `And` / `Or` / `Not` combinators and `RegionRestraint` lift
- `Handler` trait with built-ins: `NullHandler`, `ProgressHandler`, `EarlyStopHandler`, `XYZHandler`
- `Relaxer` trait with `TorsionMcRelaxer` built-in
- Per-axis periodic boundary conditions via `InsideBoxRestraint::new(min, max, [bool; 3])`
- `pbc` keyword in `.inp` scripts and `Molpack::with_periodic_box`
- Optional `rayon` feature for parallel objective evaluation
- Optional `cli` feature gating the `molpack` binary
- Python bindings via PyO3 + maturin (`molcrafts-molpack` on PyPI)
- Five canonical Packmol workloads as runnable examples (mixture, bilayer, spherical, interface, solvprotein)
- Criterion benchmark suite (e2e + hot-path microbenchmarks)
- Rustdoc chapters: `getting_started`, `concepts`, `architecture`, `extending`

[Unreleased]: https://github.com/MolCrafts/molpack/compare/HEAD...HEAD
