# Changelog

All notable changes to `molcrafts-molpack` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- All 14 restraint kinds are now reachable from `.inp` scripts and the Python
  API (previously only 5 were wired through each frontend). Newly available:
  `inside`/`outside` `cube`/`ellipsoid`/`cylinder`, `outside box`, and
  `over`/`below gaussian`; Python gains the matching `*Restraint` classes. A
  single keyword→spec table (parser) and `spec_to_restraint` map (build) are
  the one place a new kind is wired.

### Changed
- `.inp` parser now rejects `atoms 0` (atom indices are 1-based) and `inside`/
  `outside sphere` with `radius <= 0`; Python `Target`/`Packer` reject
  `count == 0`. Previously these were silently accepted (`atoms 0` mapped to the
  wrong atom; a non-positive radius fed nonsense into the distance math).

### Fixed
- Python `PackResult.positions`/`frame`/`elements` now raise `ValueError` on a
  malformed packed frame instead of panicking (which aborted the interpreter at
  the FFI boundary).
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
