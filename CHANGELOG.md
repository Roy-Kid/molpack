# Changelog

All notable changes to `molcrafts-molpack` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial implementation of the three-phase packing algorithm (per-type init → geometric pre-fit → main loop with `movebad` heuristic)
- GENCAN quasi-Newton optimizer with bound constraints (`gencan` module)
- 14 concrete `*Restraint` types covering Packmol kinds 2–15
- `Region` trait with `And` / `Or` / `Not` combinators and `FromRegion` lift
- `Handler` trait with four built-ins: `NullHandler`, `ProgressHandler`, `EarlyStopHandler`, `XYZHandler`
- `Relaxer` trait with `TorsionMcRelaxer` built-in
- Optional `rayon` feature for parallel objective evaluation
- Optional `filesystem` feature forwarded from `molrs_io`
- Python bindings via PyO3 + maturin (`molcrafts-molpack` on PyPI)
- Five canonical Packmol workloads as examples (mixture, bilayer, spherical, interface, solvprotein)
- Criterion benchmark suite (five benches covering end-to-end and hot-path microbenchmarks)
- Rustdoc chapters: `getting_started`, `concepts`, `architecture`, `extending`

[Unreleased]: https://github.com/MolCrafts/molpack/compare/HEAD...HEAD
