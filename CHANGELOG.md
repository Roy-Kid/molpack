# Changelog

All notable changes to `molcrafts-molpack` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed (breaking)

- In-loop geometry: deleted `Relaxer` / `TorsionMcRelaxer` / `LBFGSRelaxer` /
  `Target::with_relaxer`. Single entry `Molpack::with_optimizer(OptimizeSelect,
  impl molrs::optimize::Optimizer)` (`ff` feature). Selection is
  `OptimizeSelect::per_copy|joint(names)` with optional `.with_environment(rcut)`.
  Callers build molrs optimizers (`LBFGS`, `SoftSpec::into_optimizer`) or
  molpack `TorsionMcOptimizer`. No deprecated shims.
- Post-pack relaxation: deleted the `molpack.relaxer` module and its
  `LAMMPSRelaxer` façade over `molpy.engine.LAMMPSEngine`. molpack no longer
  ships any relaxer, in-loop or post-pack; drive an MD engine directly from the
  packed frame if a settle step is wanted.

### Fixed

- Python tests construct `molrs.Frame.from_dict` with both `"blocks"` and
  `"meta"` (required by molrs 0.9.x).

### Changed

- Pinned runtime deps to `molcrafts-molrs==0.9.3` and `molcrafts-molpy==0.9.3`;
  CI checks out `MolCrafts/molrs@v0.9.3` for path builds.
- Python tests via `tox -c python -e py` (config in `python/pyproject.toml`
  `[tool.tox]`; isolated non-editable wheel; no `scripts/`). Git hooks use
  **prek**; pre-push mirrors CI.

### Changed (earlier)

- Updated the Rust integration for molrs's chemical-perception reorganization:
  torsion relaxers now use `molrs::perceive::rotatable` after the removal of the
  compatibility `molrs::chem` module.
- Aligned the force-field relaxer documentation and regression coverage with
  molrs's generic MMFF94 pipeline (`typify` → typed frame →
  `intramolecular_pairs` → `ForceField::to_potentials`) and the renamed
  `MMFF94Typifier` API.
- Updated nonbonded force-field coverage for molrs's explicit `coul/cut`
  contract: force fields now provide their Coulomb constant and dielectric
  instead of relying on kernel-owned physical defaults.
- Updated post-pack relaxation references to molpy's `LAMMPSEngine` API.

## [0.1.1] - 2026-07-08

### Changed
- Bumped the `molcrafts-molrs` dependency (Rust crate + `molrs-ffi` FFI layer and
  the Python runtime pins) from 0.6.0 to 0.7.0, and the Python runtime
  `molcrafts-molpy` pin from 0.6.0 to 0.7.0. No molpack API change — a
  dependency-tracking release against the unified molrs/molpy 0.7.0 line.
- Docs build now pins the shared `molcrafts-zensical-theme` at 0.1.1.

## [0.1.0] - 2026-07-04

Inaugural release of `molcrafts-molpack`.

### Added
- Three-phase packing algorithm (per-type init → geometric pre-fit → main loop with `movebad` heuristic)
- GENCAN quasi-Newton optimizer with bound constraints
- 14 concrete `*Restraint` types covering Packmol kinds 2–15, each with an analytic `f` energy and a matching `fg` gradient (verified by finite-difference unit tests)
- `Region` trait with `And` / `Or` / `Not` combinators and `RegionRestraint` lift
- `Handler` trait with built-ins: `NullHandler`, `ProgressHandler`, `EarlyStopHandler`, `XYZHandler`
- `Relaxer` trait with `TorsionMcRelaxer` built-in
- `ForceFieldRelaxer` (feature `ff`): a `Relaxer` that relaxes a flexible
  molecule's internal geometry during packing by energy-minimizing it under a
  caller-supplied molrs force-field `Potential` (via molrs's L-BFGS optimizer),
  accepting the relaxed conformer only when it does not worsen the packer
  objective. The `ff` feature enables `molrs/ff` and re-exports the `Potential`
  trait at the crate root.
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
- Five canonical Packmol workloads as runnable examples (mixture, bilayer, spherical, interface, solvprotein), plus the `mt_scaling` measurement harness (parallel speed-up-vs-size sweep)
- Criterion regression benchmark suite (`cargo bench --benches`, no `io` feature
  required — each synthesizes its own geometry in-process): the fused
  pair/objective kernel (`compute_f` / `compute_fg`), objective dispatch
  (`via_inherent` vs `via_dyn`), one `run_iteration` step, restraint evaluation
  (`InsideBox` / `InsideSphere` `f`/`fg`), and a tiny end-to-end `Molpack::pack`
  mixture. Each uses `sample_size(10)` + ~500 ms measurement so the whole run is
  cheap; `.github/workflows/bench.yml` tracks their history on canonical pushes
  via `benchmark-action/github-action-benchmark`.
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
- Parallel objective gradient uses a half-stencil scheme with parallel projection (under the `rayon` feature), halving redundant pair evaluations on the hot path.
- Restraint code split from a single `src/restraint.rs` into the `src/restraint/` module (`geometric/`, `profile/`) as the type count grew.
- Relaxer Monte-Carlo loop reuses a single trial buffer (`mem::swap` on accept) instead of cloning per step.
- `self_avoidance_penalty` evaluates the cheap distance gate before the exclusion-set hash probe.
- CI runs the test suites plus the Packmol regression (on non-PR events); formatting and lint (`fmt`, `clippy`, `ruff`, `ty`) are enforced by pre-commit hooks.
- Depends on the unified `molcrafts-molrs` 0.6.0 crate: `core` is always-on and re-exported at the crate root, while `io` and `ff` are feature-gated modules that forward to `molrs/io` and `molrs/ff`.

[0.1.1]: https://github.com/MolCrafts/molpack/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/MolCrafts/molpack/releases/tag/v0.1.0
