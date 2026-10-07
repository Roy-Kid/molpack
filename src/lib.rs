//! # molpack
//!
//! Packmol-grade molecular packing in pure Rust. Produces a non-overlapping
//! arrangement of N molecule types with copy counts and geometric restraints,
//! using a faithful port of Packmol's GENCAN-driven three-phase algorithm
//! (Martínez et al. 2009). Correctness is checked against Packmol's reference
//! output for five canonical workloads.
//!
//! This crate was split out of the molrs workspace in 2026 and is now
//! maintained independently. It depends on the unified `molcrafts-molrs` crate
//! for shared data structures (always-on `core`) and, behind the `io` feature,
//! its file I/O module.
//!
//! ## Documentation map
//!
//! This crate is documented in four dedicated modules; start with
//! [`getting_started`] if you are new.
//!
//! - [`getting_started`] — install, load or build a molecule template,
//!   declare a target with a spatial restraint, run one pack, save the
//!   result. Written against the Python package, which is the shortest
//!   path from a loaded structure to a packed box.
//! - [`concepts`] — every abstraction defined in one place: `AtomRestraint`,
//!   the molrs `Region` lift, `Callback`, `Objective`, `Target`, `PackEngine`,
//!   `PackSystem`; the scope equivalence law; the two-scale contract;
//!   the direction-3 extension pattern.
//! - [`architecture`] — module map, dependency graph, core-type
//!   relationships, full `pack()` lifecycle diagram, hot-path
//!   `evaluate()` walkthrough, invariants, design decisions.
//! - [`extending`] — tutorials for writing your own `AtomRestraint` /
//!   molrs `Region` / `Callback` and for binding an in-loop optimizer; testing
//!   discipline; common pitfalls; contributing flow.
//!
//! Reference material (not rustdoc):
//!
//! - [Packmol parity](https://github.com/MolCrafts/molpack/blob/master/docs/packmol_parity.md)
//!   — kind-number ↔ Rust struct mapping with Fortran pointers.
//!
//! ## Quick example
//!
//! ```rust,no_run
//! use std::sync::Arc;
//! use molpack::{GencanPack, PackEngine, RegionRestraint, Target};
//! use molrs::core::Cuboid;
//! use ndarray::array;
//!
//! let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
//! let radii     = [1.52, 1.20, 1.20];
//!
//! // Geometry is a molrs region; molpack lifts it to "stay inside".
//! let cube = Cuboid::new(array![0.0, 0.0, 0.0], array![40.0, 40.0, 40.0]);
//! let target = Target::from_coords(&positions, &radii, 100)
//!     .with_name("water")
//!     .with_restraint(RegionRestraint(Arc::new(cube)));
//!
//! let result = GencanPack::new()
//!     .with_tolerance(2.0)
//!     .with_precision(0.01)
//!     .with_seed(42)
//!     .run(&[target], 200)?;
//!
//! let natoms = result.natoms();
//! println!("packed {natoms} atoms");
//! # Ok::<(), molpack::PackError>(())
//! ```
//!
//! ## Public surface at a glance
//!
//! Every public item has exactly one path. The crate root carries the
//! vocabulary below; three namespaces carry the rest — [`context`] (the
//! per-atom layout a custom objective or callback reads off a
//! [`PackSystem`]), [`grow`] (growth configuration and priors) and
//! [`script`] (the `.inp` loader). molrs types (`Frame`, `SimBox`, regions,
//! the [`molrs::optimize::Optimizer`] trait, `F`) are named at their molrs
//! home; molpack does not re-export them.
//!
//! | Category | Items |
//! |---|---|
//! | Engine entries | [`PackEngine`], [`GencanPack`], [`CbmcGrow`], [`LatticeGrow`], [`LogLevel`] |
//! | Run lifecycle | [`Pipeline`], [`StageFactory`], [`EngineSetup`], [`State`], [`IntraResidual`] |
//! | Stage combinators (Rust-only) | [`Until`], [`OnViolation`], [`Invariant`], [`Layers`], [`Violation`], [`RestraintsSatisfied`] |
//! | Shared settings + space | [`PackSettings`] |
//! | Target  | [`Target`], [`CenteringMode`], [`Angle`], [`Axis`], [`Placement`] |
//! | Run state | [`PackSystem`], [`RigidView`] |
//! | Stage seam (Rust-only) | [`Stage`], [`Requires`], [`Guarantees`], [`StageOutcome`], [`Budget`], [`PackState`], [`Placed`] |
//! | Per-atom restraints | [`AtomRestraint`], [`RegionRestraint`] (over [`molrs::core::Region`]), [`CellRestraint`] |
//! | Group restraints | [`Restraint`], [`GroupCtx`], [`GaussianPlane`], [`GaussianPoint`], [`ExponentialPlane`], [`ExponentialPoint`], [`TabulatedPlane`], [`TabulatedPoint`], [`SelfSeparation`] |
//! | Callback trait + built-ins | [`Callback`], [`LammpsLogCallback`], [`ProgressCallback`], [`EarlyStopCallback`], `XyzTrajectoryCallback` (feature `io`), [`StepReport`], [`StageProgress`], [`PhaseProgress`], [`PhaseReport`] |
//! | Objective | [`Objective`], [`EvalMode`], [`EvalOutput`] |
//! | In-loop optimizer | [`OptimizeSelect`], [`OptimizeMode`], [`GencanPack::with_optimizer`], [`TorsionMcOptimizer`], over molrs's [`molrs::optimize::Optimizer`] trait |
//! | Errors | [`PackError`] |
//!
//! ## Feature flags
//!
//! - `rayon` — opt into the parallel evaluator (also forwards to `molrs`'s
//!   `rayon`).
//! - `io` — pull in molrs's `io` module so `script::Script::build` reads the
//!   template files through the molrs reader of each one's
//!   [`script::StructureFormat`] (PDB, XYZ, SDF/MOL, LAMMPS, …) and hands back a
//!   `script::BuildResult`, and so `XyzTrajectoryCallback` can write its trajectory
//!   through molrs's extended XYZ writer. PyO3 / WASM / embedding hosts that
//!   bring their own loader leave this off and use [`script::Script::lower`]
//!   with [`script::StructurePlan::apply`] instead.
//! - `cli` — build the `molpack` binary (pulls in `clap` and implies `io`).
//!
//! A force-field optimizer bound through [`GencanPack::with_optimizer`]
//! (molrs's `Lbfgs` over a `Potential`) needs molrs's `ff` feature, which the
//! caller turns on in its own molrs dependency; molpack has no `ff` feature.
//!
//! Precision is fixed at `f64` via `molrs::op::F`.

mod assemble;
mod callback;
pub mod context;
mod error;
mod euler;
mod eval;
pub mod grow;
mod invariant;
mod objective;
mod optimizer;
mod outcome;
mod pack;
mod pack_space;
mod pipeline;
mod random;
mod restraint;
pub mod script;
mod settings;
mod stage;
mod state;
mod target;
mod template;
#[cfg(test)]
mod test_fixtures;

#[cfg(feature = "io")]
pub use callback::XyzTrajectoryCallback;
pub use callback::{
    Callback, EarlyStopCallback, LammpsLogCallback, LogLevel, PhaseProgress, PhaseReport,
    ProgressCallback, StageProgress, StepReport,
};
pub use context::pack_state::{PackState, Placed};
pub use context::pack_system::PackSystem;
pub use context::rigid_view::RigidView;
pub use error::PackError;
pub use eval::{EvalMode, EvalOutput};
pub use grow::cbmc_grow::CbmcGrow;
pub use grow::lattice::LatticeGrow;
pub use invariant::{Invariant, Layers, RestraintsSatisfied, Violation};
pub use objective::Objective;
pub use optimizer::{OptimizeMode, OptimizeSelect, TorsionMcOptimizer};
pub use pack::GencanPack;
pub use pipeline::{EngineSetup, OnViolation, PackEngine, Pipeline, StageFactory, Until};
pub use restraint::{
    AtomRestraint, CellRestraint, ExponentialPlane, ExponentialPoint, GaussianPlane, GaussianPoint,
    GroupCtx, RegionRestraint, Restraint, SelfSeparation, TabulatedPlane, TabulatedPoint,
};
pub use settings::PackSettings;
pub use stage::{Budget, Guarantees, Requires, Stage, StageOutcome};
pub use state::{IntraResidual, State};
pub use target::{Angle, Axis, CenteringMode, Placement, Target};

// ────────────────────────────────────────────────────────────────────────────
// Documentation modules (rustdoc-only; no runtime items).
// Content lives in `docs/*.md`, loaded via `include_str!` so each markdown
// file can be edited independently while rustdoc renders the whole chapter.
// ────────────────────────────────────────────────────────────────────────────

#[doc = include_str!("../docs/getting_started.md")]
pub mod getting_started {}

#[doc = include_str!("../docs/concepts.md")]
pub mod concepts {}

#[doc = include_str!("../docs/architecture.md")]
pub mod architecture {}

#[doc = include_str!("../docs/extending.md")]
pub mod extending {}
