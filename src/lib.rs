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
//!   the molrs `Region` lift, `Handler`, `Objective`, `Target`, `PackEngine`,
//!   `PackContext`; the scope equivalence law; the two-scale contract;
//!   the direction-3 extension pattern.
//! - [`architecture`] — module map, dependency graph, core-type
//!   relationships, full `pack()` lifecycle diagram, hot-path
//!   `evaluate()` walkthrough, invariants, design decisions.
//! - [`extending`] — tutorials for writing your own `AtomRestraint` /
//!   molrs `Region` / `Handler` and for binding an in-loop optimizer; testing
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
//! use molpack::{GenCanPack, PackEngine, RegionRestraint, Target};
//! use molrs::spatial::region::Cuboid;
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
//! let result = GenCanPack::new()
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
//! | Category | Items |
//! |---|---|
//! | Engine entries | [`PackEngine`], [`GenCanPack`], [`CbmcGrow`], [`LatticeGrow`], [`LogLevel`] |
//! | Run lifecycle | [`Pipeline`], [`StageFactory`], [`pipeline::EngineSetup`], [`State`], [`IntraResidual`] |
//! | Stage combinators (Rust-only) | [`Until`], [`OnViolation`], [`Invariant`], [`Layers`], [`Violation`], [`RestraintsSatisfied`] |
//! | Shared settings + space | [`PackSettings`] (`entry`) |
//! | Target  | [`Target`], [`CenteringMode`] |
//! | Rigid placement vector | [`RigidView`] |
//! | Stage seam (Rust-only) | [`Stage`], [`Requires`], [`Guarantees`], [`StageOutcome`], [`Budget`], [`PackState`], [`Placed`] |
//! | AtomRestraint trait + the region lift | [`AtomRestraint`], [`RegionRestraint`] (over [`molrs::spatial::region::Region`]), [`CellRestraint`] |
//! | Handler trait + built-ins | [`Handler`], [`LammpsLogHandler`], [`ProgressHandler`], [`EarlyStopHandler`], [`XYZHandler`], [`StepInfo`], [`handler::StageInfo`], [`PhaseInfo`], [`PhaseReport`] |
//! | In-loop optimizer | [`OptimizeSelect`], [`OptimizeMode`], [`GenCanPack::with_optimizer`], [`TorsionMcOptimizer`], and molrs's [`Optimizer`] trait |
//! | Errors | [`PackError`] |
//!
//! ## Feature flags
//!
//! - `rayon` — opt into the parallel evaluator (also forwards to `molrs`'s
//!   `rayon`).
//! - `io` — pull in molrs's `io` module so `script::Script::build` can read
//!   Protein Data Bank (PDB) / structure-data-file (SDF) / XYZ / LAMMPS files
//!   directly and hand back a `script::BuildResult`. PyO3 / WASM / embedding
//!   hosts that bring their own loader leave this off and use
//!   [`script::Script::lower`] with [`script::StructurePlan::apply`] instead.
//! - `cli` — build the `molpack` binary and its integration tests (pulls in
//!   `clap` and implies `io`).
//! - `ff` — forward molrs's `ff` module (force fields and the optimizers
//!   built on them) for callers who bind one through
//!   [`GenCanPack::with_optimizer`]. molpack compiles nothing extra under it:
//!   the in-loop optimizer seam itself is always on.
//!
//! Precision is fixed at `f64` via `molrs::types::F`.

pub mod assemble;
pub mod constraints;
pub mod context;
pub mod entry;
pub mod error;
pub mod euler;
pub mod gencan;
pub mod grow;
pub mod handler;
pub mod initial;
pub mod invariant;
pub mod movebad;
mod numerics;
pub mod objective;
pub mod optimizer;
pub mod pipeline;
mod random;
pub mod restraint;
pub mod script;
pub mod stage;
pub mod target;
mod template;
#[cfg(test)]
mod testutil;

pub use context::{PackContext, PackState, Placed, RigidView};
pub use entry::IntraResidual;
pub use entry::PackSettings;
pub use entry::State;
pub use error::PackError;
pub use gencan::entry::GenCanPack;
pub use grow::entry::CbmcGrow;
pub use grow::lattice::{LatticeConfig, LatticeGrow};
pub use handler::{
    EarlyStopHandler, Handler, LammpsLogHandler, LogLevel, PhaseInfo, PhaseReport, ProgressHandler,
    StepInfo, XYZHandler,
};
pub use invariant::{Invariant, Layers, RestraintsSatisfied, Violation};
pub use molrs::BondDistanceWeights;
pub use molrs::Element;
pub use molrs::types::F;
pub use pipeline::combinators::{OnViolation, Until};
pub use pipeline::{PackEngine, Pipeline, StageFactory};
// The in-loop optimizer seam. The trait and its report (molrs core) appear in
// molpack's own signatures (`with_optimizer`); concrete force-field optimizers
// (`LBFGS`, `Potential`) stay at their molrs home.
pub use molrs::optimize::{OptReport, Optimizer};
pub use optimizer::{OptimizeMode, OptimizeSelect, TorsionMcOptimizer};
pub use restraint::{AtomRestraint, CellRestraint, RegionRestraint};
pub use stage::{Budget, Guarantees, Requires, Stage, StageOutcome};
pub use target::{Angle, Axis, CenteringMode, Placement, Target};

// Custom-objective extension surface. An engine run drives a `dyn Objective`
// through GENCAN; downstream code that implements a bespoke objective (or wants
// to evaluate the packing energy/gradient directly) names these at the crate
// root rather than reaching into the `objective` / `constraints` modules.
pub use constraints::{Constraints, EvalMode, EvalOutput};
pub use objective::Objective;

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

// ────────────────────────────────────────────────────────────────────────────
// Prelude — bulk re-export of the vocabulary a typical packing script needs.
// ────────────────────────────────────────────────────────────────────────────

/// Bulk re-export of the items a typical packing script needs.
///
/// ```no_run
/// use std::sync::Arc;
/// use molpack::prelude::*;
/// use molrs::spatial::region::Cuboid;
/// use ndarray::array;
///
/// let cube = Cuboid::new(array![0.0, 0.0, 0.0], array![10.0, 10.0, 10.0]);
/// let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 10)
///     .with_restraint(RegionRestraint(Arc::new(cube)));
/// let result = GenCanPack::new().run(&[target], 100)?;
/// # Ok::<(), molpack::PackError>(())
/// ```
///
/// The crate root still re-exports everything for direct `use molpack::T`
/// access; the prelude exists to avoid a 20-line `use` block at the top of
/// every example.
pub mod prelude {
    pub use crate::{
        // Target + centering + angle / axis / placement
        Angle,
        // AtomRestraint trait + the molrs region lift + the cell
        AtomRestraint,
        Axis,
        CbmcGrow,
        CellRestraint,
        CenteringMode,
        // Handlers
        EarlyStopHandler,
        GenCanPack,
        Handler,
        // Core builder + result + error
        LammpsLogHandler,
        LogLevel,
        // In-loop optimizer seam
        OptimizeMode,
        OptimizeSelect,
        Optimizer,
        PackEngine,
        PackError,
        PhaseInfo,
        PhaseReport,
        Placement,
        ProgressHandler,
        RegionRestraint,
        State,
        StepInfo,
        Target,
        TorsionMcOptimizer,
        XYZHandler,
    };
}
