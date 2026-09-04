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
//! for shared data structures (always-on `core`) and, behind feature flags, its
//! file I/O (`io`) and force-field (`ff`) modules.
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
//!   `Region`, `Handler`, `Objective`, `Target`, `PackEngine`,
//!   `PackContext`; the scope equivalence law; the two-scale contract;
//!   the direction-3 extension pattern.
//! - [`architecture`] — module map, dependency graph, core-type
//!   relationships, full `pack()` lifecycle diagram, hot-path
//!   `evaluate()` walkthrough, invariants, design decisions.
//! - [`extending`] — tutorials for writing your own `AtomRestraint` /
//!   `Region` / `Handler` and for binding an in-loop optimizer; testing +
//!   benchmarking discipline; common pitfalls; contributing flow.
//!
//! Reference material (not rustdoc):
//!
//! - [Packmol parity](https://github.com/MolCrafts/molpack/blob/master/docs/packmol_parity.md)
//!   — kind-number ↔ Rust struct mapping with Fortran pointers.
//!
//! ## Quick example
//!
//! ```rust,no_run
//! use molpack::{GenCanPack, InsideBoxRestraint, PackEngine, Target};
//!
//! let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
//! let radii     = [1.52, 1.20, 1.20];
//!
//! let target = Target::from_coords(&positions, &radii, 100)
//!     .with_name("water")
//!     .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0, 40.0, 40.0], [false; 3]));
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
//! | Run lifecycle | [`Pipeline`], [`StageFactory`], [`pipeline::EngineSetup`], [`PackResult`], [`IntraResidual`] |
//! | Stage combinators (Rust-only) | [`Until`], [`OnViolation`], [`Invariant`], [`Layers`], [`Violation`], [`RestraintsSatisfied`] |
//! | Shared settings + space | [`PackSettings`] (`entry`) |
//! | Target  | [`Target`], [`CenteringMode`] |
//! | Rigid placement vector | [`RigidView`] |
//! | Stage seam (Rust-only) | [`Stage`], [`Requires`], [`Guarantees`], [`StageOutcome`], [`Budget`], [`PackState`], [`Placed`] |
//! | AtomRestraint trait + 14 concrete structs | [`AtomRestraint`] + `InsideBox` / `InsideCube` / `InsideSphere` / `InsideEllipsoid` / `InsideCylinder` / `Outside*` variants / `AbovePlane` / `BelowPlane` / `AboveGaussian` / `BelowGaussian` — each suffixed `…AtomRestraint` |
//! | Region trait + combinators + lift | [`Region`], [`RegionExt`], [`And`], [`Or`], [`Not`], [`RegionRestraint`], [`InsideBoxRegion`], [`InsideCellRegion`], [`InsideSphereRegion`], [`OutsideSphereRegion`], [`Aabb`] |
//! | Handler trait + built-ins | [`Handler`], [`NullHandler`], [`LammpsLogHandler`], [`ProgressHandler`], [`EarlyStopHandler`], [`XYZHandler`], [`StepInfo`], [`handler::StageInfo`], [`PhaseInfo`], [`PhaseReport`] |
//! | In-loop optimizer (feature `ff`) | `OptimizeSelect`, `GenCanPack::with_optimizer`, `TorsionMcOptimizer`, and molrs's `Optimizer` trait |
//! | Errors | [`PackError`] |
//! | Validation | [`validate_from_targets`], [`ValidationReport`], [`ViolationMetrics`] |
//! | Examples harness (feature `io`) | `ExampleCase`, `build_targets`, `example_dir_from_manifest`, `render_inp_script` |
//!
//! The last two rows name items that exist only when their Cargo feature is
//! enabled. A default-feature documentation build cannot resolve a link to
//! something it did not compile, so those names are written in plain code font
//! rather than as cross-references; build with `--features ff,io` to see them
//! in this crate's rustdoc.
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
//! - `ff` — pull in molrs's `ff` module (typifiers for the Merck Molecular
//!   Force Field, MMFF94 / MMFF94s, plus the limited-memory
//!   Broyden–Fletcher–Goldfarb–Shanno minimizer, L-BFGS) and enable the in-loop
//!   optimizer bindings `GenCanPack::with_optimizer` + `OptimizeSelect`.
//!
//! Precision is fixed at `f64` via `molrs::types::F`.

pub mod assemble;
#[cfg(feature = "io")]
pub mod cases;
pub mod constraints;
pub mod context;
pub mod entry;
pub mod error;
pub mod euler;
pub mod frame;
pub mod gencan;
pub mod grow;
pub mod handler;
pub mod initial;
pub mod invariant;
pub mod movebad;
mod numerics;
pub mod objective;
#[cfg(feature = "ff")]
pub mod optimizer;
pub mod pipeline;
mod random;
pub mod region;
pub mod restraint;
pub mod script;
pub mod stage;
pub mod target;
mod template;
pub mod validation;

#[cfg(feature = "io")]
pub use cases::{ExampleCase, build_targets, example_dir_from_manifest, render_inp_script};
pub use context::{PackContext, PackState, Placed, RigidView};
pub use entry::IntraResidual;
pub use entry::PackResult;
pub use entry::PackSettings;
pub use error::PackError;
pub use frame::{compute_mol_ids, context_to_frame, finalize_frame, frame_to_coords};
pub use gencan::entry::GenCanPack;
pub use grow::entry::CbmcGrow;
pub use grow::lattice::{LatticeConfig, LatticeGrow};
pub use handler::{
    EarlyStopHandler, Handler, LammpsLogHandler, LogLevel, NullHandler, PhaseInfo, PhaseReport,
    ProgressHandler, StepInfo, XYZHandler,
};
pub use invariant::{Invariant, Layers, RestraintsSatisfied, Violation};
pub use molrs::Element;
pub use molrs::types::F;
pub use pipeline::combinators::{OnViolation, Until};
pub use pipeline::{PackEngine, Pipeline, StageFactory};
pub use region::{
    Aabb, And, InsideBoxRegion, InsideCellRegion, InsideSphereRegion, Not, Or, OutsideSphereRegion,
    Region, RegionExt, RegionRestraint,
};
// In-loop optimizers require molrs `ff` (Optimizer trait + Potential).
#[cfg(feature = "ff")]
pub use molrs::ff::potential::Potential;
#[cfg(feature = "ff")]
pub use molrs::optimize::{LBFGS, OptReport, Optimizer};
#[cfg(feature = "ff")]
pub use optimizer::{OptimizeMode, OptimizeSelect, TorsionMcOptimizer};
pub use restraint::{
    AboveGaussianRestraint, AbovePlaneRestraint, AtomRestraint, BelowGaussianRestraint,
    BelowPlaneRestraint, InsideBoxRestraint, InsideCubeRestraint, InsideCylinderRestraint,
    InsideEllipsoidRestraint, InsideSphereRestraint, OutsideBoxRestraint, OutsideCubeRestraint,
    OutsideCylinderRestraint, OutsideEllipsoidRestraint, OutsideSphereRestraint,
};
pub use stage::{Budget, Guarantees, Requires, Stage, StageOutcome};
pub use target::{Angle, Axis, CenteringMode, Placement, Target};
pub use validation::{ValidationReport, ViolationMetrics, validate_from_targets};

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
/// use molpack::prelude::*;
///
/// let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 10)
///     .with_restraint(InsideBoxRestraint::new([0.0; 3], [10.0; 3], [false; 3]));
/// let result = GenCanPack::new().run(&[target], 100)?;
/// # Ok::<(), molpack::PackError>(())
/// ```
///
/// The crate root still re-exports everything for direct `use molpack::T`
/// access; the prelude exists to avoid a 20-line `use` block at the top of
/// every example.
pub mod prelude {
    pub use crate::{
        // Region + combinators + lift
        Aabb,
        // AtomRestraint trait + 14 concrete impls
        AboveGaussianRestraint,
        AbovePlaneRestraint,
        And,
        // Target + centering + angle / axis / placement
        Angle,
        AtomRestraint,
        Axis,
        BelowGaussianRestraint,
        BelowPlaneRestraint,
        CbmcGrow,
        CenteringMode,
        // Handlers
        EarlyStopHandler,
        GenCanPack,
        Handler,
        InsideBoxRegion,
        InsideBoxRestraint,
        InsideCellRegion,
        InsideCubeRestraint,
        InsideCylinderRestraint,
        InsideEllipsoidRestraint,
        InsideSphereRegion,
        InsideSphereRestraint,
        // Core builder + result + error
        LammpsLogHandler,
        LogLevel,
        Not,
        NullHandler,
        Or,
        OutsideBoxRestraint,
        OutsideCubeRestraint,
        OutsideCylinderRestraint,
        OutsideEllipsoidRestraint,
        OutsideSphereRegion,
        OutsideSphereRestraint,
        PackEngine,
        PackError,
        PackResult,
        PhaseInfo,
        PhaseReport,
        Placement,
        ProgressHandler,
        Region,
        RegionExt,
        RegionRestraint,
        StepInfo,
        Target,
        XYZHandler,
    };
    #[cfg(feature = "ff")]
    pub use crate::{LBFGS, OptimizeMode, OptimizeSelect, Optimizer, TorsionMcOptimizer};
}
