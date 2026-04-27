//! # molpack
//!
//! Packmol-grade molecular packing in pure Rust. Produces a non-overlapping
//! arrangement of N molecule types with copy counts and geometric restraints,
//! using a faithful port of Packmol's GENCAN-driven three-phase algorithm
//! (Martínez et al. 2009). Correctness is checked against Packmol's reference
//! output for five canonical workloads.
//!
//! This crate was split out of the molrs workspace in 2026 and is now
//! maintained independently. It still depends on `molcrafts-molrs-core` and
//! `molcrafts-molrs-io` for shared data structures and file I/O.
//!
//! ## Documentation map
//!
//! This crate is documented in four dedicated modules; start with
//! [`getting_started`] if you are new.
//!
//! - [`getting_started`] — install, hello-world packing, the three
//!   restraint scopes, handlers, relaxers, PBC, running the canonical
//!   examples.
//! - [`concepts`] — every abstraction defined in one place: `Restraint`,
//!   `Region`, `Relaxer`, `Handler`, `Objective`, `Target`, `Molpack`,
//!   `PackContext`; the scope equivalence law; the two-scale contract;
//!   the direction-3 extension pattern.
//! - [`architecture`] — module map, dependency graph, core-type
//!   relationships, full `pack()` lifecycle diagram, hot-path
//!   `evaluate()` walkthrough, invariants, design decisions.
//! - [`extending`] — tutorials for writing your own `Restraint` /
//!   `Region` / `Handler` / `Relaxer`; testing + benchmarking
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
//! use molpack::{InsideBoxRestraint, Molpack, Target};
//!
//! let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
//! let radii     = [1.52, 1.20, 1.20];
//!
//! let target = Target::from_coords(&positions, &radii, 100)
//!     .with_name("water")
//!     .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0, 40.0, 40.0], [false; 3]));
//!
//! let result = Molpack::new()
//!     .with_tolerance(2.0)
//!     .with_precision(0.01)
//!     .with_seed(42)
//!     .pack(&[target], 200)?;
//!
//! println!("converged = {}, natoms = {}", result.converged, result.natoms());
//! # Ok::<(), molpack::PackError>(())
//! ```
//!
//! ## Public surface at a glance
//!
//! | Category | Items |
//! |---|---|
//! | Builder | [`Molpack`], [`PackResult`] |
//! | Target  | [`Target`], [`CenteringMode`] |
//! | Restraint trait + 14 concrete structs | [`Restraint`] + `InsideBox` / `InsideCube` / `InsideSphere` / `InsideEllipsoid` / `InsideCylinder` / `Outside*` variants / `AbovePlane` / `BelowPlane` / `AboveGaussian` / `BelowGaussian` — each suffixed `…Restraint` |
//! | Region trait + combinators + lift | [`Region`], [`RegionExt`], [`And`], [`Or`], [`Not`], [`RegionRestraint`], [`InsideBoxRegion`], [`InsideSphereRegion`], [`OutsideSphereRegion`], [`Aabb`] |
//! | Handler trait + built-ins | [`Handler`], [`NullHandler`], [`ProgressHandler`], [`EarlyStopHandler`], [`XYZHandler`], [`StepInfo`], [`PhaseInfo`], [`PhaseReport`] |
//! | Relaxer trait + built-in | [`Relaxer`], [`RelaxerRunner`], [`TorsionMcRelaxer`] (aliased to `TorsionMcHook`) |
//! | Errors | [`PackError`] |
//! | Validation | [`validate_from_targets`], [`ValidationReport`], [`ViolationMetrics`] |
//! | Examples harness | [`ExampleCase`], [`build_targets`], [`example_dir_from_manifest`], [`render_packmol_input`] |
//!
//! ## Feature flags
//!
//! - `rayon` — opt into the parallel evaluator (also forwards to `molrs-core`).
//! - `io` — pull in molrs-io so [`script::Script::build`] can read PDB / SDF /
//!   XYZ / LAMMPS files directly. PyO3 / WASM / embedding hosts that bring
//!   their own loader leave this off and use [`script::Script::lower`] with
//!   [`script::StructurePlan::apply`] instead.
//! - `cli` — build the `molpack` binary and its integration tests (pulls in
//!   `clap` and implies `io`).
//!
//! Precision is fixed at `f64` via `molrs::types::F`.

pub mod api;
#[cfg(feature = "io")]
pub mod cases;
pub mod cell;
pub mod constraints;
pub mod context;
pub mod error;
pub mod euler;
pub mod frame;
pub mod gencan;
pub mod handler;
pub mod initial;
pub mod movebad;
mod numerics;
pub mod objective;
pub mod packer;
mod random;
pub mod region;
pub mod relaxer;
pub mod restraint;
pub mod script;
pub mod target;
pub mod validation;

#[cfg(feature = "io")]
pub use cases::{ExampleCase, build_targets, example_dir_from_manifest, render_packmol_input};
pub use context::PackContext;
pub use error::PackError;
pub use frame::{compute_mol_ids, context_to_frame, finalize_frame, frame_to_coords};
pub use handler::{
    EarlyStopHandler, Handler, NullHandler, PhaseInfo, PhaseReport, ProgressHandler, StepInfo,
    XYZHandler,
};
pub use molrs::Element;
pub use molrs::types::F;
pub use packer::{Molpack, PackResult};
pub use region::{
    Aabb, And, InsideBoxRegion, InsideSphereRegion, Not, Or, OutsideSphereRegion, Region,
    RegionExt, RegionRestraint,
};
#[allow(deprecated)]
pub use region::{BBox, FromRegion};
pub use relaxer::{Relaxer, RelaxerRunner, TorsionMcRelaxer};
pub use restraint::{
    AboveGaussianRestraint, AbovePlaneRestraint, BelowGaussianRestraint, BelowPlaneRestraint,
    InsideBoxRestraint, InsideCubeRestraint, InsideCylinderRestraint, InsideEllipsoidRestraint,
    InsideSphereRestraint, OutsideBoxRestraint, OutsideCubeRestraint, OutsideCylinderRestraint,
    OutsideEllipsoidRestraint, OutsideSphereRestraint, Restraint,
};
#[allow(deprecated)]
pub use target::FixedPlacement;
pub use target::{Angle, Axis, CenteringMode, Placement, Target};
pub use validation::{ValidationReport, ViolationMetrics, validate_from_targets};

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
/// let result = Molpack::new().pack(&[target], 100)?;
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
        // Restraint trait + 14 concrete impls
        AboveGaussianRestraint,
        AbovePlaneRestraint,
        And,
        // Target + centering + angle / axis / placement
        Angle,
        Axis,
        BelowGaussianRestraint,
        BelowPlaneRestraint,
        CenteringMode,
        // Handlers
        EarlyStopHandler,
        Handler,
        InsideBoxRegion,
        InsideBoxRestraint,
        InsideCubeRestraint,
        InsideCylinderRestraint,
        InsideEllipsoidRestraint,
        InsideSphereRegion,
        InsideSphereRestraint,
        // Core builder + result + error
        Molpack,
        Not,
        NullHandler,
        Or,
        OutsideBoxRestraint,
        OutsideCubeRestraint,
        OutsideCylinderRestraint,
        OutsideEllipsoidRestraint,
        OutsideSphereRegion,
        OutsideSphereRestraint,
        PackError,
        PackResult,
        PhaseInfo,
        PhaseReport,
        Placement,
        ProgressHandler,
        Region,
        RegionExt,
        RegionRestraint,
        // Relaxer
        Relaxer,
        RelaxerRunner,
        Restraint,
        StepInfo,
        Target,
        TorsionMcRelaxer,
        XYZHandler,
    };
}
