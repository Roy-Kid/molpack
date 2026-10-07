//! PyO3 bindings for molpack.
//!
//! Python API:
//!
//! | Python class     | Rust wrapper         | Purpose                            |
//! |------------------|----------------------|------------------------------------|
//! | `Target`         | [`PyTarget`]         | Molecule specification for packing |
//! | `GencanPack`     | [`PyGencanPack`]     | Rigid-body GENCAN packing entry    |
//! | `CbmcGrow`       | [`PyCbmcGrow`]       | Chain-growth entry                 |
//! | `LatticeGrow`    | [`PyLatticeGrow`]    | Lattice-growth entry               |
//! | `Pipeline`       | [`PyPipeline`]       | Entries composed as stages         |
//! | `State`     | [`PyState`]     | Frame + diagnostics from `run()`   |
//! | `IntraResidual`  | [`PyIntraResidual`]  | Nested scored/exempted intra mins  |
//! | `StepReport`       | [`PyStepReport`]       | Read-only snapshot for callbacks    |
//! | `StageProgress`      | [`PyStageProgress`]      | Which stage a callback came from   |
//! | `StepContext`    | [`PyStepContext`]    | Callback-scoped live-context guard |
//!
//! Geometric restraints are molrs region objects (`molrs.core.Sphere`, `Cuboid`,
//! `Parallelepiped`, `HalfSpace`, `Cylinder`, `Ellipsoid`, `Polyhedron`,
//! `SphereUnion`, or a `&` / `|` / `~` composition), resolved through their
//! `molrs.RegionRef/<line>` capsule and lifted by `RegionRestraint` — this
//! wheel defines no geometry class. Custom Python restraints are attached by
//! passing any object with callable `f(x, scale, scale2)` and
//! `fg(x, scale, scale2)` methods to `Target.with_restraint` — no dedicated
//! class needed.
//!
//! Custom Python progress callbacks are registered via the entries'
//! `with_callback(obj)`; see the [`callback`] module for the method contract.

use pyo3::prelude::*;

mod molrs_capsule;

mod errors;
use errors::register_errors;

mod restraint;
use restraint::{
    PyExponentialPlane, PyExponentialPoint, PyGaussianPlane, PyGaussianPoint, PySelfSeparation,
    PyTabulatedPlane, PyTabulatedPoint,
};

mod callback;
use callback::{PyStageProgress, PyStepContext, PyStepReport};

mod grow;

mod target;
use target::{PyAngle, PyAxis, PyCenteringMode, PyTarget};

mod packing_methods;
use packing_methods::{PyCbmcGrow, PyGencanPack, PyLatticeGrow, PyPipeline};
mod state;
use state::{PyIntraResidual, PyState};

mod parallel;
use parallel::{init_thread_pool, num_threads, rayon_enabled};

mod script;
use script::{PyScriptJob, load_script};

#[pymodule]
fn molpack(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // ABI handshake first: this extension exchanges molrs_ffi handle capsules
    // with the installed molcrafts-molrs wheel, so both must embed the same
    // molrs minor line (minor-line = ABI version). A mismatch must be a clear
    // ImportError here, not a capsule ValueError (or worse) mid-run.
    molrs_capsule::check_abi(m.py())?;

    m.add_class::<PyAngle>()?;
    m.add_class::<PyAxis>()?;
    m.add_class::<PyCenteringMode>()?;

    m.add_class::<PyGaussianPlane>()?;
    m.add_class::<PyGaussianPoint>()?;
    m.add_class::<PyExponentialPlane>()?;
    m.add_class::<PyExponentialPoint>()?;
    m.add_class::<PyTabulatedPlane>()?;
    m.add_class::<PyTabulatedPoint>()?;
    m.add_class::<PySelfSeparation>()?;

    m.add_class::<grow::PyTorsionPrior>()?;
    m.add_class::<grow::PyAnglePrior>()?;

    m.add_class::<PyTarget>()?;
    m.add_class::<PyGencanPack>()?;
    m.add_class::<PyCbmcGrow>()?;
    m.add_class::<PyLatticeGrow>()?;
    m.add_class::<PyPipeline>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyIntraResidual>()?;
    m.add_class::<PyStepReport>()?;
    m.add_class::<PyStageProgress>()?;
    m.add_class::<PyStepContext>()?;

    m.add_class::<PyScriptJob>()?;
    m.add_function(wrap_pyfunction!(load_script, m)?)?;

    m.add_function(wrap_pyfunction!(rayon_enabled, m)?)?;
    m.add_function(wrap_pyfunction!(num_threads, m)?)?;
    m.add_function(wrap_pyfunction!(init_thread_pool, m)?)?;

    register_errors(m.py(), m)?;

    Ok(())
}
