//! PyO3 bindings for molpack.
//!
//! Python API:
//!
//! | Python class     | Rust wrapper         | Purpose                            |
//! |------------------|----------------------|------------------------------------|
//! | `Target`         | [`PyTarget`]         | Molecule specification for packing |
//! | `GenCanPack`     | [`PyGenCanPack`]     | Rigid-body GENCAN packing entry    |
//! | `CbmcGrow`       | [`PyCbmcGrow`]       | Chain-growth entry                 |
//! | `LatticeGrow`    | [`PyLatticeGrow`]    | Lattice-growth entry               |
//! | `Pipeline`       | [`PyPipeline`]       | Entries composed as stages         |
//! | `State`     | [`PyState`]     | Frame + diagnostics from `run()`   |
//! | `IntraResidual`  | [`PyIntraResidual`]  | Nested scored/exempted intra mins  |
//! | `StepInfo`       | [`PyStepInfo`]       | Read-only snapshot for handlers    |
//! | `StageInfo`      | [`PyStageInfo`]      | Which stage a callback came from   |
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
//! Custom Python progress handlers are registered via the entries'
//! `with_handler(obj)`; see the [`handler`] module for the method contract.

use pyo3::prelude::*;

mod interop;

mod errors;
use errors::register_errors;

mod types;
use types::{PyAngle, PyAxis, PyCenteringMode};

mod restraint;
use restraint::{
    PyExponentialPlane, PyExponentialPoint, PyGaussianPlane, PyGaussianPoint, PySelfSeparation,
    PyTabulatedPlane, PyTabulatedPoint,
};

mod handler;
use handler::{PyStageInfo, PyStepContext, PyStepInfo};

mod grow;

mod target;
use target::PyTarget;

mod entry;
use entry::{PyCbmcGrow, PyGenCanPack, PyLatticeGrow, PyPipeline};
mod result;
use result::{PyIntraResidual, PyState};

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
    interop::check_abi(m.py())?;

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
    m.add_class::<PyGenCanPack>()?;
    m.add_class::<PyCbmcGrow>()?;
    m.add_class::<PyLatticeGrow>()?;
    m.add_class::<PyPipeline>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyIntraResidual>()?;
    m.add_class::<PyStepInfo>()?;
    m.add_class::<PyStageInfo>()?;
    m.add_class::<PyStepContext>()?;

    m.add_class::<PyScriptJob>()?;
    m.add_function(wrap_pyfunction!(load_script, m)?)?;

    m.add_function(wrap_pyfunction!(rayon_enabled, m)?)?;
    m.add_function(wrap_pyfunction!(num_threads, m)?)?;
    m.add_function(wrap_pyfunction!(init_thread_pool, m)?)?;

    register_errors(m.py(), m)?;

    Ok(())
}
