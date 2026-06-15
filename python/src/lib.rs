//! PyO3 bindings for molpack.
//!
//! Python API:
//!
//! | Python class     | Rust wrapper         | Purpose                            |
//! |------------------|----------------------|------------------------------------|
//! | `Target`         | [`PyTarget`]         | Molecule specification for packing |
//! | `Molpack`        | [`PyPacker`]         | Molecular packer (Packmol port)    |
//! | `PackResult`     | [`PyPackResult`]     | Diagnostics from `pack_with_report` |
//! | `StepInfo`       | [`PyStepInfo`]       | Read-only snapshot for handlers    |
//! | `InsideBox`      | [`PyInsideBox`]      | Box restraint                      |
//! | `InsideSphere`   | [`PyInsideSphere`]   | Sphere restraint (inside)          |
//! | `OutsideSphere`  | [`PyOutsideSphere`]  | Sphere restraint (outside)         |
//! | `AbovePlane`     | [`PyAbovePlane`]     | Half-space restraint               |
//! | `BelowPlane`     | [`PyBelowPlane`]     | Half-space restraint               |
//!
//! Custom Python restraints are attached by passing any object with
//! callable `f(x, scale, scale2)` and `fg(x, scale, scale2)` methods to
//! `Target.with_restraint` — no dedicated class needed.
//!
//! Custom Python progress handlers are registered via
//! `Molpack.add_handler(obj)`; see the [`handler`] module for the
//! method contract.

use pyo3::prelude::*;

mod interop;

mod helpers;
use helpers::register_errors;

mod types;
use types::{PyAngle, PyAxis, PyCenteringMode};

mod constraint;
use constraint::{
    PyAbovePlaneRestraint, PyBelowPlaneRestraint, PyExponentialPlane, PyExponentialPoint,
    PyGaussianPlane, PyGaussianPoint, PyInsideBoxRestraint, PyInsideSphereRestraint,
    PyOutsideSphereRestraint, PyTabulatedPlane, PyTabulatedPoint,
};

mod handler;
use handler::PyStepInfo;

mod target;
use target::PyTarget;

mod relaxer;

mod packer;
use packer::{PyPackResult, PyPacker};

mod parallel;
use parallel::{init_thread_pool, num_threads, rayon_enabled};

mod script;
use script::{PyScriptJob, load_script};

#[pymodule]
fn molpack(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAngle>()?;
    m.add_class::<PyAxis>()?;
    m.add_class::<PyCenteringMode>()?;

    m.add_class::<PyInsideBoxRestraint>()?;
    m.add_class::<PyInsideSphereRestraint>()?;
    m.add_class::<PyOutsideSphereRestraint>()?;
    m.add_class::<PyAbovePlaneRestraint>()?;
    m.add_class::<PyBelowPlaneRestraint>()?;
    m.add_class::<PyGaussianPlane>()?;
    m.add_class::<PyGaussianPoint>()?;
    m.add_class::<PyExponentialPlane>()?;
    m.add_class::<PyExponentialPoint>()?;
    m.add_class::<PyTabulatedPlane>()?;
    m.add_class::<PyTabulatedPoint>()?;

    m.add_class::<PyTarget>()?;
    m.add_class::<PyPacker>()?;
    m.add_class::<PyPackResult>()?;
    m.add_class::<PyStepInfo>()?;

    m.add_class::<relaxer::PyTorsionMcRelaxer>()?;
    #[cfg(feature = "ff")]
    m.add_class::<relaxer::PyLBFGSRelaxer>()?;

    m.add_class::<PyScriptJob>()?;
    m.add_function(wrap_pyfunction!(load_script, m)?)?;

    m.add_function(wrap_pyfunction!(rayon_enabled, m)?)?;
    m.add_function(wrap_pyfunction!(num_threads, m)?)?;
    m.add_function(wrap_pyfunction!(init_thread_pool, m)?)?;

    register_errors(m.py(), m)?;

    Ok(())
}
