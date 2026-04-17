//! PyO3 bindings for molpack.
//!
//! Python API:
//!
//! | Python class     | Rust wrapper         | Purpose                            |
//! |------------------|----------------------|------------------------------------|
//! | `Target`         | [`PyTarget`]         | Molecule specification for packing |
//! | `Molpack`        | [`PyPacker`]         | Molecular packer (Packmol port)    |
//! | `PackResult`     | [`PyPackResult`]     | Pack output (positions, elements)  |
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

mod helpers;

mod constraint;
use constraint::{PyAbovePlane, PyBelowPlane, PyInsideBox, PyInsideSphere, PyOutsideSphere};

mod handler;
use handler::PyStepInfo;

mod target;
use target::PyTarget;

mod packer;
use packer::{PyPackResult, PyPacker};

#[pymodule]
fn molpack(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyInsideBox>()?;
    m.add_class::<PyInsideSphere>()?;
    m.add_class::<PyOutsideSphere>()?;
    m.add_class::<PyAbovePlane>()?;
    m.add_class::<PyBelowPlane>()?;

    m.add_class::<PyTarget>()?;
    m.add_class::<PyPacker>()?;
    m.add_class::<PyPackResult>()?;
    m.add_class::<PyStepInfo>()?;

    Ok(())
}
