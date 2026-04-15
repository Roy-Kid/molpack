//! PyO3 bindings for molpack.
//!
//! Python API:
//!
//! | Python class         | Rust wrapper         | Purpose                             |
//! |----------------------|----------------------|-------------------------------------|
//! | `Target`             | [`PyTarget`]         | Molecule specification for packing  |
//! | `Packer`             | [`PyPacker`]         | Molecular packer (Packmol port)     |
//! | `PackResult`         | [`PyPackResult`]     | Pack output (positions, elements)   |
//! | `InsideBox`          | [`PyInsideBox`]      | Box restraint                       |
//! | `InsideSphere`       | [`PyInsideSphere`]   | Sphere restraint (inside)           |
//! | `OutsideSphere`      | [`PyOutsideSphere`]  | Sphere restraint (outside)          |
//! | `AbovePlane`         | [`PyAbovePlane`]     | Half-space restraint                |
//! | `BelowPlane`         | [`PyBelowPlane`]     | Half-space restraint                |
//! | `MoleculeConstraint` | [`PyMoleculeConstraint`] | Bundle of restraints            |

use pyo3::prelude::*;

mod helpers;

mod constraint;
use constraint::{
    PyAbovePlane, PyBelowPlane, PyInsideBox, PyInsideSphere, PyMoleculeConstraint, PyOutsideSphere,
};

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
    m.add_class::<PyMoleculeConstraint>()?;

    m.add_class::<PyTarget>()?;
    m.add_class::<PyPacker>()?;
    m.add_class::<PyPackResult>()?;

    Ok(())
}
