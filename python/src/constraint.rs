//! Python wrappers for molecular packing restraints.
//!
//! Each concrete restraint (``InsideBox``, ``InsideSphere``, ``OutsideSphere``,
//! ``AbovePlane``, ``BelowPlane``) can be composed via ``.and_()`` to build a
//! ``MoleculeConstraint`` — a bundle of restraints applied together.
//!
//! When a ``MoleculeConstraint`` (or a single restraint) is passed to
//! ``Target.with_constraint`` / ``with_constraint_for_atoms``, every restraint
//! in the bundle is attached to the target independently.

use crate::helpers::NpF;
use molpack::restraint::Restraint;
use molpack::restraint::{
    AbovePlaneRestraint, BelowPlaneRestraint, InsideBoxRestraint, InsideSphereRestraint,
    OutsideSphereRestraint,
};
use pyo3::prelude::*;

#[derive(Clone, Debug)]
pub(crate) enum AnyRestraint {
    InsideBox(InsideBoxRestraint),
    InsideSphere(InsideSphereRestraint),
    OutsideSphere(OutsideSphereRestraint),
    AbovePlane(AbovePlaneRestraint),
    BelowPlane(BelowPlaneRestraint),
}

impl Restraint for AnyRestraint {
    fn f(&self, x: &[molpack::F; 3], scale: molpack::F, scale2: molpack::F) -> molpack::F {
        match self {
            Self::InsideBox(r) => r.f(x, scale, scale2),
            Self::InsideSphere(r) => r.f(x, scale, scale2),
            Self::OutsideSphere(r) => r.f(x, scale, scale2),
            Self::AbovePlane(r) => r.f(x, scale, scale2),
            Self::BelowPlane(r) => r.f(x, scale, scale2),
        }
    }
    fn fg(
        &self,
        x: &[molpack::F; 3],
        scale: molpack::F,
        scale2: molpack::F,
        g: &mut [molpack::F; 3],
    ) -> molpack::F {
        match self {
            Self::InsideBox(r) => r.fg(x, scale, scale2, g),
            Self::InsideSphere(r) => r.fg(x, scale, scale2, g),
            Self::OutsideSphere(r) => r.fg(x, scale, scale2, g),
            Self::AbovePlane(r) => r.fg(x, scale, scale2, g),
            Self::BelowPlane(r) => r.fg(x, scale, scale2, g),
        }
    }
}

fn extract_single(obj: &Bound<'_, pyo3::types::PyAny>) -> PyResult<AnyRestraint> {
    if let Ok(c) = obj.extract::<PyInsideBox>() {
        return Ok(AnyRestraint::InsideBox(c.inner));
    }
    if let Ok(c) = obj.extract::<PyInsideSphere>() {
        return Ok(AnyRestraint::InsideSphere(c.inner));
    }
    if let Ok(c) = obj.extract::<PyOutsideSphere>() {
        return Ok(AnyRestraint::OutsideSphere(c.inner));
    }
    if let Ok(c) = obj.extract::<PyAbovePlane>() {
        return Ok(AnyRestraint::AbovePlane(c.inner));
    }
    if let Ok(c) = obj.extract::<PyBelowPlane>() {
        return Ok(AnyRestraint::BelowPlane(c.inner));
    }
    Err(pyo3::exceptions::PyTypeError::new_err(
        "expected a restraint (InsideBox, InsideSphere, OutsideSphere, AbovePlane, BelowPlane, MoleculeConstraint)",
    ))
}

pub(crate) fn extract_restraints(
    obj: &Bound<'_, pyo3::types::PyAny>,
) -> PyResult<Vec<AnyRestraint>> {
    if let Ok(mc) = obj.extract::<PyMoleculeConstraint>() {
        return Ok(mc.restraints);
    }
    Ok(vec![extract_single(obj)?])
}

macro_rules! restraint_and {
    ($self_:expr, $other:expr) => {{
        let mut rs = vec![AnyRestraint::from($self_.inner.clone())];
        rs.extend(extract_restraints($other)?);
        Ok(PyMoleculeConstraint { restraints: rs })
    }};
}

#[pyclass(name = "InsideBox", from_py_object)]
#[derive(Clone)]
pub struct PyInsideBox {
    pub(crate) inner: InsideBoxRestraint,
}

#[pymethods]
impl PyInsideBox {
    #[new]
    fn new(min: [NpF; 3], max: [NpF; 3]) -> Self {
        Self {
            inner: InsideBoxRestraint::new(min, max),
        }
    }

    #[pyo3(name = "and_")]
    fn and_(&self, other: &Bound<'_, pyo3::types::PyAny>) -> PyResult<PyMoleculeConstraint> {
        restraint_and!(self, other)
    }

    fn __repr__(&self) -> String {
        "InsideBox(...)".to_string()
    }
}

impl From<InsideBoxRestraint> for AnyRestraint {
    fn from(r: InsideBoxRestraint) -> Self {
        AnyRestraint::InsideBox(r)
    }
}

#[pyclass(name = "InsideSphere", from_py_object)]
#[derive(Clone)]
pub struct PyInsideSphere {
    pub(crate) inner: InsideSphereRestraint,
}

#[pymethods]
impl PyInsideSphere {
    #[new]
    fn new(radius: NpF, center: [NpF; 3]) -> Self {
        Self {
            inner: InsideSphereRestraint::new(center, radius),
        }
    }

    #[pyo3(name = "and_")]
    fn and_(&self, other: &Bound<'_, pyo3::types::PyAny>) -> PyResult<PyMoleculeConstraint> {
        restraint_and!(self, other)
    }

    fn __repr__(&self) -> String {
        "InsideSphere(...)".to_string()
    }
}

impl From<InsideSphereRestraint> for AnyRestraint {
    fn from(r: InsideSphereRestraint) -> Self {
        AnyRestraint::InsideSphere(r)
    }
}

#[pyclass(name = "OutsideSphere", from_py_object)]
#[derive(Clone)]
pub struct PyOutsideSphere {
    pub(crate) inner: OutsideSphereRestraint,
}

#[pymethods]
impl PyOutsideSphere {
    #[new]
    fn new(radius: NpF, center: [NpF; 3]) -> Self {
        Self {
            inner: OutsideSphereRestraint::new(center, radius),
        }
    }

    #[pyo3(name = "and_")]
    fn and_(&self, other: &Bound<'_, pyo3::types::PyAny>) -> PyResult<PyMoleculeConstraint> {
        restraint_and!(self, other)
    }

    fn __repr__(&self) -> String {
        "OutsideSphere(...)".to_string()
    }
}

impl From<OutsideSphereRestraint> for AnyRestraint {
    fn from(r: OutsideSphereRestraint) -> Self {
        AnyRestraint::OutsideSphere(r)
    }
}

#[pyclass(name = "AbovePlane", from_py_object)]
#[derive(Clone)]
pub struct PyAbovePlane {
    pub(crate) inner: AbovePlaneRestraint,
}

#[pymethods]
impl PyAbovePlane {
    #[new]
    fn new(normal: [NpF; 3], distance: NpF) -> Self {
        Self {
            inner: AbovePlaneRestraint::new(normal, distance),
        }
    }

    #[pyo3(name = "and_")]
    fn and_(&self, other: &Bound<'_, pyo3::types::PyAny>) -> PyResult<PyMoleculeConstraint> {
        restraint_and!(self, other)
    }

    fn __repr__(&self) -> String {
        "AbovePlane(...)".to_string()
    }
}

impl From<AbovePlaneRestraint> for AnyRestraint {
    fn from(r: AbovePlaneRestraint) -> Self {
        AnyRestraint::AbovePlane(r)
    }
}

#[pyclass(name = "BelowPlane", from_py_object)]
#[derive(Clone)]
pub struct PyBelowPlane {
    pub(crate) inner: BelowPlaneRestraint,
}

#[pymethods]
impl PyBelowPlane {
    #[new]
    fn new(normal: [NpF; 3], distance: NpF) -> Self {
        Self {
            inner: BelowPlaneRestraint::new(normal, distance),
        }
    }

    #[pyo3(name = "and_")]
    fn and_(&self, other: &Bound<'_, pyo3::types::PyAny>) -> PyResult<PyMoleculeConstraint> {
        restraint_and!(self, other)
    }

    fn __repr__(&self) -> String {
        "BelowPlane(...)".to_string()
    }
}

impl From<BelowPlaneRestraint> for AnyRestraint {
    fn from(r: BelowPlaneRestraint) -> Self {
        AnyRestraint::BelowPlane(r)
    }
}

#[pyclass(name = "MoleculeConstraint", from_py_object)]
#[derive(Clone)]
pub struct PyMoleculeConstraint {
    pub(crate) restraints: Vec<AnyRestraint>,
}

#[pymethods]
impl PyMoleculeConstraint {
    #[pyo3(name = "and_")]
    fn and_(&self, other: &Bound<'_, pyo3::types::PyAny>) -> PyResult<PyMoleculeConstraint> {
        let mut rs = self.restraints.clone();
        rs.extend(extract_restraints(other)?);
        Ok(PyMoleculeConstraint { restraints: rs })
    }

    fn __repr__(&self) -> String {
        format!("MoleculeConstraint(restraints={})", self.restraints.len())
    }
}
