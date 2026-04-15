//! Python wrappers for molecular packing restraints.
//!
//! Each concrete restraint (``InsideBox``, ``InsideSphere``, ``OutsideSphere``,
//! ``AbovePlane``, ``BelowPlane``) can be attached to a target via
//! ``Target.with_restraint`` or ``Target.with_restraint_for_atoms``.
//! Multiple restraints may be stacked by calling these methods repeatedly.

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

pub(crate) fn extract_restraint(
    obj: &Bound<'_, pyo3::types::PyAny>,
) -> PyResult<AnyRestraint> {
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
        "expected a restraint (InsideBox, InsideSphere, OutsideSphere, AbovePlane, BelowPlane)",
    ))
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

    fn __repr__(&self) -> String {
        "InsideBox(...)".to_string()
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

    fn __repr__(&self) -> String {
        "InsideSphere(...)".to_string()
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

    fn __repr__(&self) -> String {
        "OutsideSphere(...)".to_string()
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

    fn __repr__(&self) -> String {
        "AbovePlane(...)".to_string()
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

    fn __repr__(&self) -> String {
        "BelowPlane(...)".to_string()
    }
}
