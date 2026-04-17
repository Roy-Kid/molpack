//! Python wrappers for molecular packing restraints.
//!
//! Each built-in restraint (``InsideBox``, ``InsideSphere``, ``OutsideSphere``,
//! ``AbovePlane``, ``BelowPlane``) is a ``#[pyclass]`` with its own
//! parameters. Custom Python-defined restraints are supported via **duck
//! typing**: any object exposing callable ``f(x, scale, scale2)`` and
//! ``fg(x, scale, scale2)`` attributes may be passed to
//! ``Target.with_restraint``; see [`PyCallableRestraint`] for the contract.

use std::sync::Arc;

use crate::helpers::{NpF, stash_err};
use molpack::F;
use molpack::restraint::Restraint;
use molpack::restraint::{
    AbovePlaneRestraint, BelowPlaneRestraint, InsideBoxRestraint, InsideSphereRestraint,
    OutsideSphereRestraint,
};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

// Pass-through wrapper: an owned `Arc<dyn Restraint>` that itself
// implements `Restraint`, so it can be fed into
// `Target::with_restraint(impl Restraint)`. Adding a new restraint type
// only requires a new arm in `extract_restraint`.

#[derive(Clone)]
pub(crate) struct SharedRestraint(pub Arc<dyn Restraint>);

impl std::fmt::Debug for SharedRestraint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SharedRestraint({})", self.0.name())
    }
}

impl Restraint for SharedRestraint {
    #[inline]
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F {
        self.0.f(x, scale, scale2)
    }
    #[inline]
    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        self.0.fg(x, scale, scale2, g)
    }
    #[inline]
    fn is_parallel_safe(&self) -> bool {
        self.0.is_parallel_safe()
    }
    #[inline]
    fn name(&self) -> &'static str {
        self.0.name()
    }
}

// ============================================================================
// Extractor: try each built-in `#[pyclass]`, else duck-type on `f`/`fg`.
// ============================================================================

pub(crate) fn extract_restraint(
    obj: &Bound<'_, pyo3::types::PyAny>,
) -> PyResult<SharedRestraint> {
    if let Ok(c) = obj.extract::<PyInsideBox>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyInsideSphere>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyOutsideSphere>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyAbovePlane>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyBelowPlane>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }

    // Duck-typed Python restraint: object with callable `f` and `fg`
    // methods. Bound methods are resolved once here so the hot path
    // skips per-call attribute lookups.
    if let (Ok(f_method), Ok(fg_method)) = (obj.getattr("f"), obj.getattr("fg")) {
        return Ok(SharedRestraint(Arc::new(PyCallableRestraint {
            f_method: f_method.unbind(),
            fg_method: fg_method.unbind(),
        })));
    }

    Err(PyTypeError::new_err(
        "expected a restraint: one of InsideBox / InsideSphere / OutsideSphere / \
         AbovePlane / BelowPlane, or an object with callable `f(x, scale, scale2)` \
         and `fg(x, scale, scale2)` methods",
    ))
}

// PyCallableRestraint — bridge from the Rust `Restraint` trait to a
// Python object. Stores the **bound methods** directly (resolved at
// attach time) instead of the host object, because each restraint
// evaluation goes through `fg` inside the GENCAN inner loop — a per-atom
// string lookup there is measurable on larger systems.
//
// Python contract:
//   obj.f(x, scale, scale2)  -> float
//   obj.fg(x, scale, scale2) -> (float, (gx, gy, gz))
//
// `x` is passed as a 3-tuple; `fg`'s returned gradient is a flat tuple
// that Rust accumulates into `g` with `+=` on the caller's behalf.
//
// `is_parallel_safe() -> false` — the GIL serializes callbacks, and
// pretending otherwise would deadlock rayon reductions.

pub(crate) struct PyCallableRestraint {
    f_method: Py<pyo3::types::PyAny>,
    fg_method: Py<pyo3::types::PyAny>,
}

impl std::fmt::Debug for PyCallableRestraint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PyCallableRestraint")
    }
}

impl Restraint for PyCallableRestraint {
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F {
        Python::attach(|py| {
            let args = ((x[0], x[1], x[2]), scale, scale2);
            match self.f_method.bind(py).call1(args) {
                Ok(res) => match res.extract::<F>() {
                    Ok(v) => v,
                    Err(e) => {
                        stash_err(e);
                        0.0
                    }
                },
                Err(e) => {
                    stash_err(e);
                    0.0
                }
            }
        })
    }

    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        Python::attach(|py| {
            let args = ((x[0], x[1], x[2]), scale, scale2);
            match self.fg_method.bind(py).call1(args) {
                Ok(res) => match res.extract::<(F, (F, F, F))>() {
                    Ok((v, (gx, gy, gz))) => {
                        g[0] += gx;
                        g[1] += gy;
                        g[2] += gz;
                        v
                    }
                    Err(e) => {
                        stash_err(PyTypeError::new_err(format!(
                            "callable restraint `fg` must return (float, (gx, gy, gz)); {e}",
                        )));
                        0.0
                    }
                },
                Err(e) => {
                    stash_err(e);
                    0.0
                }
            }
        })
    }

    fn is_parallel_safe(&self) -> bool {
        false
    }

    fn name(&self) -> &'static str {
        "PyCallableRestraint"
    }
}

// ============================================================================
// Built-in restraint `#[pyclass]` wrappers (unchanged from prior version —
// they each hold a concrete `*Restraint` struct and expose a `__new__`).
// ============================================================================

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
