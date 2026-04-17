//! Python wrappers for cross-cutting typed values: `Angle`, `Axis`,
//! `CenteringMode`. All three mirror the Rust types defined in
//! `molpack::target` 1:1.

use molpack::target::{Angle, Axis, CenteringMode};
use pyo3::prelude::*;

// ── Angle ──────────────────────────────────────────────────────────────────

/// Angular quantity with explicit units at the call site.
///
/// ```python
/// from molpack import Angle
/// Angle.from_degrees(30.0)
/// Angle.from_radians(3.14159 / 6)
/// ```
#[pyclass(name = "Angle", frozen, eq, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub struct PyAngle {
    pub(crate) inner: Angle,
}

#[pymethods]
impl PyAngle {
    /// Construct an angle from degrees.
    #[classmethod]
    #[pyo3(signature = (deg))]
    fn from_degrees(_cls: &Bound<'_, pyo3::types::PyType>, deg: f64) -> Self {
        Self {
            inner: Angle::from_degrees(deg),
        }
    }

    /// Construct an angle from radians.
    #[classmethod]
    #[pyo3(signature = (rad))]
    fn from_radians(_cls: &Bound<'_, pyo3::types::PyType>, rad: f64) -> Self {
        Self {
            inner: Angle::from_radians(rad),
        }
    }

    /// Zero rotation. Exposed as `Angle.ZERO`.
    #[classattr]
    #[allow(non_snake_case)]
    fn ZERO() -> Self {
        Self { inner: Angle::ZERO }
    }

    #[getter]
    fn degrees(&self) -> f64 {
        self.inner.degrees()
    }

    #[getter]
    fn radians(&self) -> f64 {
        self.inner.radians()
    }

    fn __repr__(&self) -> String {
        format!("Angle.from_degrees({})", self.inner.degrees())
    }
}

// ── Axis ───────────────────────────────────────────────────────────────────

/// Cartesian axis selector.
///
/// ```python
/// from molpack import Axis
/// Axis.X / Axis.Y / Axis.Z
/// ```
#[pyclass(name = "Axis", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum PyAxis {
    X,
    Y,
    Z,
}

impl From<PyAxis> for Axis {
    fn from(v: PyAxis) -> Self {
        match v {
            PyAxis::X => Axis::X,
            PyAxis::Y => Axis::Y,
            PyAxis::Z => Axis::Z,
        }
    }
}

// ── CenteringMode ──────────────────────────────────────────────────────────

/// Centering behavior for a target's reference coordinates.
///
/// - ``AUTO``  : free targets centered, fixed targets kept in place (default).
/// - ``CENTER``: always center.
/// - ``OFF``   : keep input coordinates unchanged.
#[pyclass(name = "CenteringMode", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
#[allow(clippy::upper_case_acronyms)]
pub enum PyCenteringMode {
    AUTO,
    CENTER,
    OFF,
}

impl From<PyCenteringMode> for CenteringMode {
    fn from(v: PyCenteringMode) -> Self {
        match v {
            PyCenteringMode::AUTO => CenteringMode::Auto,
            PyCenteringMode::CENTER => CenteringMode::Center,
            PyCenteringMode::OFF => CenteringMode::Off,
        }
    }
}
