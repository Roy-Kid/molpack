//! Typed bindings for the growth statistics inputs: `TorsionPrior` and
//! `AnglePrior`. The growth knobs themselves live on the `CbmcGrow` entry.

use molpack::grow::{AnglePrior, TorsionPrior};
use pyo3::prelude::*;

use molrs::op::F;

/// Geometric torsion prior — the mandatory statistics input of growth.
#[pyclass(name = "TorsionPrior", frozen, from_py_object)]
#[derive(Clone)]
pub struct PyTorsionPrior {
    pub(crate) inner: TorsionPrior,
}

#[pymethods]
impl PyTorsionPrior {
    /// Uniform on (−π, π]. Freely-rotating-chain statistics — negative
    /// control only; quantitatively wrong for melts.
    #[staticmethod]
    fn uniform() -> Self {
        Self {
            inner: TorsionPrior::Uniform,
        }
    }

    /// von-Mises-like spread of concentration `kappa` around the template's
    /// own torsion values.
    #[staticmethod]
    fn template(kappa: F) -> Self {
        Self {
            inner: TorsionPrior::Template { kappa },
        }
    }

    /// RIS-style discrete states as `(angle_rad, weight)` pairs.
    #[staticmethod]
    fn states(states: Vec<(F, F)>) -> Self {
        Self {
            inner: TorsionPrior::States(
                states.into_iter().map(|(a, w)| (a as F, w as F)).collect(),
            ),
        }
    }

    /// Three-state trans/gauche± prior calibrated from a target
    /// characteristic ratio (PEO: `three_state_from_c_inf(5.5, 1.9106)`).
    #[staticmethod]
    fn three_state_from_c_inf(c_inf: F, theta_rad: F) -> Self {
        Self {
            inner: TorsionPrior::three_state_from_c_inf(c_inf, theta_rad),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

/// Placement-angle prior. `template()` (the default) copies angles verbatim
/// — the all-atom behavior; `wlc*` is the CG persistence control.
#[pyclass(name = "AnglePrior", frozen, from_py_object)]
#[derive(Clone)]
pub struct PyAnglePrior {
    pub(crate) inner: AnglePrior,
}

#[pymethods]
impl PyAnglePrior {
    #[staticmethod]
    fn template() -> Self {
        Self {
            inner: AnglePrior::Template,
        }
    }

    /// Discrete worm-like chain with tilt `kappa`.
    #[staticmethod]
    fn wlc(kappa: F) -> Self {
        Self {
            inner: AnglePrior::Wlc { kappa },
        }
    }

    /// WLC tilt calibrated from a target characteristic ratio
    /// (Kremer–Grest melts: `wlc_from_c_inf(1.76)`).
    #[staticmethod]
    fn wlc_from_c_inf(c_inf: F) -> Self {
        Self {
            inner: AnglePrior::wlc_from_c_inf(c_inf),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
