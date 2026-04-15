//! Shared helpers for the molpack PyO3 bindings.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

/// Numpy float type matching molpack's `F = f64`.
pub type NpF = f64;

/// Convert a [`molpack::PackError`] to a Python `RuntimeError`.
pub fn pack_error_to_pyerr(e: molpack::PackError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}
