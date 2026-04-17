//! Shared helpers for the molpack PyO3 bindings.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::sync::Mutex;

/// Numpy float type matching molpack's `F = f64`.
pub type NpF = f64;

// ── Typed exception hierarchy ──────────────────────────────────────────────
//
// Rooted at `PackError` so callers can catch any packing failure with a
// single ``except molpack.PackError`` clause. Leaf types mirror the Rust
// `PackError` variants and are the canonical exception classes users should
// match against.

pyo3::create_exception!(molpack, PackError, PyRuntimeError);
pyo3::create_exception!(molpack, ConstraintsFailedError, PackError);
pyo3::create_exception!(molpack, MaxIterationsError, PackError);
pyo3::create_exception!(molpack, NoTargetsError, PackError);
pyo3::create_exception!(molpack, EmptyMoleculeError, PackError);
pyo3::create_exception!(molpack, InvalidPBCBoxError, PackError);
pyo3::create_exception!(molpack, ConflictingPeriodicBoxesError, PackError);

/// Register all `PackError` subclasses on a module.
pub fn register_errors(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("PackError", py.get_type::<PackError>())?;
    m.add(
        "ConstraintsFailedError",
        py.get_type::<ConstraintsFailedError>(),
    )?;
    m.add("MaxIterationsError", py.get_type::<MaxIterationsError>())?;
    m.add("NoTargetsError", py.get_type::<NoTargetsError>())?;
    m.add("EmptyMoleculeError", py.get_type::<EmptyMoleculeError>())?;
    m.add("InvalidPBCBoxError", py.get_type::<InvalidPBCBoxError>())?;
    m.add(
        "ConflictingPeriodicBoxesError",
        py.get_type::<ConflictingPeriodicBoxesError>(),
    )?;
    Ok(())
}

/// Convert a [`molpack::PackError`] to the matching typed Python exception.
pub fn pack_error_to_pyerr(e: molpack::PackError) -> PyErr {
    let msg = e.to_string();
    match e {
        molpack::PackError::ConstraintsFailed(_) => ConstraintsFailedError::new_err(msg),
        molpack::PackError::MaxIterations => MaxIterationsError::new_err(msg),
        molpack::PackError::NoTargets => NoTargetsError::new_err(msg),
        molpack::PackError::EmptyMolecule(_) => EmptyMoleculeError::new_err(msg),
        molpack::PackError::InvalidPBCBox { .. } => InvalidPBCBoxError::new_err(msg),
        molpack::PackError::ConflictingPeriodicBoxes { .. } => {
            ConflictingPeriodicBoxesError::new_err(msg)
        }
    }
}

/// Sink for Python exceptions raised inside Rust-invoked callbacks
/// (`PyCallableRestraint::fg`, `PyHandlerWrapper::on_*`). The Rust trait
/// signatures can't surface `PyErr` in-band, so callbacks stash the
/// first error here and set their stop-flag; `Molpack.pack()` drains
/// the slot at return time and re-raises.
///
/// `PyErr` is `!Send` alone (needs GIL to drop), but `Mutex<Option<PyErr>>`
/// is `Send + Sync`, so this works as a global even though a plain
/// `Cell` would not. Concurrent `pack()` calls from different threads
/// are not supported.
static PACK_ERR: Mutex<Option<PyErr>> = Mutex::new(None);

/// Record a Python error raised inside a Rust-invoked callback. Only the
/// first error per `pack()` invocation is kept so the user sees the root
/// cause, not a cascade.
pub fn stash_err(e: PyErr) {
    let mut slot = PACK_ERR.lock().expect("PACK_ERR poisoned");
    if slot.is_none() {
        *slot = Some(e);
    }
}

/// Take the stashed error, leaving the slot empty.
pub fn take_err() -> Option<PyErr> {
    PACK_ERR.lock().expect("PACK_ERR poisoned").take()
}
