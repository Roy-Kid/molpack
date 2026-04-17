//! Shared helpers for the molpack PyO3 bindings.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::sync::Mutex;

/// Numpy float type matching molpack's `F = f64`.
pub type NpF = f64;

/// Convert a [`molpack::PackError`] to a Python `RuntimeError`.
pub fn pack_error_to_pyerr(e: molpack::PackError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
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
