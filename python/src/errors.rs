//! Error mapping for the molpack PyO3 bindings: the typed `PackError`
//! exception hierarchy, the Rust → Python error conversions, and the slot
//! that carries a Python exception raised inside a Rust-invoked callback.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::sync::Mutex;

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
    Ok(())
}

/// Convert a [`molpack::script::ScriptError`] to a Python exception.
///
/// Packing failures surfaced through the script loader fan out to the
/// same typed hierarchy as direct `pack()` calls; parse and I/O errors
/// become `ValueError` and `OSError` respectively.
pub fn script_error_to_pyerr(e: molpack::script::ScriptError) -> PyErr {
    use molpack::script::ScriptError;
    match e {
        ScriptError::Pack(p) => pack_error_to_pyerr(p),
        ScriptError::Io { .. } => pyo3::exceptions::PyOSError::new_err(e.to_string()),
        _ => pyo3::exceptions::PyValueError::new_err(e.to_string()),
    }
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
        // Triclinic / cell validation errors surface as ValueError until
        // dedicated Python exception types are added.
        molpack::PackError::InvalidCell { .. } => pyo3::exceptions::PyValueError::new_err(msg),
        molpack::PackError::RestraintAcrossPeriodicAxis { .. } => {
            pyo3::exceptions::PyValueError::new_err(msg)
        }
        // A short radius that is not shorter than the packing radius is a bad
        // input value, not a packing failure.
        molpack::PackError::ShortRadiusNotShorter { .. } => {
            pyo3::exceptions::PyValueError::new_err(msg)
        }
        // Growth and density declarations are input contracts: a target that
        // cannot be grown, a density fighting an explicit box, or a mass the
        // elements cannot resolve are all bad input values.
        // Templates whose shared columns disagree on dtype are bad input too.
        molpack::PackError::Grow { .. }
        | molpack::PackError::TemplateColumns { .. }
        | molpack::PackError::DensityConflictsWithBox
        | molpack::PackError::SeedMismatch { .. }
        | molpack::PackError::UnknownMass { .. } => pyo3::exceptions::PyValueError::new_err(msg),
        // A badly composed run — a stage chained where its precondition
        // cannot hold, a preset carrying a second set of shared settings, a
        // pipeline with no stages at all, or a guarded stage that left an
        // invariant broken. The wheel exposes no pipeline surface yet, so
        // these are unreachable from Python today; the arm exists so the
        // mapping stays total and the message is never lost.
        molpack::PackError::StageOrder { .. }
        | molpack::PackError::PresetSettingsInsidePipeline { .. }
        | molpack::PackError::InvariantViolated { .. }
        | molpack::PackError::NoStages => pyo3::exceptions::PyValueError::new_err(msg),
    }
}

/// Sink for Python exceptions raised inside Rust-invoked callbacks
/// (`PyCallableRestraint::fg`, `PyHandlerWrapper::on_*`). The Rust trait
/// signatures can't surface `PyErr` in-band, so callbacks stash the
/// first error here and set their stop-flag; the entry's `run()` drains
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
