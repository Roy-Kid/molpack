//! Zero-copy interop between molrs / molpy Python objects and molrs-ffi handles.
//!
//! molrs and molpack are **separate** PyO3 extensions, so a `molrs.Frame`
//! pyclass cannot be `.extract()`d into a Rust value here. Instead molrs-python
//! exposes a stable-FFI capsule — `Frame._ffi_frameref_capsule()` and
//! `ForceField._ffi_forcefield_capsule()` — and molpack resolves it to the
//! shared `molrs_ffi::{FrameRef, ForceFieldRef}` handle. This is the exact
//! pattern molrs-cxxapi uses (`frame_clone_from_addr`): **no dict marshalling,
//! no consumer-side data type** (there is no `mpk.Frame`).
//!
//! Soundness: both wheels link the same `molcrafts-molrs-ffi` and the same
//! always-on `molcrafts-molrs` core, whose `Frame` / `Block` / `SimBox` layout
//! is feature-independent. So a handle minted by molrs-python (built with the
//! `full` feature set) and the `molrs::Frame` it lends are layout-identical to
//! what molpack (built `ff`-only) sees across the extension boundary.
//!
//! The version contract is **minor-line = ABI version** (`molrs_ffi::abi`):
//! layout is frozen within a molrs minor line (enforced by molrs-ffi's layout
//! snapshot gate), so both wheels must embed the same `major.minor` — patch
//! may drift. Two gates enforce it here: [`check_abi`] compares
//! `molrs._ffi_abi_token()` against the embedded line at import (clear
//! `ImportError`), and the capsule names carry the line
//! (`molrs.FrameRef/<major.minor>`), so even a stale consumer fails the name
//! check instead of dereferencing a drifted layout.

use molrs::Frame;
use molrs::spatial::simbox::SimBox;
use molrs_ffi::{FfiError, FrameRef};
use ndarray::Array1;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule, PyModule};

use molpack::F;

/// Map a molrs-ffi handle error into a Python exception.
fn ffi_err(e: FfiError) -> PyErr {
    PyTypeError::new_err(format!("molrs FFI error: {e}"))
}

/// Resolve a `molrs.Frame` / `molpy.Frame` to a shared [`FrameRef`] (zero-copy).
///
/// Clones the handle the capsule carries (two `Rc` bumps) onto the same store,
/// so reads/writes through the returned handle are visible in the originating
/// Python frame. The object must expose `_ffi_frameref_capsule()` — i.e. be a
/// real molrs/molpy `Frame` (a plain `dict` is no longer accepted).
pub fn frame_from_py(obj: &Bound<'_, PyAny>) -> PyResult<FrameRef> {
    let capsule = capsule_from(obj, "_ffi_frameref_capsule")?;
    // The expected name carries the ABI line of the molrs this wheel embeds
    // (`molrs.FrameRef/<major.minor>`), so a producer on another minor line
    // fails here cleanly instead of being dereferenced.
    let expected = molrs_ffi::abi::frameref_capsule_name();
    let ptr = capsule.pointer_checked(Some(expected)).map_err(|err| {
        PyValueError::new_err(format!(
            "{err} — molpack embeds molrs ABI line {line} (capsule name \
             {expected:?}); the producing molrs/molpy wheel is on a different \
             minor line. Align molcrafts-molrs and molcrafts-molpack on one \
             minor line. / molpack 与 molrs 的 minor 版本线不一致，请对齐后重装。",
            line = molrs_ffi::abi::abi_line(),
        ))
    })?;
    let pp = ptr.as_ptr() as *const *const FrameRef;
    // SAFETY: the versioned capsule's void* is `*mut *mut FrameRef` (the
    // exporter boxes a `*mut FrameRef`); deref twice to reach the cloned handle
    // and `.clone()` it (Rc bumps). The capsule is only touched under the GIL.
    let fref = unsafe { (**pp).clone() };
    Ok(fref)
}

/// Resolve a frame-like Python object to an **owned** core [`Frame`], deep-copied
/// out of the shared store.
///
/// The common case for consumers that keep the frame (target assembly, potential
/// compilation) — equivalent to `frame_from_py(obj)?.clone_frame()`.
pub fn owned_frame_from_py(obj: &Bound<'_, PyAny>) -> PyResult<Frame> {
    frame_from_py(obj)?.clone_frame().map_err(ffi_err)
}

/// Stamp an orthorhombic periodic box from `(min, max)` corners onto *frame*.
///
/// No-op when *box_bounds* is `None`. Call once before [`frame_to_py`].
pub fn stamp_box_bounds(frame: &mut Frame, box_bounds: Option<([F; 3], [F; 3])>) -> PyResult<()> {
    let Some((lo, hi)) = box_bounds else {
        return Ok(());
    };
    let lengths = Array1::from_vec(vec![hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]]);
    let origin = Array1::from_vec(lo.to_vec());
    let simbox = SimBox::ortho(lengths, origin, [true, true, true])
        .map_err(|e| PyValueError::new_err(format!("building periodic box: {e:?}")))?;
    frame.simbox = Some(simbox);
    Ok(())
}

/// One-shot Rust ``Frame`` → Python ``molrs.Frame``.
///
/// Two cdylibs cannot share a pyclass, so this clones into a standalone
/// ``FrameRef`` capsule. Do it once and keep the Python object.
pub fn frame_to_py<'py>(py: Python<'py>, frame: &Frame) -> PyResult<Bound<'py, PyAny>> {
    let fref = FrameRef::new_standalone();
    fref.with_mut(|slot| *slot = frame.clone())
        .map_err(ffi_err)?;
    let capsule = export_frame_capsule(py, fref)?;
    let molrs = PyModule::import(py, "molrs")?;
    molrs
        .getattr("Frame")?
        .call_method1("_from_ffi_frameref_capsule", (capsule,))
}

/// Box a `FrameRef` into a versioned `molrs.FrameRef/<major.minor>` PyCapsule
/// — the exporter side of the return path, mirroring molrs-python's
/// `Frame._ffi_frameref_capsule`.
fn export_frame_capsule<'py>(py: Python<'py>, fref: FrameRef) -> PyResult<Bound<'py, PyCapsule>> {
    // Same constructor and pointer shape as molrs-python's
    // ``Frame._ffi_frameref_capsule``: PyO3 boxes the ``FrameRefPtr`` payload,
    // so the capsule void* is ``*mut *mut FrameRef``. The shared abi module is
    // the single source of the (``&'static CStr``) name — never hard-code it.
    let raw = FrameRefPtr(Box::into_raw(Box::new(fref)));
    let name = molrs_ffi::abi::frameref_capsule_name();
    PyCapsule::new_with_value_and_destructor(py, raw, name, |ptr: FrameRefPtr, _ctx| {
        // SAFETY: `ptr.0` came from `Box::into_raw` above and is reclaimed
        // exactly once when the capsule dies.
        drop(unsafe { Box::from_raw(ptr.0) });
    })
}

/// Import-time ABI handshake against the installed `molcrafts-molrs` wheel.
///
/// Calls `molrs._ffi_abi_token()` and compares its ABI line against the line
/// molpack embeds. Runs once from the `#[pymodule]` init so a minor-line
/// mismatch is a clear `ImportError` naming both versions, not a later
/// capsule-name `ValueError` deep inside a pack run.
pub fn check_abi(py: Python<'_>) -> PyResult<()> {
    let embedded = molrs_ffi::abi::abi_line();
    let molrs = PyModule::import(py, "molrs")?;
    let token = match molrs.getattr("_ffi_abi_token") {
        Ok(f) => f.call0()?,
        Err(_) => {
            // Pre-0.14 wheels have no handshake — they are on an older line
            // by definition (the token and the versioned capsule names were
            // introduced together).
            return Err(pyo3::exceptions::PyImportError::new_err(format!(
                "molpack embeds molrs ABI line {embedded}, but the installed \
                 molcrafts-molrs predates the ABI handshake (≤0.13). Install \
                 a matching wheel: pip install 'molcrafts-molrs>={embedded}.0,\
                 <{next}' / 已安装的 molcrafts-molrs 过旧，请安装 {embedded}.* 版本。",
                next = next_minor(embedded),
            )));
        }
    };
    let (line, version): (String, String) = token
        .cast::<pyo3::types::PyTuple>()
        .map_err(|_| {
            pyo3::exceptions::PyImportError::new_err("molrs._ffi_abi_token() returned a non-tuple")
        })
        .and_then(|t| Ok((t.get_item(0)?.extract()?, t.get_item(1)?.extract()?)))?;
    if line != embedded {
        return Err(pyo3::exceptions::PyImportError::new_err(format!(
            "Minor-line mismatch: molpack embeds molrs ABI line {embedded}, \
             but the installed molcrafts-molrs is {version} (line {line}). \
             Handles cannot cross minor lines — install a matching wheel: \
             pip install 'molcrafts-molrs>={embedded}.0,<{next}' / molpack \
             与已安装的 molcrafts-molrs({version})minor 版本线不一致，请对齐。",
            next = next_minor(embedded),
        )));
    }
    Ok(())
}

/// `"0.14"` → `"0.15"` — the exclusive upper bound of a minor line, for pip
/// range hints in handshake errors.
fn next_minor(line: &str) -> String {
    let (major, minor) = line.split_once('.').unwrap_or((line, "0"));
    let bumped = minor.parse::<u64>().map(|m| m + 1).unwrap_or(0);
    format!("{major}.{bumped}")
}

/// `Send` wrapper around a `*mut FrameRef` for the capsule payload (mirrors
/// molrs-python's `FrameRefPtr`).
///
/// `FrameRef` is `!Send` (holds an `Rc`); the capsule is only ever created, read,
/// and destroyed under the Python GIL, so the single-threaded discipline holds.
/// `#[repr(transparent)]` makes the capsule's `void*` a `*mut *mut FrameRef`,
/// the shape `Frame._from_ffi_frameref_capsule` resolves.
#[repr(transparent)]
struct FrameRefPtr(*mut FrameRef);

// SAFETY: GIL-guarded, single-threaded use only — see the type-level doc.
unsafe impl Send for FrameRefPtr {}

/// Call `obj.<method>()` and downcast the result to a `PyCapsule`, with a clear
/// error when the object is not a molrs/molpy `Frame` / `ForceField`.
fn capsule_from<'py>(obj: &Bound<'py, PyAny>, method: &str) -> PyResult<Bound<'py, PyCapsule>> {
    let cap = obj.call_method0(method).map_err(|e| {
        PyTypeError::new_err(format!(
            "expected a molrs/molpy object exposing {method}(): {e}"
        ))
    })?;
    cap.cast_into::<PyCapsule>()
        .map_err(|_| PyTypeError::new_err(format!("{method}() did not return a PyCapsule")))
}
