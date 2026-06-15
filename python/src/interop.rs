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

use molrs::Frame;
use molrs::spatial::region::simbox::SimBox;
use molrs_ffi::{FfiError, FrameRef};
use ndarray::Array1;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule, PyModule};

#[cfg(feature = "ff")]
use molrs_ffi::ForceFieldRef;

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
    let ptr = capsule.pointer_checked(Some(c"molrs.FrameRef"))?;
    let pp = ptr.as_ptr() as *const *const FrameRef;
    // SAFETY: a "molrs.FrameRef" capsule's void* is `*mut *mut FrameRef` (the
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

/// Resolve a `molrs.ForceField` / `molpy.ForceField` to a shared
/// [`ForceFieldRef`] (zero-copy). Consumed by `LBFGSRelaxer(ff)`.
#[cfg(feature = "ff")]
pub fn forcefield_from_py(obj: &Bound<'_, PyAny>) -> PyResult<ForceFieldRef> {
    let capsule = capsule_from(obj, "_ffi_forcefield_capsule")?;
    let ptr = capsule.pointer_checked(Some(c"molrs.ForceFieldRef"))?;
    let pp = ptr.as_ptr() as *const *const ForceFieldRef;
    // SAFETY: see `frame_from_py`; the payload is `*mut *mut ForceFieldRef`.
    let ffref = unsafe { (**pp).clone() };
    Ok(ffref)
}

/// Build a Python `molrs.Frame` from a Rust [`Frame`] — the **return path**.
///
/// Stamps an orthorhombic periodic box from `box_bounds` (`(min, max)` corners)
/// when present, wraps the frame in a fresh `FrameRef`, exports a
/// `"molrs.FrameRef"` capsule, and rebuilds it as a `molrs.Frame` through
/// `Frame._from_ffi_frameref_capsule`. No column marshalling.
pub fn frame_to_py<'py>(
    py: Python<'py>,
    frame: &Frame,
    box_bounds: Option<([F; 3], [F; 3])>,
) -> PyResult<Bound<'py, PyAny>> {
    let mut frame = frame.clone();
    if let Some((lo, hi)) = box_bounds {
        let lengths = Array1::from_vec(vec![hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]]);
        let origin = Array1::from_vec(lo.to_vec());
        let simbox = SimBox::ortho(lengths, origin, [true, true, true])
            .map_err(|e| PyValueError::new_err(format!("building periodic box: {e:?}")))?;
        frame.simbox = Some(simbox);
    }
    // Move the result frame into a fresh standalone store and export its capsule.
    let fref = FrameRef::new_standalone();
    fref.with_mut(|slot| *slot = frame).map_err(ffi_err)?;
    let capsule = export_frame_capsule(py, fref)?;
    let molrs = PyModule::import(py, "molrs")?;
    molrs
        .getattr("Frame")?
        .call_method1("_from_ffi_frameref_capsule", (capsule,))
}

/// Box a `FrameRef` into a `"molrs.FrameRef"` PyCapsule — the exporter side of
/// the return path, mirroring molrs-python's `Frame._ffi_frameref_capsule`.
fn export_frame_capsule<'py>(py: Python<'py>, fref: FrameRef) -> PyResult<Bound<'py, PyCapsule>> {
    let raw = FrameRefPtr(Box::into_raw(Box::new(fref)));
    let name = std::ffi::CString::new("molrs.FrameRef").expect("static capsule name");
    PyCapsule::new_with_destructor(py, raw, Some(name), |ptr: FrameRefPtr, _ctx| {
        // SAFETY: `ptr.0` came from `Box::into_raw` above and is reclaimed
        // exactly once when the capsule dies.
        drop(unsafe { Box::from_raw(ptr.0) });
    })
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
