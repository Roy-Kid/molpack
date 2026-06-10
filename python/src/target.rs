//! Python wrapper for packing `Target`.
//!
//! [`PyTarget`] describes one type of molecule to pack: its template
//! geometry, topology, and the number of copies.
//!
//! The constructor accepts any frame-like object with an ``"atoms"`` block —
//! ``molrs.Frame``, ``molpy.Frame``, or a plain ``dict`` — and converts it to a
//! Rust [`molrs::Frame`] (see [`crate::frame_marshal`]). The full frame, with
//! topology, is handed to the core [`Target`], which owns the assembly. The
//! binding never re-derives the result frame; it only marshals data.

use crate::constraint::extract_restraint;
use crate::helpers::NpF;
use crate::types::{PyAngle, PyAxis, PyCenteringMode};
use molpack::F;
use molpack::target::Target;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

/// Build a [`Target`] from any frame-like Python object plus a copy count.
///
/// Shared by [`PyTarget::new`] and the script loader. The frame is converted to
/// a Rust [`molrs::Frame`] so the core retains its full topology.
pub(crate) fn target_from_frame(frame: &Bound<'_, PyAny>, count: usize) -> PyResult<Target> {
    let rust_frame = crate::frame_marshal::pyframe_to_rust(frame)?;
    let atoms = rust_frame
        .get("atoms")
        .ok_or_else(|| PyValueError::new_err(r#"frame must have an "atoms" block"#))?;
    if atoms.get("x").is_none() {
        return Err(PyValueError::new_err(
            r#"atoms block must have "x" / "y" / "z" columns"#,
        ));
    }
    Ok(Target::new(rust_frame, count))
}

#[pyclass(name = "Target", from_py_object)]
#[derive(Clone)]
pub struct PyTarget {
    pub(crate) inner: Target,
}

#[pymethods]
impl PyTarget {
    /// Create a packing target from a molecule frame.
    ///
    /// Parameters
    /// ----------
    /// frame : Frame-like
    ///     Any object with ``frame["atoms"]`` returning a block.
    ///     Supports ``molrs.Frame``, ``molpy.Frame``, and plain dicts.
    /// count : int
    ///     Number of copies to pack.
    ///
    /// Add a display name via :meth:`with_name`.
    #[new]
    #[pyo3(signature = (frame, count))]
    fn new(frame: &Bound<'_, PyAny>, count: usize) -> PyResult<Self> {
        Ok(PyTarget {
            inner: target_from_frame(frame, count)?,
        })
    }

    fn with_name(&self, name: &str) -> Self {
        PyTarget {
            inner: self.inner.clone().with_name(name),
        }
    }

    fn with_restraint(&self, restraint: &Bound<'_, pyo3::types::PyAny>) -> PyResult<Self> {
        let r = extract_restraint(restraint)?;
        Ok(PyTarget {
            inner: self.inner.clone().with_restraint(r),
        })
    }

    /// Attach a restraint to selected atoms of every copy.
    ///
    /// ``indices`` are **0-based** (Rust/Python native). Porting from a
    /// Packmol ``.inp`` file? Subtract 1 at the call site.
    fn with_atom_restraint(
        &self,
        indices: Vec<usize>,
        restraint: &Bound<'_, pyo3::types::PyAny>,
    ) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        let r = extract_restraint(restraint)?;
        Ok(PyTarget {
            inner: self.inner.clone().with_atom_restraint(&indices, r),
        })
    }

    fn with_perturb_budget(&self, budget: usize) -> Self {
        PyTarget {
            inner: self.inner.clone().with_perturb_budget(budget),
        }
    }

    fn with_centering(&self, mode: PyCenteringMode) -> Self {
        PyTarget {
            inner: self.inner.clone().with_centering(mode.into()),
        }
    }

    fn with_rotation_bound(&self, axis: PyAxis, center: PyAngle, half_width: PyAngle) -> Self {
        PyTarget {
            inner: self.inner.clone().with_rotation_bound(
                axis.into(),
                center.inner,
                half_width.inner,
            ),
        }
    }

    fn fixed_at(&self, position: [NpF; 3]) -> Self {
        PyTarget {
            inner: self.inner.clone().fixed_at(position),
        }
    }

    /// Apply an Euler orientation to a previously-fixed target.
    ///
    /// Takes a 3-tuple of :class:`Angle` — e.g.
    /// ``(Angle.from_degrees(45), Angle.ZERO, Angle.from_degrees(90))``.
    /// Must be called after :meth:`fixed_at`.
    fn with_orientation(&self, orientation: (PyAngle, PyAngle, PyAngle)) -> Self {
        PyTarget {
            inner: self.inner.clone().with_orientation([
                orientation.0.inner,
                orientation.1.inner,
                orientation.2.inner,
            ]),
        }
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }

    #[getter]
    fn natoms(&self) -> usize {
        self.inner.natoms()
    }

    #[getter]
    fn count(&self) -> usize {
        self.inner.count
    }

    #[getter]
    fn elements(&self) -> Vec<String> {
        self.inner.elements.clone()
    }

    #[getter]
    fn radii(&self) -> Vec<F> {
        self.inner.radii.clone()
    }

    #[getter]
    fn is_fixed(&self) -> bool {
        self.inner.fixed_at.is_some()
    }

    fn __repr__(&self) -> String {
        format!(
            "Target(natoms={}, count={}, name={:?})",
            self.inner.natoms(),
            self.inner.count,
            self.inner.name
        )
    }
}

fn validate_atom_indices(indices: &[usize], natoms: usize) -> PyResult<()> {
    for &index in indices {
        if index >= natoms {
            return Err(PyValueError::new_err(format!(
                "atom indices are 0-based and must be in 0..{natoms}, got {index}",
            )));
        }
    }
    Ok(())
}
