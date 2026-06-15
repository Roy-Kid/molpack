//! Python wrapper for packing `Target`.
//!
//! [`PyTarget`] describes one type of molecule to pack: its template
//! geometry, topology, and the number of copies.
//!
//! The constructor accepts a real molrs/molpy `Frame` (``molrs.Frame`` or
//! ``molpy.Frame``) carrying an ``"atoms"`` block. The frame crosses the
//! language boundary **zero-copy** through its stable-FFI capsule (see
//! [`crate::interop`]) — no dict marshalling, no consumer-side data type. The
//! full frame, with topology, is handed to the core [`Target`], which owns the
//! assembly.

use crate::constraint::{extract_collective_restraint, extract_restraint, try_atom_builtin};
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
    let rust_frame = crate::interop::owned_frame_from_py(frame)?;
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
    /// frame : molrs.Frame | molpy.Frame
    ///     A molrs/molpy frame with an ``"atoms"`` block (``x`` / ``y`` / ``z``
    ///     columns). Resolved zero-copy via its FFI capsule — a plain ``dict``
    ///     is no longer accepted; build a ``molrs.Frame`` first.
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

    /// Attach a restraint to this target — the single unified extension point.
    ///
    /// Accepts:
    ///
    /// * a built-in **geometric** restraint (:class:`InsideBoxRestraint`,
    ///   :class:`InsideSphereRestraint`, :class:`OutsideSphereRestraint`,
    ///   :class:`AbovePlaneRestraint`, :class:`BelowPlaneRestraint`) — its
    ///   ``f`` / ``fg`` see **one atom** at a time;
    /// * a built-in **distribution** restraint (:class:`GaussianPlane`,
    ///   :class:`GaussianPoint`, :class:`ExponentialPlane`,
    ///   :class:`ExponentialPoint`, :class:`TabulatedPlane`,
    ///   :class:`TabulatedPoint`);
    /// * any object with callable ``f`` / ``fg`` — the duck-typed extension
    ///   point.
    ///
    /// For the geometric built-ins ``f(x, scale, scale2)`` /
    /// ``fg(x, scale, scale2)`` see a single atom's ``(x, y, z)``. For the
    /// distribution and custom (duck-typed) restraints,
    /// ``f(coords, scale, scale2)`` / ``fg(coords, scale, scale2)`` see **every
    /// copy's** ``(x, y, z)`` (``coords`` is the full list) and ``fg`` returns
    /// ``(energy, [(gx, gy, gz), ...])`` — one gradient triple per copy.
    fn with_restraint(&self, restraint: &Bound<'_, pyo3::types::PyAny>) -> PyResult<Self> {
        if let Some(atom_r) = try_atom_builtin(restraint) {
            Ok(PyTarget {
                inner: self.inner.clone().with_restraint(atom_r),
            })
        } else {
            let group_r = extract_collective_restraint(restraint)?;
            Ok(PyTarget {
                inner: self.inner.clone().with_collective_restraint(group_r),
            })
        }
    }

    /// Attach an in-loop geometry relaxer (relaxation-assisted packing).
    ///
    /// Accepts either built-in relaxer:
    ///
    /// * :class:`TorsionMcRelaxer` — engine-free Monte-Carlo torsion sampling
    ///   (always available);
    /// * :class:`LBFGSRelaxer` — force-field L-BFGS minimization (`ff` feature).
    ///
    /// Requires ``count == 1`` — every copy shares the reference geometry the
    /// relaxer rewrites, so a relaxed target packs one molecule.
    fn with_relaxer(&self, relaxer: &Bound<'_, PyAny>) -> PyResult<Self> {
        if self.inner.count != 1 {
            return Err(PyValueError::new_err(format!(
                "with_relaxer requires count == 1 (all copies share the reference \
                 geometry the relaxer rewrites), got count = {}",
                self.inner.count
            )));
        }
        // Torsion-MC relaxer is core (no feature gate).
        if let Ok(tm) = relaxer.extract::<crate::relaxer::PyTorsionMcRelaxer>() {
            return Ok(PyTarget {
                inner: self.inner.clone().with_relaxer(tm.inner),
            });
        }
        #[cfg(feature = "ff")]
        if let Ok(lb) = relaxer.extract::<crate::relaxer::PyLBFGSRelaxer>() {
            return Ok(PyTarget {
                inner: self.inner.clone().with_relaxer(lb.inner),
            });
        }
        Err(PyValueError::new_err(
            "with_relaxer expects a TorsionMcRelaxer or LBFGSRelaxer",
        ))
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
