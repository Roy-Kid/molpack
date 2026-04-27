//! Python wrapper for packing `Target`.
//!
//! [`PyTarget`] describes one type of molecule to pack: its template
//! geometry, element symbols, and the number of copies.
//!
//! The constructor accepts any ``frame``-like object with a ``"atoms"`` block:
//!
//! - ``molrs.Frame`` (from ``molrs.read_pdb`` / ``read_xyz``) — block columns
//!   are accessed via ``.view(name)`` since ``molrs.Block`` has no ``__getitem__``.
//! - ``molpy.Frame`` / plain dicts — columns accessed via ``block[name]``.
//!
//! Rust calls the appropriate Python method by probing ``.view()`` first, then
//! falling back to ``[]``.  The element column is tried as ``"element"`` first
//! (molpy, XYZ), then ``"symbol"`` (molrs PDB).

use crate::constraint::extract_restraint;
use crate::helpers::NpF;
use crate::types::{PyAngle, PyAxis, PyCenteringMode};
use molpack::F;
use molpack::target::Target;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

/// Bondi (1964) van der Waals radius for a given element symbol (Å).
/// Falls back to 1.50 Å for unknowns.
pub(crate) fn vdw_radius_for(element: &str) -> NpF {
    match element.to_uppercase().as_str() {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "P" => 1.80,
        "S" => 1.80,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        "NA" => 2.27,
        "K" => 2.75,
        "CA" => 1.73,
        "MG" => 1.73,
        "ZN" => 1.39,
        "FE" => 1.94,
        _ => 1.50,
    }
}

/// Build a [`Target`] from any frame-like Python object plus a copy count.
///
/// Used by both [`PyTarget::new`] and the script loader so they share one
/// duck-typed extraction path. The frame may be a [`molrs.Frame`] (block
/// columns via `.view(name)`), a `molpy.Frame`, or a plain dict.
pub(crate) fn target_from_frame(frame: &Bound<'_, PyAny>, count: usize) -> PyResult<Target> {
    let atoms = frame
        .get_item("atoms")
        .map_err(|_| PyValueError::new_err(r#"frame must have an "atoms" block"#))?;

    let x: Vec<NpF> = read_column(&atoms, "x")?;
    let y: Vec<NpF> = read_column(&atoms, "y")?;
    let z: Vec<NpF> = read_column(&atoms, "z")?;
    let elements: Vec<String> = read_element_column(&atoms)?;

    let n = x.len();
    if y.len() != n || z.len() != n || elements.len() != n {
        return Err(PyValueError::new_err(
            "x, y, z, and element arrays must all have the same length",
        ));
    }

    let coords: Vec<[F; 3]> = x
        .iter()
        .zip(y.iter())
        .zip(z.iter())
        .map(|((xi, yi), zi)| [*xi, *yi, *zi])
        .collect();
    let radii: Vec<F> = elements.iter().map(|e| vdw_radius_for(e)).collect();

    let mut target = Target::from_coords(&coords, &radii, count);
    target.elements = elements;
    Ok(target)
}

/// Read a numeric column from an atoms block.
///
/// Tries `block.view(name)` first (``molrs.Block``), then ``block[name]``
/// (``molpy.Block`` / plain dict).
fn read_column(block: &Bound<'_, PyAny>, name: &str) -> PyResult<Vec<NpF>> {
    if let Ok(col) = block.call_method1("view", (name,)) {
        return col.extract();
    }
    block
        .get_item(name)
        .map_err(|_| PyValueError::new_err(format!(r#"atoms block has no "{name}" column"#)))?
        .extract()
}

/// Read the element-symbol column from an atoms block.
///
/// Tries the columns ``"element"`` and ``"symbol"`` in order, using both
/// ``.view()`` (molrs) and ``[]`` (molpy / dict) access per column.
///
/// | Source          | Column    | Access  |
/// |-----------------|-----------|---------|
/// | molrs PDB       | ``symbol``  | `.view()` |
/// | molrs XYZ       | ``element`` | `.view()` |
/// | molpy / dict    | ``element`` | `[]`    |
fn read_element_column(block: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
    for col in ["element", "symbol"] {
        if let Ok(arr) = block.call_method1("view", (col,)) {
            return arr.extract();
        }
        if let Ok(arr) = block.get_item(col) {
            return arr.extract();
        }
    }
    Err(PyValueError::new_err(
        r#"atoms block must have an "element" or "symbol" column"#,
    ))
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
