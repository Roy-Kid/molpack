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
        .map_err(|_| {
            PyValueError::new_err(format!(r#"atoms block has no "{name}" column"#))
        })?
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
    /// name : str
    ///     Molecule name (used in diagnostics and output).
    /// frame : Frame-like
    ///     Any object with ``frame["atoms"]`` returning a block.
    ///     Supports ``molrs.Frame``, ``molpy.Frame``, and plain dicts.
    /// count : int
    ///     Number of copies to pack.
    #[new]
    #[pyo3(signature = (name, frame, count))]
    fn new(
        name: &str,
        frame: &Bound<'_, PyAny>,
        count: usize,
    ) -> PyResult<Self> {
        let atoms = frame.get_item("atoms").map_err(|_| {
            PyValueError::new_err(r#"frame must have an "atoms" block"#)
        })?;

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
        let target = target.with_name(name);

        Ok(PyTarget { inner: target })
    }

    fn with_name(&self, name: &str) -> Self {
        PyTarget {
            inner: self.inner.clone().with_name(name),
        }
    }

    fn with_restraint(&self, constraint: &Bound<'_, pyo3::types::PyAny>) -> PyResult<Self> {
        let restraint = extract_restraint(constraint)?;
        Ok(PyTarget {
            inner: self.inner.clone().with_restraint(restraint),
        })
    }

    fn with_restraint_for_atoms(
        &self,
        indices: Vec<usize>,
        constraint: &Bound<'_, pyo3::types::PyAny>,
    ) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        let restraint = extract_restraint(constraint)?;
        Ok(PyTarget {
            inner: self
                .inner
                .clone()
                .with_restraint_for_atoms(&indices, restraint),
        })
    }

    fn with_maxmove(&self, maxmove: usize) -> Self {
        PyTarget {
            inner: self.inner.clone().with_maxmove(maxmove),
        }
    }

    fn with_center(&self) -> Self {
        PyTarget {
            inner: self.inner.clone().with_center(),
        }
    }

    fn without_centering(&self) -> Self {
        PyTarget {
            inner: self.inner.clone().without_centering(),
        }
    }

    fn constrain_rotation_x(&self, center_deg: NpF, half_width_deg: NpF) -> Self {
        PyTarget {
            inner: self
                .inner
                .clone()
                .constrain_rotation_x(center_deg, half_width_deg),
        }
    }

    fn constrain_rotation_y(&self, center_deg: NpF, half_width_deg: NpF) -> Self {
        PyTarget {
            inner: self
                .inner
                .clone()
                .constrain_rotation_y(center_deg, half_width_deg),
        }
    }

    fn constrain_rotation_z(&self, center_deg: NpF, half_width_deg: NpF) -> Self {
        PyTarget {
            inner: self
                .inner
                .clone()
                .constrain_rotation_z(center_deg, half_width_deg),
        }
    }

    fn fixed_at(&self, position: [NpF; 3]) -> Self {
        PyTarget {
            inner: self.inner.clone().fixed_at(position),
        }
    }

    fn fixed_at_with_euler(&self, position: [NpF; 3], euler: [NpF; 3]) -> Self {
        PyTarget {
            inner: self.inner.clone().fixed_at_with_euler(position, euler),
        }
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
        if index == 0 || index > natoms {
            return Err(PyValueError::new_err(format!(
                "atom indices must use Packmol 1-based indexing in 1..={natoms}, got {index}",
            )));
        }
    }
    Ok(())
}
