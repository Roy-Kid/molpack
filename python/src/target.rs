//! Python wrapper for packing `Target`.
//!
//! A [`PyTarget`] describes one type of molecule to pack: its template
//! coordinates, radii, element symbols, and the number of copies.
//!
//! Unlike the old molrs-python binding, this standalone molpack binding does
//! not accept a `molrs.Frame` object directly (PyO3 types do not cross
//! extension modules). Construct targets from numpy arrays via
//! [`PyTarget::from_coords`]; callers with a `molrs.Frame` should extract
//! positions/elements on the Python side and forward them here.

use crate::constraint::extract_restraints;
use crate::helpers::NpF;
use molpack::F;
use molpack::target::Target;
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "Target", from_py_object)]
#[derive(Clone)]
pub struct PyTarget {
    pub(crate) inner: Target,
}

#[pymethods]
impl PyTarget {
    /// Create a target from coordinate and radii arrays.
    ///
    /// Parameters
    /// ----------
    /// positions : numpy.ndarray, shape (N, 3), dtype float64
    ///     Template atom positions.
    /// radii : numpy.ndarray, shape (N,), dtype float64
    ///     Van der Waals radii.
    /// count : int
    ///     Number of copies to pack.
    /// elements : list[str], optional
    ///     Element symbols, one per atom. Defaults to ``["X"] * N``.
    ///
    /// Returns
    /// -------
    /// Target
    #[staticmethod]
    #[pyo3(signature = (positions, radii, count, elements=None))]
    fn from_coords(
        positions: PyReadonlyArray2<'_, NpF>,
        radii: PyReadonlyArray1<'_, NpF>,
        count: usize,
        elements: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let pos = positions.as_array();
        let rad = radii.as_array();
        if pos.ncols() != 3 {
            return Err(PyValueError::new_err("positions must have shape (N, 3)"));
        }
        if pos.nrows() != rad.len() {
            return Err(PyValueError::new_err(
                "positions and radii must have the same length",
            ));
        }
        if let Some(ref els) = elements
            && els.len() != pos.nrows()
        {
            return Err(PyValueError::new_err(
                "elements must have the same length as positions",
            ));
        }
        let coords: Vec<[F; 3]> = pos.rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect();
        let radii_vec: Vec<F> = rad.to_vec();
        let mut target = Target::from_coords(&coords, &radii_vec, count);
        if let Some(els) = elements {
            target.elements = els;
        }
        Ok(PyTarget { inner: target })
    }

    fn with_name(&self, name: &str) -> Self {
        PyTarget {
            inner: self.inner.clone().with_name(name),
        }
    }

    fn with_constraint(&self, constraint: &Bound<'_, pyo3::types::PyAny>) -> PyResult<Self> {
        let restraints = extract_restraints(constraint)?;
        let mut target = self.inner.clone();
        for r in restraints {
            target = target.with_restraint(r);
        }
        Ok(PyTarget { inner: target })
    }

    fn with_constraint_for_atoms(
        &self,
        indices: Vec<usize>,
        constraint: &Bound<'_, pyo3::types::PyAny>,
    ) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        let restraints = extract_restraints(constraint)?;
        let mut target = self.inner.clone();
        for r in restraints {
            target = target.with_restraint_for_atoms(&indices, r);
        }
        Ok(PyTarget { inner: target })
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
