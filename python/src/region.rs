//! Python wrapper for [`molpack::StlRegion`].
//!
//! Construction is `from_file(path, scale=1.0)` with `scale` in Å per file
//! unit. The object is a Region; `Target.with_restraint` lifts it through
//! [`molpack::RegionRestraint`].
//!
//! `contains` and `signed_distance` answer the region's two questions for a
//! batch of points, so a caller can ask what the packer was told to enforce
//! without re-deriving mesh geometry on the Python side.

use std::path::PathBuf;

use molpack::{F, Region, StlError, StlRegion};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;

/// Closed triangle mesh as a packing region (file units × `scale` → Å).
#[pyclass(name = "StlRegion", from_py_object)]
#[derive(Clone)]
pub struct PyStlRegion {
    pub(crate) inner: StlRegion,
}

fn stl_error_to_pyerr(err: StlError) -> PyErr {
    match err {
        StlError::Io { .. } => PyOSError::new_err(err.to_string()),
        other => PyValueError::new_err(other.to_string()),
    }
}

#[pymethods]
impl PyStlRegion {
    /// Load a watertight STL. `scale=1.0` means the file is already in Å.
    #[classmethod]
    #[pyo3(signature = (path, scale=1.0))]
    fn from_file(
        _cls: &Bound<'_, pyo3::types::PyType>,
        path: PathBuf,
        scale: f64,
    ) -> PyResult<Self> {
        // `from_file` gates `scale` itself and `InvalidScale` already maps to
        // exactly this `ValueError`; re-checking here would be a second home
        // for one rule.
        let inner = StlRegion::from_file(&path, scale).map_err(stl_error_to_pyerr)?;
        Ok(Self { inner })
    }

    /// Whether each point of an ``(n, 3)`` array is inside the mesh.
    fn contains<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, F>,
    ) -> PyResult<Bound<'py, PyArray1<bool>>> {
        let inside: Vec<bool> = rows(&points)?
            .into_iter()
            .map(|x| self.inner.contains(&x))
            .collect();
        Ok(inside.into_pyarray(py))
    }

    /// Signed distance to the mesh for each point: negative inside, Å.
    fn signed_distance<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, F>,
    ) -> PyResult<Bound<'py, PyArray1<F>>> {
        let d: Vec<F> = rows(&points)?
            .into_iter()
            .map(|x| self.inner.signed_distance(&x))
            .collect();
        Ok(d.into_pyarray(py))
    }

    fn __repr__(&self) -> &'static str {
        "StlRegion"
    }
}

/// `(n, 3)` array to owned points; anything else is a named `ValueError`.
fn rows(points: &PyReadonlyArray2<'_, F>) -> PyResult<Vec<[F; 3]>> {
    let view = points.as_array();
    if view.ncols() != 3 {
        return Err(PyValueError::new_err(format!(
            "points must be (n, 3), got (_, {})",
            view.ncols()
        )));
    }
    Ok(view
        .rows()
        .into_iter()
        .map(|r| [r[0], r[1], r[2]])
        .collect())
}
