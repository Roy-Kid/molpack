//! Python wrapper for the pack result.
//!
//! [`PyState`] is returned by every engine entry's ``run()``
//! (`GenCanPack`, `CbmcGrow`): the packed ``molrs.Frame`` plus structured
//! diagnostics (`converged` / `fdist` / `frest` / `softened` / `intra`).

use crate::helpers::NpF;
use molpack::F;
use molpack::State;
use numpy::IntoPyArray;
use numpy::PyArray2;
use pyo3::prelude::*;

/// Intra-molecular residual of one configuration: same-copy **scored** vs
/// **exempted** minima, in Å (minimum image). An empty class is ``+inf``.
///
/// Nested on [`PyState`] as ``result.intra.scored`` / ``.exempted``.
/// Forwards the values assembled in Rust; does not recompute from positions.
#[pyclass(name = "IntraResidual", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct PyIntraResidual {
    /// Minimum same-copy pair distance among pairs the table scores (Å).
    #[pyo3(get)]
    pub scored: F,
    /// Minimum same-copy pair distance among pairs the table exempts (Å).
    #[pyo3(get)]
    pub exempted: F,
}

#[pymethods]
impl PyIntraResidual {
    fn __repr__(&self) -> String {
        format!(
            "IntraResidual(scored={}, exempted={})",
            self.scored, self.exempted
        )
    }
}

#[pyclass(name = "State", from_py_object)]
pub struct PyState {
    pub(crate) inner: State,
    /// The Python ``molrs.Frame`` exported once at pack time.
    pub(crate) py_frame: Py<PyAny>,
}

impl Clone for PyState {
    fn clone(&self) -> Self {
        Python::attach(|py| Self {
            inner: self.inner.clone(),
            py_frame: self.py_frame.clone_ref(py),
        })
    }
}

#[pymethods]
impl PyState {
    /// Packed atom positions as a numpy array of shape ``(N, 3)``.
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<NpF>> {
        let pos = self.inner.positions();
        let n = pos.len();
        let flat: Vec<F> = pos.iter().flat_map(|p| [p[0], p[1], p[2]]).collect();
        let arr = ndarray::Array2::from_shape_vec((n, 3), flat).expect("positions shape");
        arr.into_pyarray(py)
    }

    /// Packed ``molrs.Frame`` (same object every access).
    ///
    /// Topology is replayed onto the packed coordinates. A periodic cell is
    /// present only if the engine declared one (``with_periodic_box``);
    /// otherwise ``frame.box`` is ``None`` and the caller assigns one.
    #[getter]
    fn frame(&self, py: Python<'_>) -> Py<PyAny> {
        self.py_frame.clone_ref(py)
    }

    #[getter]
    fn elements(&self) -> Vec<String> {
        let atoms = self.inner.frame.get("atoms").expect("no atoms block");
        atoms
            .get_string("element")
            .expect("no element column")
            .iter()
            .cloned()
            .collect()
    }

    #[getter]
    fn converged(&self) -> bool {
        self.inner.converged
    }

    #[getter]
    fn fdist(&self) -> F {
        self.inner.fdist
    }

    #[getter]
    fn frest(&self) -> F {
        self.inner.frest
    }

    /// How many times the growth solver had to relax its constructive
    /// hard-core guarantee (always 0 on the GENCAN path).
    #[getter]
    fn softened(&self) -> usize {
        self.inner.softened
    }

    /// Same-copy scored vs exempted minima (Å, minimum image).
    ///
    /// Forwards [`molpack::State::intra`]; does not recompute from
    /// positions. An empty class is ``float('inf')``.
    #[getter]
    fn intra(&self) -> PyIntraResidual {
        PyIntraResidual {
            scored: self.inner.intra.scored,
            exempted: self.inner.intra.exempted,
        }
    }

    #[getter]
    fn natoms(&self) -> usize {
        self.inner.natoms()
    }

    fn __repr__(&self) -> String {
        format!(
            "State(converged={}, fdist={:.4}, frest={:.4}, natoms={})",
            self.inner.converged,
            self.inner.fdist,
            self.inner.frest,
            self.inner.natoms()
        )
    }
}
