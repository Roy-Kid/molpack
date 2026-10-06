//! Python wrapper for the pack result.
//!
//! [`PyState`] is returned by every engine entry's ``run()``
//! (`GenCanPack`, `CbmcGrow`): the packed ``molrs.store.Frame`` plus structured
//! diagnostics (`converged` / `fdist` / `frest` / `degraded` / `intra`).

use molpack::State;
use molrs::op::types::F;
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
    /// The Python ``molrs.store.Frame`` exported once at pack time.
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
    fn positions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<F>> {
        self.inner
            .frame
            .coords()
            .expect("an assembled frame always carries x / y / z")
            .into_pyarray(py)
    }

    /// Packed ``molrs.store.Frame`` (same object every access).
    ///
    /// Topology is replayed onto the packed coordinates. A periodic cell is
    /// present only if the engine declared one (``with_periodic_box``);
    /// otherwise ``frame.box`` is ``None`` and the caller assigns one.
    #[getter]
    fn frame(&self, py: Python<'_>) -> Py<PyAny> {
        self.py_frame.clone_ref(py)
    }

    /// Element symbol per atom. Raises ``KeyError`` when the templates carry
    /// no ``element`` column (a coarse-grained model, for one).
    #[getter]
    fn elements(&self) -> PyResult<Vec<String>> {
        self.inner
            .frame
            .get("atoms")
            .and_then(|atoms| atoms.get(molrs::store::keys::ELEMENT))
            .and_then(molrs::store::Column::as_string)
            .map(|column| column.iter().cloned().collect())
            .ok_or_else(|| {
                pyo3::exceptions::PyKeyError::new_err("the packed frame has no element column")
            })
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
    fn degraded(&self) -> usize {
        self.inner.degraded
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
