//! Python wrapper for the Molpack molecular packer.
//!
//! [`PyPacker`] is the Python face of [`molpack::Molpack`]. Builders mirror
//! the Rust names 1:1. [`PyPackResult`] is the output container — callers
//! get a Frame-compatible `.frame` (dict with an `atoms` block) and a
//! numpy `positions` view.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use crate::constraint::extract_restraint;
use crate::handler::PyHandlerWrapper;
use crate::helpers::{NpF, pack_error_to_pyerr, take_err};
use crate::target::PyTarget;
use molpack::F;
use molpack::handler::ProgressHandler;
use molpack::packer::{Molpack, PackResult};
use numpy::IntoPyArray;
use numpy::PyArray2;
use pyo3::prelude::*;
use pyo3::types::PyDict;

// ── PackResult ─────────────────────────────────────────────────────────────

#[pyclass(name = "PackResult", from_py_object)]
#[derive(Clone)]
pub struct PyPackResult {
    inner: PackResult,
}

#[pymethods]
impl PyPackResult {
    /// Packed atom positions as a numpy array of shape ``(N, 3)``.
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<NpF>> {
        let pos = self.inner.positions();
        let n = pos.len();
        let flat: Vec<F> = pos.iter().flat_map(|p| [p[0], p[1], p[2]]).collect();
        let arr = ndarray::Array2::from_shape_vec((n, 3), flat).expect("positions shape");
        arr.into_pyarray(py)
    }

    /// Packed frame as a dict compatible with ``molrs.Frame``.
    ///
    /// Contains a single ``"atoms"`` block with at least the columns
    /// ``x``, ``y``, ``z``, ``element``. Users persist the result by
    /// passing this frame to a writer — e.g. ``molrs.write_pdb``. molpack
    /// does not provide writers.
    #[getter]
    fn frame<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let atoms_block = self
            .inner
            .frame
            .get("atoms")
            .expect("frame has no 'atoms' block");

        let x = atoms_block.get_float("x").expect("no 'x' column");
        let y = atoms_block.get_float("y").expect("no 'y' column");
        let z = atoms_block.get_float("z").expect("no 'z' column");
        let elements: Vec<String> = atoms_block
            .get_string("element")
            .expect("no 'element' column")
            .iter()
            .cloned()
            .collect();

        let atoms = PyDict::new(py);
        atoms.set_item(
            "x",
            ndarray::Array1::from_iter(x.iter().copied()).into_pyarray(py),
        )?;
        atoms.set_item(
            "y",
            ndarray::Array1::from_iter(y.iter().copied()).into_pyarray(py),
        )?;
        atoms.set_item(
            "z",
            ndarray::Array1::from_iter(z.iter().copied()).into_pyarray(py),
        )?;
        atoms.set_item("element", elements)?;

        let frame = PyDict::new(py);
        frame.set_item("atoms", atoms)?;
        Ok(frame)
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

    #[getter]
    fn natoms(&self) -> usize {
        self.inner.natoms()
    }

    fn __repr__(&self) -> String {
        format!(
            "PackResult(converged={}, fdist={:.4}, frest={:.4}, natoms={})",
            self.inner.converged,
            self.inner.fdist,
            self.inner.frest,
            self.inner.natoms()
        )
    }
}

// ── Molpack ────────────────────────────────────────────────────────────────

/// Molecular packer — configured via the fluent `with_*` API.
///
/// ```python
/// from molpack import Molpack
/// packer = (
///     Molpack()
///     .with_tolerance(2.0)
///     .with_precision(0.01)
///     .with_seed(42)
/// )
/// result = packer.pack(targets, max_loops=200)
/// ```
#[pyclass(name = "Molpack")]
pub struct PyPacker {
    pub(crate) precision: Option<F>,
    pub(crate) tolerance: Option<F>,
    pub(crate) inner_iterations: Option<usize>,
    pub(crate) init_passes: Option<usize>,
    pub(crate) init_box_half_size: Option<F>,
    pub(crate) perturb_fraction: Option<F>,
    pub(crate) random_perturb: Option<bool>,
    pub(crate) perturb: Option<bool>,
    pub(crate) seed: Option<u64>,
    pub(crate) parallel_eval: Option<bool>,
    pub(crate) progress: bool,
    /// Global periodic-boundary box — `(min, max)` in Å. Every axis is
    /// treated as periodic.
    pub(crate) periodic_box: Option<([F; 3], [F; 3])>,
    pub(crate) py_handlers: Vec<Py<pyo3::types::PyAny>>,
    pub(crate) global_restraints: Vec<Py<pyo3::types::PyAny>>,
}

impl Default for PyPacker {
    fn default() -> Self {
        PyPacker {
            precision: None,
            tolerance: None,
            inner_iterations: None,
            init_passes: None,
            init_box_half_size: None,
            perturb_fraction: None,
            random_perturb: None,
            perturb: None,
            seed: None,
            parallel_eval: None,
            progress: false,
            periodic_box: None,
            py_handlers: Vec::new(),
            global_restraints: Vec::new(),
        }
    }
}

#[pymethods]
impl PyPacker {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    fn with_tolerance(&self, tolerance: NpF) -> Self {
        let mut cloned = self.clone_fields();
        cloned.tolerance = Some(tolerance);
        cloned
    }

    fn with_precision(&self, precision: NpF) -> Self {
        let mut cloned = self.clone_fields();
        cloned.precision = Some(precision);
        cloned
    }

    fn with_inner_iterations(&self, n: usize) -> Self {
        let mut cloned = self.clone_fields();
        cloned.inner_iterations = Some(n);
        cloned
    }

    fn with_init_passes(&self, n: usize) -> Self {
        let mut cloned = self.clone_fields();
        cloned.init_passes = Some(n);
        cloned
    }

    fn with_init_box_half_size(&self, half_size: NpF) -> Self {
        let mut cloned = self.clone_fields();
        cloned.init_box_half_size = Some(half_size);
        cloned
    }

    /// Declare a global periodic-boundary box (Packmol ``pbc``). Every
    /// axis is periodic. See :meth:`Molpack.with_periodic_box` in the
    /// Rust API for the full semantics.
    fn with_periodic_box(&self, min: [NpF; 3], max: [NpF; 3]) -> Self {
        let mut cloned = self.clone_fields();
        cloned.periodic_box = Some((min, max));
        cloned
    }

    fn with_perturb_fraction(&self, f: NpF) -> Self {
        let mut cloned = self.clone_fields();
        cloned.perturb_fraction = Some(f);
        cloned
    }

    fn with_random_perturb(&self, enabled: bool) -> Self {
        let mut cloned = self.clone_fields();
        cloned.random_perturb = Some(enabled);
        cloned
    }

    fn with_perturb(&self, enabled: bool) -> Self {
        let mut cloned = self.clone_fields();
        cloned.perturb = Some(enabled);
        cloned
    }

    fn with_seed(&self, seed: u64) -> Self {
        let mut cloned = self.clone_fields();
        cloned.seed = Some(seed);
        cloned
    }

    fn with_parallel_eval(&self, enabled: bool) -> Self {
        let mut cloned = self.clone_fields();
        cloned.parallel_eval = Some(enabled);
        cloned
    }

    fn with_progress(&self, enabled: bool) -> Self {
        let mut cloned = self.clone_fields();
        cloned.progress = enabled;
        cloned
    }

    /// Append a Python handler. See :class:`StepInfo` for the callback
    /// contract and :mod:`molpack` for the ``Handler`` protocol.
    fn with_handler(&self, handler: Py<pyo3::types::PyAny>) -> Self {
        let mut cloned = self.clone_fields();
        cloned.py_handlers.push(handler);
        cloned
    }

    /// Append a global restraint — applied to every atom of every target
    /// at ``pack()`` time. Equivalent to calling
    /// ``target.with_restraint(r)`` on every target.
    fn with_global_restraint(&self, restraint: Py<pyo3::types::PyAny>) -> Self {
        let mut cloned = self.clone_fields();
        cloned.global_restraints.push(restraint);
        cloned
    }

    #[pyo3(signature = (targets, max_loops=200))]
    fn pack(
        &self,
        py: Python<'_>,
        targets: Vec<PyTarget>,
        max_loops: usize,
    ) -> PyResult<PyPackResult> {
        let rust_targets: Vec<_> = targets.into_iter().map(|t| t.inner).collect();

        let mut packer = Molpack::new();
        if let Some(v) = self.tolerance {
            packer = packer.with_tolerance(v);
        }
        if let Some(v) = self.precision {
            packer = packer.with_precision(v);
        }
        if let Some(v) = self.inner_iterations {
            packer = packer.with_inner_iterations(v);
        }
        if let Some(v) = self.init_passes {
            packer = packer.with_init_passes(v);
        }
        if let Some(v) = self.init_box_half_size {
            packer = packer.with_init_box_half_size(v);
        }
        if let Some(v) = self.perturb_fraction {
            packer = packer.with_perturb_fraction(v);
        }
        if let Some(v) = self.random_perturb {
            packer = packer.with_random_perturb(v);
        }
        if let Some(v) = self.perturb {
            packer = packer.with_perturb(v);
        }
        if let Some(v) = self.seed {
            packer = packer.with_seed(v);
        }
        if let Some(v) = self.parallel_eval {
            packer = packer.with_parallel_eval(v);
        }
        if let Some((min, max)) = self.periodic_box {
            packer = packer.with_periodic_box(min, max);
        }

        for gr in &self.global_restraints {
            let shared = extract_restraint(gr.bind(py))?;
            packer = packer.with_global_restraint(shared);
        }

        if self.progress {
            packer = packer.with_handler(ProgressHandler::new());
        }

        // Any handler stashing an error or returning `True` from `on_step`
        // flips this shared flag, which terminates the loop at the next check.
        let stop_flag = Arc::new(AtomicBool::new(false));
        for py_h in &self.py_handlers {
            let wrapper = PyHandlerWrapper::new(py_h.clone_ref(py), Arc::clone(&stop_flag));
            packer = packer.with_handler(wrapper);
        }

        let result = packer.pack(&rust_targets, max_loops);

        // Python exceptions take priority: any PackError from this run is
        // a downstream symptom of the callback that raised.
        if let Some(py_err) = take_err() {
            return Err(py_err);
        }

        let result = result.map_err(pack_error_to_pyerr)?;
        Ok(PyPackResult { inner: result })
    }

    fn __repr__(&self) -> String {
        format!(
            "Molpack(handlers={}, global_restraints={}, seed={:?})",
            self.py_handlers.len(),
            self.global_restraints.len(),
            self.seed,
        )
    }
}

impl PyPacker {
    fn clone_fields(&self) -> PyPacker {
        // `Py<PyAny>::clone_ref` needs a GIL token; the builder is
        // called from Python so `Python::attach` is a cheap no-op here.
        Python::attach(|py| PyPacker {
            precision: self.precision,
            tolerance: self.tolerance,
            inner_iterations: self.inner_iterations,
            init_passes: self.init_passes,
            init_box_half_size: self.init_box_half_size,
            perturb_fraction: self.perturb_fraction,
            random_perturb: self.random_perturb,
            perturb: self.perturb,
            seed: self.seed,
            parallel_eval: self.parallel_eval,
            progress: self.progress,
            periodic_box: self.periodic_box,
            py_handlers: self.py_handlers.iter().map(|h| h.clone_ref(py)).collect(),
            global_restraints: self
                .global_restraints
                .iter()
                .map(|r| r.clone_ref(py))
                .collect(),
        })
    }
}
