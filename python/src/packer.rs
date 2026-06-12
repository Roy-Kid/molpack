//! Python wrapper for the Molpack molecular packer.
//!
//! [`PyPacker`] is the Python face of [`molpack::Molpack`]. Builders mirror
//! the Rust names 1:1. `Molpack.pack()` returns a ready-to-use `molrs.Frame`;
//! [`PyPackResult`] remains available through `pack_with_report()` for callers
//! that need structured diagnostics.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use crate::constraint::extract_restraint;
use crate::handler::PyHandlerWrapper;
use crate::helpers::{NpF, pack_error_to_pyerr, take_err};
use crate::parallel::rayon_compiled;
use crate::target::PyTarget;
use molpack::F;
use molpack::handler::MolpackLogLevel;
use molpack::packer::{Molpack, PackResult};
use numpy::IntoPyArray;
use numpy::PyArray2;
use pyo3::prelude::*;
// ── PackResult ─────────────────────────────────────────────────────────────

#[pyclass(name = "PackResult", from_py_object)]
#[derive(Clone)]
pub struct PyPackResult {
    inner: PackResult,
    /// The periodic box `(min, max)` the packer was configured with, if any.
    /// Stamped onto `frame.box` so callers don't re-create it. `None` when no
    /// box was declared (e.g. restraint-only packing).
    box_bounds: Option<([F; 3], [F; 3])>,
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

    /// Packed result as a ready-to-use ``molrs.Frame``.
    ///
    /// The frame — full topology replayed onto the packed coordinates, plus the
    /// periodic box if one was declared via :meth:`Molpack.with_periodic_box` —
    /// is assembled by the Rust core (:mod:`molpack.assemble`), so every
    /// language binding returns an identical result. This getter only marshals
    /// that frame across the language boundary. Force fields are out of scope —
    /// merge them separately.
    #[getter]
    fn frame<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        crate::frame_marshal::rust_frame_to_py(py, &self.inner.frame, self.box_bounds)
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
    pub(crate) log_level: Option<MolpackLogLevel>,
    pub(crate) log_frequency: Option<usize>,
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
            log_level: None,
            log_frequency: None,
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

    /// Run the pair-kernel reductions on rayon worker threads.
    ///
    /// Fail-fast: enabling parallel evaluation in a wheel that was **not**
    /// built with the `rayon` feature raises ``RuntimeError`` rather than
    /// silently running serially — a silent no-op here is what makes
    /// thread-count scaling studies show a flat curve. Build with
    /// ``maturin develop --release`` (the wheel ships `rayon` by default).
    ///
    /// Parallelism also only applies to global evaluations (`!move_flag`);
    /// per-molecule GENCAN move evaluations stay serial by design.
    fn with_parallel_eval(&self, enabled: bool) -> PyResult<Self> {
        if enabled && !rayon_compiled() {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "parallel evaluation requested but this molpack wheel was built \
                 without the `rayon` feature; rebuild with `maturin develop --release` \
                 (rayon is enabled by default) or check `molpack.rayon_enabled()`",
            ));
        }
        let mut cloned = self.clone_fields();
        cloned.parallel_eval = Some(enabled);
        Ok(cloned)
    }

    fn with_progress(&self, enabled: bool) -> Self {
        let mut cloned = self.clone_fields();
        cloned.progress = enabled;
        cloned
    }

    fn with_lammps_output(&self, enabled: bool) -> Self {
        let mut cloned = self.clone_fields();
        cloned.log_level = Some(if enabled {
            MolpackLogLevel::Progress
        } else {
            MolpackLogLevel::Quiet
        });
        cloned
    }

    fn with_log_level(&self, level: &str) -> PyResult<Self> {
        let parsed = match level.to_ascii_lowercase().as_str() {
            "quiet" | "off" | "none" => MolpackLogLevel::Quiet,
            "summary" => MolpackLogLevel::Summary,
            "progress" | "thermo" => MolpackLogLevel::Progress,
            "verbose" | "debug" => MolpackLogLevel::Verbose,
            other => {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "unknown log level {other:?}; expected quiet, summary, progress, or verbose"
                )));
            }
        };
        let mut cloned = self.clone_fields();
        cloned.log_level = Some(parsed);
        Ok(cloned)
    }

    fn with_log_frequency(&self, n: usize) -> Self {
        let mut cloned = self.clone_fields();
        cloned.log_frequency = Some(n.max(1));
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
    fn pack<'py>(
        &self,
        py: Python<'py>,
        targets: Vec<PyTarget>,
        max_loops: usize,
    ) -> PyResult<Bound<'py, PyAny>> {
        let result = self.pack_report_inner(py, targets, max_loops)?;
        crate::frame_marshal::rust_frame_to_py(py, &result.inner.frame, self.periodic_box)
    }

    #[pyo3(signature = (targets, max_loops=200))]
    fn pack_with_report(
        &self,
        py: Python<'_>,
        targets: Vec<PyTarget>,
        max_loops: usize,
    ) -> PyResult<PyPackResult> {
        self.pack_report_inner(py, targets, max_loops)
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
    fn pack_report_inner(
        &self,
        py: Python<'_>,
        targets: Vec<PyTarget>,
        max_loops: usize,
    ) -> PyResult<PyPackResult> {
        // The Rust core retains each target's frame and assembles the result —
        // the binding just unwraps to the core `Target`.
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
        let level = self.log_level.unwrap_or(if self.progress {
            MolpackLogLevel::Progress
        } else {
            MolpackLogLevel::Quiet
        });
        packer = packer.with_log_level(level);
        if let Some(every) = self.log_frequency {
            packer = packer.with_log_frequency(every);
        }
        if let Some((min, max)) = self.periodic_box {
            packer = packer.with_periodic_box(min, max);
        }

        for gr in &self.global_restraints {
            let shared = extract_restraint(gr.bind(py))?;
            packer = packer.with_global_restraint(shared);
        }

        // Any handler stashing an error or returning `True` from `on_step`
        // flips this shared flag, which terminates the loop at the next check.
        let stop_flag = Arc::new(AtomicBool::new(false));
        for py_h in &self.py_handlers {
            let wrapper = PyHandlerWrapper::new(py_h.clone_ref(py), Arc::clone(&stop_flag));
            packer = packer.with_handler(wrapper);
        }

        let result = packer.pack_with_report(&rust_targets, max_loops);

        // Python exceptions take priority: any PackError from this run is
        // a downstream symptom of the callback that raised.
        if let Some(py_err) = take_err() {
            return Err(py_err);
        }

        let result = result.map_err(pack_error_to_pyerr)?;
        Ok(PyPackResult {
            inner: result,
            // The declared periodic box (if any) is stamped onto the result
            // frame's `box` so callers don't re-create it.
            box_bounds: self.periodic_box,
        })
    }

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
            log_level: self.log_level,
            log_frequency: self.log_frequency,
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
