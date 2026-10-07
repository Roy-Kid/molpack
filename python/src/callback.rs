//! Python-defined [`Callback`] hooks.
//!
//! A Python object attached via an entry's ``with_callback`` may
//! implement any subset of three optional methods:
//!
//! ```python
//! class MyHook:
//!     def on_start(self, ntotat: int, ntotmol: int) -> None: ...
//!     def on_step(self, step: StepReport, ctx: StepContext) -> bool | None: ...  # True → stop
//!     def on_finish(self) -> None: ...
//! ```
//!
//! `on_step` mirrors the Rust trait's `(step, sys)` pair: `ctx` is a borrow
//! guard over the live packing context, valid only for the duration of the
//! callback (accessing it afterwards raises `RuntimeError`). Its
//! `positions` property materialises an owned `(ntotat, 3)` float64 NumPy
//! array on demand — callbacks that never touch it pay nothing.
//!
//! Missing methods are silently skipped (matching the Rust trait's default
//! no-op impls). Exceptions raised inside any method are stashed in
//! [`errors::PACK_ERR`][crate::errors] and trigger early termination;
//! the entry's `run()` re-raises after the loop exits.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::errors::stash_err;
use molpack::PackContext;
use molpack::{Callback, StepReport};
use molrs::op::F;
use numpy::IntoPyArray;
use pyo3::prelude::*;
use pyo3::types::PyAny;

// ============================================================================
// PyStageProgress — which stage of the run a callback came from. Kept nested
// (unlike `PhaseProgress`, flattened below) because a pipeline's stage identity
// reads as one thing: `step.stage.index` / `.total` / `.name`.
// ============================================================================

/// Identifies the stage a callback comes from. A single-stage run reports
/// ``index == 0`` and ``total == 1``.
#[pyclass(name = "StageProgress", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct PyStageProgress {
    /// 0-based index of this stage in the run.
    #[pyo3(get)]
    pub index: usize,
    /// How many stages the run has.
    #[pyo3(get)]
    pub total: usize,
    /// The stage's own name.
    #[pyo3(get)]
    pub name: String,
}

#[pymethods]
impl PyStageProgress {
    fn __repr__(&self) -> String {
        format!(
            "StageProgress({}, {}/{})",
            self.name,
            self.index + 1,
            self.total
        )
    }
}

// ============================================================================
// PyStepReport — read-only snapshot passed to `on_step`. Flattens Rust's
// nested `PhaseProgress` for Python ergonomics.
// ============================================================================

#[pyclass(name = "StepReport", frozen)]
pub struct PyStepReport {
    /// Which stage of the run this callback came from.
    #[pyo3(get)]
    pub stage: PyStageProgress,
    /// 0-based outer-loop iteration within the current phase.
    #[pyo3(get)]
    pub loop_idx: usize,
    /// Maximum outer loops budgeted for this phase.
    #[pyo3(get)]
    pub max_loops: usize,
    /// 0-based phase index.
    #[pyo3(get)]
    pub phase: usize,
    /// Total number of phases (ntype + 1).
    #[pyo3(get)]
    pub total_phases: usize,
    /// `Some(itype)` for a per-type compaction phase; `None` for the
    /// final all-types phase.
    #[pyo3(get)]
    pub molecule_type: Option<usize>,
    /// Max inter-molecular overlap violation (0.0 = no overlap).
    #[pyo3(get)]
    pub fdist: F,
    /// Max restraint violation (0.0 = all restraints satisfied).
    #[pyo3(get)]
    pub frest: F,
    /// GENCAN objective at the user's radii (Packmol's ``fx``); growth
    /// reports ``0.0``.
    #[pyo3(get)]
    pub f: F,
    /// Improvement from last iteration, as percentage (positive = improving).
    #[pyo3(get)]
    pub improvement_pct: F,
    /// Current radius scaling factor (starts at `discale`, decays to 1.0).
    #[pyo3(get)]
    pub radscale: F,
    /// Convergence precision target.
    #[pyo3(get)]
    pub precision: F,
}

impl PyStepReport {
    fn from_report(step: &StepReport) -> Self {
        Self {
            stage: PyStageProgress {
                index: step.stage.index,
                total: step.stage.total,
                name: step.stage.name.to_owned(),
            },
            loop_idx: step.loop_idx,
            max_loops: step.max_loops,
            phase: step.phase.phase,
            total_phases: step.phase.total_phases,
            molecule_type: step.phase.molecule_type,
            fdist: step.fdist,
            frest: step.frest,
            f: step.f,
            improvement_pct: step.improvement_pct,
            radscale: step.radscale,
            precision: step.precision,
        }
    }
}

#[pymethods]
impl PyStepReport {
    fn __repr__(&self) -> String {
        format!(
            "StepReport(stage={} {}/{}, phase={}/{}, loop={}/{}, \
             fdist={:.3e}, frest={:.3e}, improvement={:.2}%)",
            self.stage.name,
            self.stage.index + 1,
            self.stage.total,
            self.phase + 1,
            self.total_phases,
            self.loop_idx + 1,
            self.max_loops,
            self.fdist,
            self.frest,
            self.improvement_pct,
        )
    }
}

// PythonCallback — bridges a Python object to the Rust `Callback` trait.
//
// Optional methods are resolved once at construction; missing ones stay
// `None` and are silently skipped. `stop_flag` lives behind an atomic
// because `Callback::should_stop` is `&self` while writes happen via the
// mutating `on_*` methods.

// PyStepContext — borrow guard over the live `PackContext`, handed to
// `on_step` and invalidated the moment the callback returns (spec: FFI
// stale-handle invalidation). Data properties copy on access, so a stashed
// guard can never dangle — it only errors.

#[pyclass(name = "StepContext", unsendable)]
pub struct PyStepContext {
    sys: std::cell::Cell<Option<*const PackContext>>,
}

impl PyStepContext {
    fn expired() -> pyo3::PyErr {
        pyo3::exceptions::PyRuntimeError::new_err(
            "StepContext expired: it is only valid inside the on_step callback \
             it was passed to — copy what you need (e.g. ctx.positions) there",
        )
    }

    fn live(&self) -> PyResult<&PackContext> {
        match self.sys.get() {
            // SAFETY: the pointer is set right before the callback and
            // cleared right after it returns, on the same thread; while it
            // is Some the borrow in `Callback::on_step` is still alive.
            Some(p) => Ok(unsafe { &*p }),
            None => Err(Self::expired()),
        }
    }
}

#[pymethods]
impl PyStepContext {
    /// Owned ``(ntotat, 3)`` float64 array of the live coordinates.
    /// Atoms the growth solver has not placed yet sit at their sentinel.
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
        Ok(ndarray::Array2::from(self.live()?.xcart.clone()).into_pyarray(py))
    }

    /// Total atom count of the packing system.
    #[getter]
    fn natoms(&self) -> PyResult<usize> {
        Ok(self.live()?.xcart.len())
    }

    fn __repr__(&self) -> String {
        match self.sys.get() {
            Some(_) => format!(
                "StepContext(natoms={})",
                self.live().map(|s| s.xcart.len()).unwrap_or(0)
            ),
            None => "StepContext(expired)".to_owned(),
        }
    }
}

pub(crate) struct PythonCallback {
    on_start: Option<Py<PyAny>>,
    on_step: Option<Py<PyAny>>,
    on_finish: Option<Py<PyAny>>,
    stop_flag: Arc<AtomicBool>,
}

impl PythonCallback {
    pub(crate) fn new(obj: Py<PyAny>, stop_flag: Arc<AtomicBool>) -> Self {
        Python::attach(|py| {
            let bound = obj.bind(py);
            Self {
                on_start: bound.getattr("on_start").ok().map(|m| m.unbind()),
                on_step: bound.getattr("on_step").ok().map(|m| m.unbind()),
                on_finish: bound.getattr("on_finish").ok().map(|m| m.unbind()),
                stop_flag,
            }
        })
    }

    fn fail(&self, err: PyErr) {
        stash_err(err);
        self.stop_flag.store(true, Ordering::SeqCst);
    }
}

impl Callback for PythonCallback {
    fn on_start(&mut self, ntotat: usize, ntotmol: usize) {
        let Some(m) = &self.on_start else { return };
        Python::attach(|py| {
            if let Err(e) = m.bind(py).call1((ntotat, ntotmol)) {
                self.fail(e);
            }
        });
    }

    fn on_step(&mut self, step: &StepReport, sys: &PackContext) {
        let Some(m) = &self.on_step else { return };
        Python::attach(|py| {
            let py_step = match Py::new(py, PyStepReport::from_report(step)) {
                Ok(v) => v,
                Err(e) => {
                    self.fail(e);
                    return;
                }
            };
            let guard = match Py::new(
                py,
                PyStepContext {
                    sys: std::cell::Cell::new(Some(sys as *const PackContext)),
                },
            ) {
                Ok(v) => v,
                Err(e) => {
                    self.fail(e);
                    return;
                }
            };
            let ret = m.bind(py).call1((py_step, guard.clone_ref(py)));
            // Stale-handle invalidation: the borrow of `sys` ends here.
            guard.borrow(py).sys.set(None);
            match ret {
                // `True` requests early stop; other return values continue.
                Ok(ret) => {
                    if let Ok(true) = ret.extract::<bool>() {
                        self.stop_flag.store(true, Ordering::SeqCst);
                    }
                }
                Err(e) => self.fail(e),
            }
        });
    }

    fn on_finish(&mut self, _sys: &PackContext) {
        let Some(m) = &self.on_finish else { return };
        Python::attach(|py| {
            if let Err(e) = m.bind(py).call1(()) {
                self.fail(e);
            }
        });
    }

    fn should_stop(&self) -> bool {
        self.stop_flag.load(Ordering::SeqCst)
    }
}
