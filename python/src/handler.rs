//! Python-defined [`Handler`] hooks.
//!
//! A Python object attached via [`crate::packer::PyPacker::add_handler`] may
//! implement any subset of three optional methods:
//!
//! ```python
//! class MyHook:
//!     def on_start(self, ntotat: int, ntotmol: int) -> None: ...
//!     def on_step(self, info: StepInfo) -> bool | None: ...   # True → stop
//!     def on_finish(self) -> None: ...
//! ```
//!
//! Missing methods are silently skipped (matching the Rust trait's default
//! no-op impls). Exceptions raised inside any method are stashed in
//! [`helpers::PACK_ERR`][crate::helpers] and trigger early termination;
//! `Molpack.pack()` re-raises after the loop exits.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::helpers::stash_err;
use molpack::F;
use molpack::context::PackContext;
use molpack::handler::{Handler, StepInfo};
use pyo3::prelude::*;
use pyo3::types::PyAny;

// ============================================================================
// PyStepInfo — read-only snapshot passed to `on_step`. Flattens Rust's
// nested `PhaseInfo` for Python ergonomics.
// ============================================================================

#[pyclass(name = "StepInfo", frozen)]
pub struct PyStepInfo {
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
    /// Improvement from last iteration, as percentage (positive = improving).
    #[pyo3(get)]
    pub improvement_pct: F,
    /// Current radius scaling factor (starts at `discale`, decays to 1.0).
    #[pyo3(get)]
    pub radscale: F,
    /// Convergence precision target.
    #[pyo3(get)]
    pub precision: F,
    /// Per-relaxer acceptance rate this iteration, as
    /// ``[(target_index, rate), ...]``. Empty when no relaxers are
    /// attached.
    #[pyo3(get)]
    pub relaxer_acceptance: Vec<(usize, F)>,
}

impl PyStepInfo {
    fn from_info(info: &StepInfo) -> Self {
        Self {
            loop_idx: info.loop_idx,
            max_loops: info.max_loops,
            phase: info.phase.phase,
            total_phases: info.phase.total_phases,
            molecule_type: info.phase.molecule_type,
            fdist: info.fdist,
            frest: info.frest,
            improvement_pct: info.improvement_pct,
            radscale: info.radscale,
            precision: info.precision,
            relaxer_acceptance: info.relaxer_acceptance.clone(),
        }
    }
}

#[pymethods]
impl PyStepInfo {
    fn __repr__(&self) -> String {
        format!(
            "StepInfo(phase={}/{}, loop={}/{}, fdist={:.3e}, frest={:.3e}, improvement={:.2}%)",
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

// PyHandlerWrapper — bridges a Python object to the Rust `Handler` trait.
//
// Optional methods are resolved once at construction; missing ones stay
// `None` and are silently skipped. `stop_flag` lives behind an atomic
// because `Handler::should_stop` is `&self` while writes happen via the
// mutating `on_*` methods.

pub(crate) struct PyHandlerWrapper {
    on_start: Option<Py<PyAny>>,
    on_step: Option<Py<PyAny>>,
    on_finish: Option<Py<PyAny>>,
    stop_flag: Arc<AtomicBool>,
}

impl PyHandlerWrapper {
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

impl Handler for PyHandlerWrapper {
    fn on_start(&mut self, ntotat: usize, ntotmol: usize) {
        let Some(m) = &self.on_start else { return };
        Python::attach(|py| {
            if let Err(e) = m.bind(py).call1((ntotat, ntotmol)) {
                self.fail(e);
            }
        });
    }

    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        let Some(m) = &self.on_step else { return };
        Python::attach(|py| {
            let py_info = match Py::new(py, PyStepInfo::from_info(info)) {
                Ok(v) => v,
                Err(e) => {
                    self.fail(e);
                    return;
                }
            };
            match m.bind(py).call1((py_info,)) {
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
