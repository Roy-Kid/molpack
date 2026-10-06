//! Engine-entry bindings (engine-entry-split): `GenCanPack` and `CbmcGrow`
//! as 1:1 mirrors of the Rust entries — one terminal verb, one-shot
//! semantics, no sugar aliases.
//!
//! The shared `with_*` builders are bound ONCE, in `entry_pymethods!` —
//! the macro stamps the same forwarders into each entry's `#[pymethods]`
//! block, so the shared surface cannot drift between entries.
//!
//! `Pipeline` takes those same entry objects as its stages. How an entry
//! becomes one is stamped by the same macro, next to the entry
//! (`IntoStageFactory`), and *which* objects may be one is the single
//! `stage_entry_registry!` table below — both the dispatch and the
//! `TypeError` text come from it, so `PyPipeline` never names a concrete
//! entry type and a fourth entry registers itself in one line.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use molpack::{
    CbmcGrow, EarlyStopHandler, GenCanPack, LatticeGrow, LogLevel, PackEngine, StageFactory,
};
use pyo3::prelude::*;

use crate::errors::{pack_error_to_pyerr, take_err};
use crate::grow::{PyAnglePrior, PyTorsionPrior};
use crate::handler::PyHandlerWrapper;
use crate::parallel::rayon_compiled;
use crate::restraint::extract_restraint;
use crate::result::PyState;
use crate::target::PyTarget;

use molrs::op::types::F;

/// Shared entry knobs mirrored on the Python side; the Rust entry is built
/// at `run()` time.
#[derive(Default)]
struct SharedKnobs {
    tolerance: Option<F>,
    precision: Option<F>,
    seed: Option<u64>,
    periodic_box: Option<([F; 3], [F; 3])>,
    density: Option<F>,
    parallel_eval: Option<bool>,
    progress: bool,
    log_level: Option<LogLevel>,
    log_frequency: Option<usize>,
    py_handlers: Vec<Py<pyo3::types::PyAny>>,
    global_restraints: Vec<Py<pyo3::types::PyAny>>,
    consumed: bool,
}

impl SharedKnobs {
    fn clone_ref(&self, py: Python<'_>) -> Self {
        Self {
            tolerance: self.tolerance,
            precision: self.precision,
            seed: self.seed,
            periodic_box: self.periodic_box,
            density: self.density,
            parallel_eval: self.parallel_eval,
            progress: self.progress,
            log_level: self.log_level,
            log_frequency: self.log_frequency,
            py_handlers: self.py_handlers.iter().map(|h| h.clone_ref(py)).collect(),
            global_restraints: self
                .global_restraints
                .iter()
                .map(|r| r.clone_ref(py))
                .collect(),
            consumed: self.consumed,
        }
    }

    fn apply<E: PackEngine>(&self, py: Python<'_>, mut engine: E) -> PyResult<E> {
        if let Some(v) = self.tolerance {
            engine = engine.with_tolerance(v);
        }
        if let Some(v) = self.precision {
            engine = engine.with_precision(v);
        }
        if let Some(v) = self.seed {
            engine = engine.with_seed(v);
        }
        if let Some((min, max)) = self.periodic_box {
            engine = engine.with_periodic_box(min, max, [true; 3]);
        }
        if let Some(v) = self.density {
            engine = engine.with_density(v);
        }
        if let Some(on) = self.parallel_eval {
            engine = engine.with_parallel_eval(on);
        }
        for gr in &self.global_restraints {
            engine = engine.with_global_restraint(extract_restraint(gr.bind(py))?);
        }
        // An explicit level wins; otherwise the progress flag decides.
        let level = self.log_level.unwrap_or(if self.progress {
            LogLevel::Progress
        } else {
            LogLevel::Quiet
        });
        engine = engine.with_log_level(level);
        if let Some(n) = self.log_frequency {
            engine = engine.with_log_frequency(n);
        }
        // One flag per `apply` call — shared by the handlers this call
        // attaches, whether they land on a standalone entry or on one stage
        // of a pipeline. The lifecycle polls every handler
        // (`handlers.iter().any(|h| h.should_stop())`), so a stage handler
        // asking to stop still stops the whole run, and a raised exception
        // still reaches `finish_run` through the shared `take_err` stash.
        let stop_flag = Arc::new(AtomicBool::new(false));
        for py_h in &self.py_handlers {
            let wrapper = PyHandlerWrapper::new(py_h.clone_ref(py), Arc::clone(&stop_flag));
            engine = engine.with_handler(Box::new(wrapper));
        }
        Ok(engine)
    }

    fn guard_one_shot(&mut self) -> PyResult<()> {
        if self.consumed {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "this engine has already run: one engine, one run — build a new one",
            ));
        }
        self.consumed = true;
        Ok(())
    }
}

fn parse_log_level(level: &str) -> PyResult<LogLevel> {
    Ok(match level.to_ascii_lowercase().as_str() {
        "quiet" | "off" | "none" => LogLevel::Quiet,
        "summary" => LogLevel::Summary,
        "progress" | "thermo" => LogLevel::Progress,
        "verbose" | "debug" => LogLevel::Verbose,
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown log level {other:?}; expected quiet, summary, progress, or verbose"
            )));
        }
    })
}

fn finish_run(
    py: Python<'_>,
    result: Result<molpack::State, molpack::PackError>,
    periodic_box: Option<([F; 3], [F; 3])>,
) -> PyResult<PyState> {
    if let Some(py_err) = take_err() {
        return Err(py_err);
    }
    let mut result = result.map_err(pack_error_to_pyerr)?;
    crate::interop::stamp_box_bounds(&mut result.frame, periodic_box)?;
    let py_frame = crate::interop::frame_to_py(py, &result.frame)?.unbind();
    Ok(PyState {
        inner: result,
        py_frame,
    })
}

/// How an entry object becomes a pipeline stage.
///
/// Stamped by [`entry_pymethods!`] next to each entry, so the conversion
/// lives with the entry it converts and a new entry cannot forget to
/// provide one. The body is the entry's own `build_engine` — the very
/// construction its `run` uses, so a stage and a standalone run are the
/// same engine — followed by [`SharedKnobs::apply`], which attaches this
/// object's Python handlers (the pipeline *adopts* them) and carries any
/// non-default shared knob through so the Rust side refuses it **by name**.
/// The binding neither strips a knob nor drops a handler.
///
/// `to_` and not `into_`: the entry object stays owned by the interpreter,
/// so the conversion borrows it and the caller's stage object is still
/// usable afterwards.
trait IntoStageFactory {
    fn to_stage_factory(&self, py: Python<'_>) -> PyResult<Box<dyn StageFactory>>;
}

/// One `#[pymethods]` block per entry: the shared builders (stamped
/// identically for every entry), the entry's own methods, and its
/// [`IntoStageFactory`] body — which is why every entry must define
/// `fn build_engine(&self, py: Python<'_>) -> PyResult<impl PackEngine>`.
macro_rules! entry_pymethods {
    ($ty:ty { $($extra:tt)* }) => {
        #[pymethods]
        impl $ty {
            fn with_tolerance(&self, v: F) -> Self {
                let mut c = self.clone_fields();
                c.shared.tolerance = Some(v);
                c
            }
            fn with_precision(&self, v: F) -> Self {
                let mut c = self.clone_fields();
                c.shared.precision = Some(v);
                c
            }
            fn with_seed(&self, v: u64) -> Self {
                let mut c = self.clone_fields();
                c.shared.seed = Some(v);
                c
            }
            fn with_periodic_box(&self, min: [F; 3], max: [F; 3]) -> Self {
                let mut c = self.clone_fields();
                c.shared.periodic_box = Some((min, max));
                c
            }
            fn with_density(&self, rho: F) -> Self {
                let mut c = self.clone_fields();
                c.shared.density = Some(rho);
                c
            }
            /// Run pair-kernel reductions on rayon worker threads. Errors if
            /// this wheel was built without the ``rayon`` feature.
            fn with_parallel_eval(&self, enabled: bool) -> PyResult<Self> {
                if enabled && !rayon_compiled() {
                    return Err(pyo3::exceptions::PyRuntimeError::new_err(
                        "parallel evaluation requested but this molpack wheel was built \
                         without the `rayon` feature; rebuild with `maturin develop --release` \
                         (rayon is enabled by default) or check `molpack.rayon_enabled()`",
                    ));
                }
                let mut c = self.clone_fields();
                c.shared.parallel_eval = Some(enabled);
                Ok(c)
            }
            /// Append a global restraint — broadcast to every target at run
            /// time.
            fn with_global_restraint(&self, restraint: Py<pyo3::types::PyAny>) -> Self {
                let mut c = self.clone_fields();
                c.shared.global_restraints.push(restraint);
                c
            }
            #[pyo3(signature = (on = true))]
            fn with_progress(&self, on: bool) -> Self {
                let mut c = self.clone_fields();
                c.shared.progress = on;
                c
            }
            /// Screen-log detail: ``"quiet"`` / ``"summary"`` / ``"progress"``
            /// / ``"verbose"``. An explicit level wins over ``with_progress``.
            fn with_log_level(&self, level: &str) -> PyResult<Self> {
                let parsed = parse_log_level(level)?;
                let mut c = self.clone_fields();
                c.shared.log_level = Some(parsed);
                Ok(c)
            }
            fn with_log_frequency(&self, n: usize) -> Self {
                let mut c = self.clone_fields();
                c.shared.log_frequency = Some(n.max(1));
                c
            }
            /// Append a Python handler. See :class:`StepInfo` for the
            /// callback contract.
            fn with_handler(&self, handler: Py<pyo3::types::PyAny>) -> Self {
                let mut c = self.clone_fields();
                c.shared.py_handlers.push(handler);
                c
            }

            $($extra)*
        }

        impl IntoStageFactory for $ty {
            fn to_stage_factory(&self, py: Python<'_>) -> PyResult<Box<dyn StageFactory>> {
                Ok(Box::new(self.shared.apply(py, self.build_engine(py)?)?))
            }
        }
    };
}

/// The one table of entry pyclasses a `Pipeline` accepts as a stage.
///
/// It generates both halves of the answer — the ordered dispatch and the
/// `TypeError` text listing what was expected — so the two cannot drift and
/// a new entry costs exactly one line, next to its own pyclass.
macro_rules! stage_entry_registry {
    ($($ty:ty => $name:literal),+ $(,)?) => {
        /// The accepted entries, in dispatch order.
        const STAGE_ENTRIES: &[&str] = &[$($name),+];

        fn not_a_stage(obj: &Bound<'_, pyo3::types::PyAny>) -> PyErr {
            let got = obj
                .get_type()
                .name()
                .and_then(|n| n.extract::<String>())
                .unwrap_or_else(|_| "?".to_owned());
            pyo3::exceptions::PyTypeError::new_err(format!(
                "stage must be one of {} (got {got})",
                STAGE_ENTRIES.join(", ")
            ))
        }

        /// Type-check a candidate stage where the caller hands it over, so a
        /// pipeline holding something it cannot run is never constructible.
        fn check_stage(obj: &Bound<'_, pyo3::types::PyAny>) -> PyResult<()> {
            $(if obj.is_instance_of::<$ty>() {
                return Ok(());
            })+
            Err(not_a_stage(obj))
        }

        /// Convert a stage object into the factory the Rust pipeline takes.
        fn stage_from_py(
            py: Python<'_>,
            obj: &Bound<'_, pyo3::types::PyAny>,
        ) -> PyResult<Box<dyn StageFactory>> {
            $(if let Ok(entry) = obj.cast::<$ty>() {
                return entry.borrow().to_stage_factory(py);
            })+
            Err(not_a_stage(obj))
        }
    };
}

/// Rigid-body GENCAN packing entry (1:1 mirror of the Rust `GenCanPack`).
#[pyclass(name = "GenCanPack")]
pub struct PyGenCanPack {
    shared: SharedKnobs,
    seed: Option<molpack::State>,
    inner_iterations: Option<usize>,
    init_passes: Option<usize>,
    init_box_half_size: Option<F>,
    perturb: Option<(F, bool, bool)>,
    avoid_overlap: Option<bool>,
    /// `None`: the entry's default early stop; `Some(None)`: switched off.
    early_stop: Option<Option<EarlyStopHandler>>,
}

entry_pymethods!(PyGenCanPack {
    #[new]
    fn new() -> Self {
        Self {
            shared: SharedKnobs::default(),
            seed: None,
            inner_iterations: None,
            init_passes: None,
            init_box_half_size: None,
            perturb: None,
            avoid_overlap: None,
            early_stop: None,
        }
    }

    /// Continue on a previous run's placement solution (the explicit
    /// push-off chain): the free copies start EXACTLY where ``result``
    /// left them; the cell travels with the seed — do not declare a box,
    /// density, or cell on a seeded engine.
    fn with_restart(&self, result: &PyState) -> Self {
        let mut c = self.clone_fields();
        c.seed = Some(result.inner.clone());
        c
    }
    fn with_inner_iterations(&self, n: usize) -> Self {
        let mut c = self.clone_fields();
        c.inner_iterations = Some(n);
        c
    }
    fn with_init_passes(&self, n: usize) -> Self {
        let mut c = self.clone_fields();
        c.init_passes = Some(n);
        c
    }
    fn with_init_box_half_size(&self, half_size: F) -> Self {
        let mut c = self.clone_fields();
        c.init_box_half_size = Some(half_size);
        c
    }
    #[pyo3(signature = (fraction, random = false, enabled = true))]
    fn with_perturb(&self, fraction: F, random: bool, enabled: bool) -> Self {
        let mut c = self.clone_fields();
        c.perturb = Some((fraction, random, enabled));
        c
    }
    #[pyo3(signature = (on = true))]
    fn with_avoid_overlap(&self, on: bool) -> Self {
        let mut c = self.clone_fields();
        c.avoid_overlap = Some(on);
        c
    }
    /// End a GENCAN phase once, at the user's radii (``radscale == 1``),
    /// Packmol's ``bestf`` has improved by less than ``threshold_pct``
    /// percent over ``patience`` loops. On by default with these values
    /// (10 % is Packmol's movebad threshold); ``enabled=False`` runs every
    /// phase to ``max_loops`` as Packmol does. A stopped final phase leaves
    /// ``State.converged`` false.
    #[pyo3(signature = (enabled = true, *, threshold_pct = 10.0, patience = 10))]
    fn with_early_stop(&self, enabled: bool, threshold_pct: F, patience: usize) -> Self {
        let mut c = self.clone_fields();
        c.early_stop =
            Some(enabled.then(|| EarlyStopHandler::new(threshold_pct).with_patience(patience)));
        c
    }

    /// Run the packing. One engine, one run. ``max_loops`` defaults to
    /// Packmol's ``nloop``: ``200 * len(targets)``.
    #[pyo3(signature = (targets, max_loops = None))]
    fn run(
        &mut self,
        py: Python<'_>,
        targets: Vec<PyTarget>,
        max_loops: Option<usize>,
    ) -> PyResult<PyState> {
        self.shared.guard_one_shot()?;
        let max_loops = max_loops.unwrap_or(GenCanPack::default_max_loops(targets.len()));
        let rust_targets: Vec<_> = targets.into_iter().map(|t| t.inner).collect();
        let engine = self.shared.apply(py, self.build_engine(py)?)?;
        let pb = self.shared.periodic_box;
        finish_run(py, engine.run(&rust_targets, max_loops), pb)
    }

    fn __repr__(&self) -> String {
        "GenCanPack(...)".to_string()
    }
});

impl PyGenCanPack {
    /// Script-loader constructor: shared knobs pre-set from a parsed `.inp`.
    pub(crate) fn from_script(
        tolerance: Option<F>,
        seed: Option<u64>,
        periodic_box: Option<([F; 3], [F; 3])>,
    ) -> Self {
        Self {
            shared: SharedKnobs {
                tolerance,
                seed,
                periodic_box,
                ..SharedKnobs::default()
            },
            seed: None,
            inner_iterations: None,
            init_passes: None,
            init_box_half_size: None,
            perturb: None,
            avoid_overlap: None,
            early_stop: None,
        }
    }

    /// The Rust engine these knobs describe — one home for that
    /// construction, read by `run` and by the pipeline conversion alike.
    fn build_engine(&self, _py: Python<'_>) -> PyResult<GenCanPack> {
        let mut engine = GenCanPack::new();
        if let Some(seed) = &self.seed {
            engine = engine.with_restart(seed);
        }
        if let Some(n) = self.inner_iterations {
            engine = engine.with_inner_iterations(n);
        }
        if let Some(n) = self.init_passes {
            engine = engine.with_init_passes(n);
        }
        if let Some(h) = self.init_box_half_size {
            engine = engine.with_init_box_half_size(h);
        }
        if let Some((f, r, e)) = self.perturb {
            engine = engine.with_perturb(f, r, e);
        }
        if let Some(on) = self.avoid_overlap {
            engine = engine.with_avoid_overlap(on);
        }
        if let Some(early_stop) = &self.early_stop {
            engine = engine.with_early_stop(early_stop.clone());
        }
        Ok(engine)
    }

    fn clone_fields(&self) -> Self {
        Python::attach(|py| Self {
            shared: self.shared.clone_ref(py),
            seed: self.seed.clone(),
            inner_iterations: self.inner_iterations,
            init_passes: self.init_passes,
            init_box_half_size: self.init_box_half_size,
            perturb: self.perturb,
            avoid_overlap: self.avoid_overlap,
            early_stop: self.early_stop.clone(),
        })
    }
}

/// Configurational-bias chain-growth entry (1:1 mirror of `CbmcGrow`).
/// The torsion prior is the one mandatory constructor argument.
#[pyclass(name = "CbmcGrow")]
pub struct PyCbmcGrow {
    shared: SharedKnobs,
    prior: molpack::grow::TorsionPrior,
    trials: Option<usize>,
    retract: Option<usize>,
    relax: Option<(usize, usize)>,
    selectivity: Option<F>,
    soften_after: Option<usize>,
    min_hard_scale: Option<F>,
    angle_prior: Option<molpack::grow::AnglePrior>,
    soft_shell: Option<F>,
    serial: bool,
    void_bias: bool,
}

entry_pymethods!(PyCbmcGrow {
    #[new]
    fn new(torsion_prior: &PyTorsionPrior) -> Self {
        Self {
            shared: SharedKnobs::default(),
            prior: torsion_prior.inner.clone(),
            trials: None,
            retract: None,
            relax: None,
            selectivity: None,
            soften_after: None,
            min_hard_scale: None,
            angle_prior: None,
            soft_shell: None,
            serial: false,
            void_bias: false,
        }
    }

    fn with_trials(&self, n: usize) -> Self {
        let mut c = self.clone_fields();
        c.trials = Some(n);
        c
    }
    fn with_retract(&self, n: usize) -> Self {
        let mut c = self.clone_fields();
        c.retract = Some(n);
        c
    }
    fn with_relax(&self, every: usize, window: usize) -> Self {
        let mut c = self.clone_fields();
        c.relax = Some((every, window));
        c
    }
    fn with_selectivity(&self, beta: F) -> Self {
        let mut c = self.clone_fields();
        c.selectivity = Some(beta);
        c
    }
    fn with_soften_after(&self, n: usize) -> Self {
        let mut c = self.clone_fields();
        c.soften_after = Some(n);
        c
    }
    fn with_min_hard_scale(&self, s: F) -> Self {
        let mut c = self.clone_fields();
        c.min_hard_scale = Some(s);
        c
    }
    fn with_angle_prior(&self, prior: &PyAnglePrior) -> Self {
        let mut c = self.clone_fields();
        c.angle_prior = Some(prior.inner.clone());
        c
    }
    fn with_soft_shell(&self, width: F) -> Self {
        let mut c = self.clone_fields();
        c.soft_shell = Some(width);
        c
    }
    #[pyo3(signature = (serial = true))]
    fn with_serial(&self, serial: bool) -> Self {
        let mut c = self.clone_fields();
        c.serial = serial;
        c
    }
    #[pyo3(signature = (void_bias = true))]
    fn with_void_bias(&self, void_bias: bool) -> Self {
        let mut c = self.clone_fields();
        c.void_bias = void_bias;
        c
    }

    /// Run the growth. One engine, one run. Reports honestly — nothing else
    /// runs on non-convergence. For the rigid push-off, chain explicitly:
    /// ``GenCanPack().with_restart(result).run(same_targets, ...)``.
    fn run(
        &mut self,
        py: Python<'_>,
        targets: Vec<PyTarget>,
        max_loops: usize,
    ) -> PyResult<PyState> {
        self.shared.guard_one_shot()?;
        let rust_targets: Vec<_> = targets.into_iter().map(|t| t.inner).collect();
        let engine = self.shared.apply(py, self.build_engine(py)?)?;
        let pb = self.shared.periodic_box;
        finish_run(py, engine.run(&rust_targets, max_loops), pb)
    }

    fn __repr__(&self) -> String {
        "CbmcGrow(...)".to_string()
    }
});

/// Diamond-lattice growth entry (1:1 mirror of `LatticeGrow`): melt-density
/// chain generation as an on-lattice SAW, decorated back to all-atom
/// geometry. Residual contacts are honest; chain a seeded ``GenCanPack``.
#[pyclass(name = "LatticeGrow")]
pub struct PyLatticeGrow {
    shared: SharedKnobs,
    prior: molpack::grow::TorsionPrior,
    occupancy_guard: Option<bool>,
}

entry_pymethods!(PyLatticeGrow {
    #[new]
    fn new(torsion_prior: &PyTorsionPrior) -> Self {
        Self {
            shared: SharedKnobs::default(),
            prior: torsion_prior.inner.clone(),
            occupancy_guard: None,
        }
    }

    /// Nearest-neighbour site exclusion (default on).
    #[pyo3(signature = (on = true))]
    fn with_occupancy_guard(&self, on: bool) -> Self {
        let mut c = self.clone_fields();
        c.occupancy_guard = Some(on);
        c
    }

    /// Run the lattice growth. One engine, one run.
    fn run(
        &mut self,
        py: Python<'_>,
        targets: Vec<PyTarget>,
        max_loops: usize,
    ) -> PyResult<PyState> {
        self.shared.guard_one_shot()?;
        let rust_targets: Vec<_> = targets.into_iter().map(|t| t.inner).collect();
        let engine = self.shared.apply(py, self.build_engine(py)?)?;
        let pb = self.shared.periodic_box;
        finish_run(py, engine.run(&rust_targets, max_loops), pb)
    }

    fn __repr__(&self) -> String {
        "LatticeGrow(...)".to_string()
    }
});

impl PyLatticeGrow {
    /// The Rust engine these knobs describe — one home for that
    /// construction, read by `run` and by the pipeline conversion alike.
    fn build_engine(&self, _py: Python<'_>) -> PyResult<LatticeGrow> {
        let mut engine = LatticeGrow::new(self.prior.clone());
        if let Some(on) = self.occupancy_guard {
            engine = engine.with_occupancy_guard(on);
        }
        Ok(engine)
    }

    fn clone_fields(&self) -> Self {
        Python::attach(|py| Self {
            shared: self.shared.clone_ref(py),
            prior: self.prior.clone(),
            occupancy_guard: self.occupancy_guard,
        })
    }
}

impl PyCbmcGrow {
    /// The Rust engine these knobs describe — one home for that
    /// construction, read by `run` and by the pipeline conversion alike.
    fn build_engine(&self, _py: Python<'_>) -> PyResult<CbmcGrow> {
        let mut engine = CbmcGrow::new(self.prior.clone());
        if let Some(n) = self.trials {
            engine = engine.with_trials(n);
        }
        if let Some(n) = self.retract {
            engine = engine.with_retract(n);
        }
        if let Some((e, w)) = self.relax {
            engine = engine.with_relax(e, w);
        }
        if let Some(b) = self.selectivity {
            engine = engine.with_selectivity(b);
        }
        if let Some(n) = self.soften_after {
            engine = engine.with_soften_after(n);
        }
        if let Some(s) = self.min_hard_scale {
            engine = engine.with_min_hard_scale(s);
        }
        if let Some(ref ap) = self.angle_prior {
            engine = engine.with_angle_prior(ap.clone());
        }
        if let Some(w) = self.soft_shell {
            engine = engine.with_soft_shell(w);
        }
        Ok(engine
            .with_serial(self.serial)
            .with_void_bias(self.void_bias))
    }

    fn clone_fields(&self) -> Self {
        Python::attach(|py| Self {
            shared: self.shared.clone_ref(py),
            prior: self.prior.clone(),
            trials: self.trials,
            retract: self.retract,
            relax: self.relax,
            selectivity: self.selectivity,
            soften_after: self.soften_after,
            min_hard_scale: self.min_hard_scale,
            angle_prior: self.angle_prior.clone(),
            soft_shell: self.soft_shell,
            serial: self.serial,
            void_bias: self.void_bias,
        })
    }
}

// The single table: which entry objects may be a pipeline stage. A fourth
// entry adds one line here, next to its own pyclass — the dispatch and the
// `TypeError` text follow, and `PyPipeline` below stays untouched.
stage_entry_registry!(
    PyGenCanPack => "GenCanPack",
    PyCbmcGrow => "CbmcGrow",
    PyLatticeGrow => "LatticeGrow",
);

/// A sequence of stages behind one lifecycle, one settings set and one
/// handler set (1:1 mirror of the Rust `Pipeline`).
///
/// The shared knobs are the run's — declare them here, not on a stage. A
/// stage's own Python handlers travel with it and are **adopted**: they
/// observe the whole run. A stage carrying a non-default shared knob is
/// refused by name at `run` (`ValueError`), never silently stripped.
///
/// The stage objects are kept as the Python objects the caller handed over
/// and turned into Rust factories at `run` — `Box<dyn StageFactory>` is not
/// clonable, while these consuming builders hand back a fresh pipeline on
/// every call. The type check does not wait for `run`: it fires where the
/// stage is handed over, so an unrunnable pipeline is not constructible.
#[pyclass(name = "Pipeline")]
pub struct PyPipeline {
    shared: SharedKnobs,
    stages: Vec<Py<pyo3::types::PyAny>>,
}

// The shared macro also stamps `IntoStageFactory` on `PyPipeline` (a Rust
// `Pipeline` is itself a `StageFactory`), but `PyPipeline` is deliberately
// NOT in `stage_entry_registry!`: nesting pipelines stays a Rust-only
// extension point in v1, so that impl has no caller by design.
entry_pymethods!(PyPipeline {
    #[new]
    #[pyo3(signature = (stages = None))]
    fn new(py: Python<'_>, stages: Option<Vec<Py<pyo3::types::PyAny>>>) -> PyResult<Self> {
        let stages = stages.unwrap_or_default();
        for stage in &stages {
            check_stage(stage.bind(py))?;
        }
        Ok(Self {
            shared: SharedKnobs::default(),
            stages,
        })
    }

    /// Append a stage — one of the entry objects the registry accepts.
    fn with_stage(&self, py: Python<'_>, stage: Py<pyo3::types::PyAny>) -> PyResult<Self> {
        check_stage(stage.bind(py))?;
        let mut c = self.clone_fields();
        c.stages.push(stage);
        Ok(c)
    }

    /// Run the chain. One pipeline, one run.
    fn run(
        &mut self,
        py: Python<'_>,
        targets: Vec<PyTarget>,
        max_loops: usize,
    ) -> PyResult<PyState> {
        self.shared.guard_one_shot()?;
        let rust_targets: Vec<_> = targets.into_iter().map(|t| t.inner).collect();
        let engine = self.shared.apply(py, self.build_engine(py)?)?;
        let pb = self.shared.periodic_box;
        finish_run(py, engine.run(&rust_targets, max_loops), pb)
    }

    fn __repr__(&self) -> String {
        format!("Pipeline({} stages)", self.stages.len())
    }
});

impl PyPipeline {
    /// The Rust engine these knobs describe: one factory per stage object,
    /// in the order they were added, each converted through the registry.
    /// Nothing here names a concrete entry type.
    fn build_engine(&self, py: Python<'_>) -> PyResult<molpack::Pipeline> {
        let mut pipeline = molpack::Pipeline::new();
        for stage in &self.stages {
            pipeline = pipeline.with_stage(stage_from_py(py, stage.bind(py))?);
        }
        Ok(pipeline)
    }

    fn clone_fields(&self) -> Self {
        Python::attach(|py| Self {
            shared: self.shared.clone_ref(py),
            stages: self.stages.iter().map(|s| s.clone_ref(py)).collect(),
        })
    }
}
