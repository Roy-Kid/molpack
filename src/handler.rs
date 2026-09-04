//! Handler trait and built-in handlers for packing progress callbacks.

use molrs::types::F;
use std::io::BufWriter;
use std::path::PathBuf;
use std::time::Instant;

use crate::context::PackContext;
use crate::frame::compute_mol_ids;
use crate::numerics::objective_small_floor;
use crate::stage::StageOutcome;

// ── Info structs ─────────────────────────────────────────────────────────────

/// Identifies the stage a callback comes from.
///
/// The field shape follows [`PhaseInfo`] on purpose — a stage is to a run
/// what a phase is to the GENCAN loop, so the two identities read the same
/// way — and, like `PhaseInfo`, this is a plain `Copy` record a caller may
/// build by literal (a handler test drives the two stage hooks with one).
#[derive(Debug, Clone, Copy)]
pub struct StageInfo {
    /// 0-based index of this stage in the run.
    pub index: usize,
    /// How many stages the run has. A single-stage run reports `1`.
    pub total: usize,
    /// The stage's own [`Stage::name`](crate::Stage::name).
    pub name: &'static str,
}

/// Information about the current packing phase.
#[derive(Debug, Clone, Copy)]
pub struct PhaseInfo {
    /// 0-based phase index.
    pub phase: usize,
    /// Total number of phases (ntype + 1).
    pub total_phases: usize,
    /// If `Some(itype)`, this is a per-type compaction phase.
    /// If `None`, this is the final all-types phase.
    pub molecule_type: Option<usize>,
}

/// Summary report emitted at the end of each packing phase.
///
/// Passed to [`Handler::on_phase_end`]. Fields mirror the counters
/// reported by `ProgressHandler::on_phase_start`; when no summary is
/// needed (e.g. early termination), callers may omit the handler hook.
#[derive(Debug, Clone, Copy, Default)]
pub struct PhaseReport {
    /// Number of outer iterations actually run in this phase.
    pub iterations: usize,
    /// Final `fdist` (max inter-molecular overlap) at phase exit.
    pub fdist: F,
    /// Final `frest` (max restraint violation) at phase exit.
    pub frest: F,
    /// Whether the phase converged under `precision`.
    pub converged: bool,
}

/// Screen-log detail level for LAMMPS-style packer output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum LogLevel {
    /// Print nothing.
    #[default]
    Quiet,
    /// Print system setup, phase summaries, and final summary.
    Summary,
    /// Print per-step thermo-style progress lines.
    Progress,
    /// Print the same progress lines plus extra diagnostic columns.
    Verbose,
}

impl LogLevel {
    #[inline]
    pub const fn is_enabled(self) -> bool {
        !matches!(self, Self::Quiet)
    }
}

/// Per-iteration progress snapshot.
///
/// Both solvers fill this in, and three fields carry a different meaning on
/// each path — see `loop_idx`, `max_loops` and `radscale` below. Growth also
/// reports `fdist` and `frest` as a constructive `0.0` on every round, because
/// it only ever commits a placement that already clears the hard core and the
/// restraints; the measured end-of-run numbers live in
/// [`PackResult`](crate::PackResult).
///
/// `#[non_exhaustive]`: the crate builds this in exactly three places (the
/// GENCAN iteration, the two growth drivers), and every stage that lands
/// later adds a field. An external construction site would turn each of
/// those additions into a breaking change for a struct nobody outside this
/// crate emits — readers are unaffected.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct StepInfo {
    /// The stage that emitted this step.
    pub stage: StageInfo,
    /// GENCAN: 0-based loop iteration within the current phase. Growth: the
    /// 1-based index of the current round.
    pub loop_idx: usize,
    /// GENCAN: maximum loops for this phase. Growth: the caller's `max_loops`
    /// verbatim — *not* the driver's round cap, which is `max_loops × (the
    /// longest chain's steps + 1)`; see [`grow::driver`](crate::grow::driver).
    pub max_loops: usize,
    /// Current phase info.
    pub phase: PhaseInfo,
    /// Max inter-molecular overlap violation (0.0 = no overlap).
    pub fdist: F,
    /// Max constraint violation (0.0 = all constraints satisfied).
    pub frest: F,
    /// Improvement from last iteration, as percentage (positive = improving).
    pub improvement_pct: F,
    /// GENCAN: current radius scaling factor (starts at discale, decays to
    /// 1.0). Growth: the dimensionless hard-core scale — the factor multiplying
    /// the pair contact distance a placement must clear (`1.0` = full declared
    /// contact) — which starts at 1.0 and steps down in rungs of 0.97 to the
    /// floor set by
    /// [`GrowConfig::with_min_hard_scale`](crate::grow::GrowConfig::with_min_hard_scale).
    pub radscale: F,
    /// Convergence precision target.
    pub precision: F,
    /// Relaxer acceptance rates: `(type_index, acceptance_rate)`.
    pub relaxer_acceptance: Vec<(usize, F)>,
}

// ── Trait ─────────────────────────────────────────────────────────────────────

/// Callback interface called by the [`PackEngine`](crate::PackEngine) lifecycle during packing.
pub trait Handler: Send {
    /// Called immediately at the start of [`run`][crate::PackEngine::run],
    /// before any computation. Use this for immediate user feedback.
    fn on_start(&mut self, _ntotat: usize, _ntotmol: usize) {}

    /// Called once after initialization completes, with valid `xcart` positions.
    /// Use this to write the initial conformation (e.g. [`XYZHandler`]).
    fn on_initialized(&mut self, _sys: &PackContext) {}

    /// Called after each outer optimization loop iteration.
    fn on_step(&mut self, info: &StepInfo, sys: &PackContext);

    /// Called at the start of each packing phase (per-type and all-types).
    /// Allows stateful handlers to reset between phases.
    fn on_phase_start(&mut self, _info: &PhaseInfo) {}

    /// Called once after the packing loop finishes (convergence or max loops).
    fn on_finish(&mut self, _sys: &PackContext) {}

    /// Return `true` to request early termination of the packing loop.
    fn should_stop(&self) -> bool {
        false
    }

    // ── v2 additions (spec §6.6) — default no-op, backward compatible ──

    /// Called after each inner GENCAN iteration (more granular than [`on_step`]).
    ///
    /// Default: no-op. Implement when you need per-inner-iteration feedback
    /// (e.g. plotting objective evolution, adaptive stop criteria).
    ///
    /// - `iter` — 0-based inner iteration counter within the current outer step
    /// - `f` — current objective value (pair + restraint)
    /// - `sys` — read-only view of the packing context
    ///
    /// [`on_step`]: Handler::on_step
    fn on_inner_iter(&mut self, _iter: u32, _f: F, _sys: &PackContext) {}

    /// Called at the end of each packing phase, with a summary report.
    ///
    /// Default: no-op. Paired with [`on_phase_start`] for symmetric
    /// setup / teardown hooks.
    ///
    /// [`on_phase_start`]: Handler::on_phase_start
    fn on_phase_end(&mut self, _info: &PhaseInfo, _report: &PhaseReport) {}

    // ── stage-pipeline-04-stage additions — default no-op, backward compatible ──

    /// Called before a stage starts, with the stage's identity.
    ///
    /// Default: no-op. The **call** belongs to the pipeline that chains
    /// stages; a single-stage run reports `index == 0` and `total == 1`.
    fn on_stage_start(&mut self, _info: &StageInfo) {}

    /// Called after a stage returns, with its identity, its outcome, and the
    /// state it just finished writing.
    ///
    /// Default: no-op. The **call** belongs to the pipeline that chains
    /// stages; a single-stage run reports `index == 0` and `total == 1`.
    ///
    /// [`StageOutcome`] deliberately carries no verdict. A handler that wants
    /// the run's violation maxima reads them off `sys` — `sys.fdist` and
    /// `sys.frest`, the shared objective's numbers on the post-run state,
    /// exactly as [`on_finish`] does. Same shape, same authority.
    ///
    /// [`on_finish`]: Handler::on_finish
    fn on_stage_end(&mut self, _info: &StageInfo, _outcome: &StageOutcome, _sys: &PackContext) {}
}

// ── NullHandler ───────────────────────────────────────────────────────────────

/// A no-op handler.
pub struct NullHandler;

impl Handler for NullHandler {
    fn on_step(&mut self, _info: &StepInfo, _sys: &PackContext) {}
}

// ── XYZHandler ────────────────────────────────────────────────────────────────

/// Writes packing snapshots as a multi-frame XYZ trajectory.
///
/// Writes a frame on every `every`-th step (starting from step 0).
/// No automatic initial or final writes.
pub struct XYZHandler {
    path: PathBuf,
    /// Write every `n` steps (must be >= 1).
    every: usize,
    file: Option<BufWriter<std::fs::File>>,
    /// Precomputed global molecule ID per atom (constant across all frames).
    mol_ids: Vec<usize>,
}

impl XYZHandler {
    /// Create a new XYZ trajectory handler.
    ///
    /// `every` controls writing frequency: a frame is written on every step
    /// where `loop_idx % every == 0` (step 0 is always included).
    /// `every` must be >= 1.
    pub fn new(path: impl Into<PathBuf>, every: usize) -> Self {
        assert!(every >= 1, "every must be >= 1");
        Self {
            path: path.into(),
            every,
            file: None,
            mol_ids: Vec::new(),
        }
    }

    fn open(&mut self) {
        if self.file.is_some() {
            return;
        }
        match std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&self.path)
        {
            Ok(f) => self.file = Some(BufWriter::new(f)),
            Err(e) => log::warn!("XYZHandler: cannot open {}: {e}", self.path.display()),
        }
    }

    fn write_frame(&mut self, comment: &str, sys: &PackContext) {
        self.open();
        let Some(ref mut w) = self.file else { return };
        use std::io::Write;
        let nat = sys.xcart.len();
        let _ = writeln!(w, "{nat}");
        let _ = writeln!(w, "Properties=species:S:1:pos:R:3:mol:I:1  {comment}");
        for (icart, pos) in sys.xcart.iter().enumerate() {
            let elem = sys
                .elements
                .get(icart)
                .and_then(|e| *e)
                .map(|e| e.symbol())
                .unwrap_or("X");
            let mol_id = self.mol_ids.get(icart).copied().unwrap_or(0);
            let _ = writeln!(
                w,
                "{elem}  {:.6}  {:.6}  {:.6}  {mol_id}",
                pos[0], pos[1], pos[2]
            );
        }
        let _ = w.flush();
    }
}

impl Handler for XYZHandler {
    fn on_initialized(&mut self, sys: &PackContext) {
        self.mol_ids = compute_mol_ids(sys);
    }

    fn on_step(&mut self, info: &StepInfo, sys: &PackContext) {
        if info.loop_idx.is_multiple_of(self.every) {
            self.write_frame(&format!("step {}", info.loop_idx), sys);
        }
    }
}

// ── ProgressHandler ───────────────────────────────────────────────────────────

/// Prints human-readable progress lines to `stderr`.
///
/// Installed by the engine lifecycle according to [`LogSpec`](crate::entry::LogSpec).
pub struct ProgressHandler {
    start: Option<Instant>,
}

impl ProgressHandler {
    pub fn new() -> Self {
        Self { start: None }
    }
}

impl Default for ProgressHandler {
    fn default() -> Self {
        Self::new()
    }
}

impl Handler for ProgressHandler {
    fn on_start(&mut self, ntotat: usize, ntotmol: usize) {
        self.start = Some(Instant::now());
        eprintln!("Packing {ntotmol} molecules ({ntotat} atoms)...");
    }

    fn on_initialized(&mut self, sys: &PackContext) {
        let elapsed = self.start.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
        eprintln!(
            "  Initializing... done ({:.1}s)  overlap: {:.4e}  constraints: {:.4e}",
            elapsed, sys.fdist, sys.frest
        );
    }

    fn on_phase_start(&mut self, info: &PhaseInfo) {
        let desc = match info.molecule_type {
            Some(itype) => format!("Compacting type {itype}"),
            None => "Optimizing all types together".to_string(),
        };
        eprintln!("  Phase [{}/{}] {desc}", info.phase + 1, info.total_phases,);
    }

    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        let elapsed = self.start.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
        eprintln!(
            "    Step [{}/{}]  overlap: {:.2e}  constraints: {:.2e}  improved {:.1}%  ({:.1}s)",
            info.loop_idx + 1,
            info.max_loops,
            info.fdist,
            info.frest,
            info.improvement_pct,
            elapsed,
        );
    }

    fn on_finish(&mut self, sys: &PackContext) {
        let elapsed = self.start.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
        if sys.fdist < 0.01 && sys.frest < 0.01 {
            eprintln!(
                "  Converged in {:.1}s — overlap: {:.2e}  constraints: {:.2e}",
                elapsed, sys.fdist, sys.frest,
            );
        } else {
            eprintln!(
                "  Did not converge ({:.1}s) — overlap: {:.2e}  constraints: {:.2e}",
                elapsed, sys.fdist, sys.frest,
            );
        }
    }
}

// ── LammpsLogHandler ─────────────────────────────────────────────────────────

/// LAMMPS-style screen log for packing runs.
///
/// Most users should enable this through
/// [`PackEngine::with_log_level`][crate::PackEngine::with_log_level]
/// or [`with_log_frequency`][crate::PackEngine::with_log_frequency]
/// instead of attaching the handler manually.
pub struct LammpsLogHandler {
    level: LogLevel,
    every: usize,
    tolerance: F,
    precision: F,
    seed: u64,
    max_loops: usize,
    ntypes: usize,
    periodic_box: Option<([F; 3], [F; 3], [bool; 3])>,
    start: Option<Instant>,
    phase_start: Option<Instant>,
}

impl LammpsLogHandler {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        level: LogLevel,
        every: usize,
        tolerance: F,
        precision: F,
        seed: u64,
        max_loops: usize,
        ntypes: usize,
        periodic_box: Option<([F; 3], [F; 3], [bool; 3])>,
    ) -> Self {
        Self {
            level,
            every: every.max(1),
            tolerance,
            precision,
            seed,
            max_loops,
            ntypes,
            periodic_box,
            start: None,
            phase_start: None,
        }
    }

    fn elapsed(&self) -> F {
        self.start.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0)
    }

    fn phase_elapsed(&self) -> F {
        self.phase_start
            .map(|t| t.elapsed().as_secs_f64())
            .unwrap_or(0.0)
    }
}

impl Handler for LammpsLogHandler {
    fn on_start(&mut self, ntotat: usize, ntotmol: usize) {
        if !self.level.is_enabled() {
            return;
        }

        let now = Instant::now();
        self.start = Some(now);
        self.phase_start = Some(now);

        eprintln!("Molpack screen log");
        eprintln!("System information:");
        eprintln!("  molecule types = {}", self.ntypes);
        eprintln!("  molecules      = {ntotmol}");
        eprintln!("  atoms          = {ntotat}");
        eprintln!("Settings:");
        eprintln!("  tolerance      = {:.6} A", self.tolerance);
        eprintln!("  precision      = {:.6}", self.precision);
        eprintln!("  seed           = {}", self.seed);
        eprintln!("  nloop          = {}", self.max_loops);
        match self.periodic_box {
            Some((min, max, flags)) => eprintln!(
                "  pbc            = [{:.6}, {:.6}, {:.6}] -> [{:.6}, {:.6}, {:.6}]  flags={:?}",
                min[0], min[1], min[2], max[0], max[1], max[2], flags
            ),
            None => eprintln!("  pbc            = off"),
        }
    }

    fn on_initialized(&mut self, sys: &PackContext) {
        if !self.level.is_enabled() {
            return;
        }
        eprintln!(
            "Initialization: time={:.3}s overlap={:.4e} restraints={:.4e}",
            self.elapsed(),
            sys.fdist,
            sys.frest
        );
    }

    fn on_phase_start(&mut self, info: &PhaseInfo) {
        if !self.level.is_enabled() {
            return;
        }
        self.phase_start = Some(Instant::now());
        let desc = match info.molecule_type {
            Some(itype) => format!("type {itype} compaction"),
            None => "all-type optimization".to_string(),
        };
        eprintln!("Phase {}/{}: {desc}", info.phase + 1, info.total_phases);
        if self.level >= LogLevel::Progress {
            if self.level >= LogLevel::Verbose {
                eprintln!(
                    "{:>8} {:>14} {:>14} {:>10} {:>10} {:>10}",
                    "Step", "Overlap", "Restraint", "Improve%", "RadScale", "Time"
                );
            } else {
                eprintln!(
                    "{:>8} {:>14} {:>14} {:>10} {:>10}",
                    "Step", "Overlap", "Restraint", "Improve%", "Time"
                );
            }
        }
    }

    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        if self.level < LogLevel::Progress || !info.loop_idx.is_multiple_of(self.every) {
            return;
        }
        if self.level >= LogLevel::Verbose {
            eprintln!(
                "{:>8} {:>14.6e} {:>14.6e} {:>10.3} {:>10.4} {:>10.3}",
                info.loop_idx + 1,
                info.fdist,
                info.frest,
                info.improvement_pct,
                info.radscale,
                self.elapsed(),
            );
        } else {
            eprintln!(
                "{:>8} {:>14.6e} {:>14.6e} {:>10.3} {:>10.3}",
                info.loop_idx + 1,
                info.fdist,
                info.frest,
                info.improvement_pct,
                self.elapsed(),
            );
        }
    }

    fn on_phase_end(&mut self, info: &PhaseInfo, report: &PhaseReport) {
        if !self.level.is_enabled() {
            return;
        }
        eprintln!(
            "Phase {}/{} summary: steps={} converged={} overlap={:.6e} restraints={:.6e} time={:.3}s",
            info.phase + 1,
            info.total_phases,
            report.iterations,
            report.converged,
            report.fdist,
            report.frest,
            self.phase_elapsed(),
        );
    }

    fn on_finish(&mut self, sys: &PackContext) {
        if !self.level.is_enabled() {
            return;
        }
        let converged = sys.fdist < self.precision && sys.frest < self.precision;
        eprintln!(
            "Final summary: converged={} overlap={:.6e} restraints={:.6e} elapsed={:.3}s",
            converged,
            sys.fdist,
            sys.frest,
            self.elapsed(),
        );
    }
}

// ── EarlyStopHandler ──────────────────────────────────────────────────────────

/// Requests early termination when improvement stalls.
///
/// Tracks `fdist + frest` total violation. After `warmup` iterations, if the
/// relative improvement drops below `threshold` for `patience` consecutive
/// steps, sets the stop flag.
///
/// Installed by the engine lifecycle according to [`LogSpec`](crate::entry::LogSpec).
pub struct EarlyStopHandler {
    /// Relative improvement threshold.
    pub threshold: F,
    /// Iterations to skip before tracking. Default: `5`.
    pub warmup: usize,
    /// Consecutive stall iterations before stopping. Default: `3`.
    pub patience: usize,
    prev_violation: F,
    stall_count: usize,
    stop: bool,
}

impl EarlyStopHandler {
    pub fn new(threshold: F) -> Self {
        Self {
            threshold,
            warmup: 5,
            patience: 3,
            prev_violation: F::INFINITY,
            stall_count: 0,
            stop: false,
        }
    }

    pub fn with_warmup(mut self, warmup: usize) -> Self {
        self.warmup = warmup;
        self
    }

    pub fn with_patience(mut self, patience: usize) -> Self {
        self.patience = patience.max(1);
        self
    }
}

impl Default for EarlyStopHandler {
    /// Default is intentionally conservative (effectively disabled),
    /// so Packmol-sized examples are not stopped before convergence.
    fn default() -> Self {
        Self::new(F::NEG_INFINITY)
    }
}

impl Handler for EarlyStopHandler {
    fn on_initialized(&mut self, _sys: &PackContext) {
        self.prev_violation = F::INFINITY;
        self.stall_count = 0;
        self.stop = false;
    }

    fn on_phase_start(&mut self, _info: &PhaseInfo) {
        self.prev_violation = F::INFINITY;
        self.stall_count = 0;
        self.stop = false;
    }

    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        let v = info.fdist + info.frest;
        if info.loop_idx <= self.warmup {
            self.prev_violation = v;
            return;
        }
        let rel_change = if self.prev_violation > 0.0 {
            (self.prev_violation - v) / self.prev_violation
        } else if v < objective_small_floor() {
            1.0 // already converged
        } else {
            F::INFINITY
        };

        if rel_change < self.threshold {
            self.stall_count += 1;
            if self.stall_count >= self.patience {
                log::debug!(
                    "EarlyStop: stalled for {} iters (rel_change={:.2e})",
                    self.stall_count,
                    rel_change
                );
                self.stop = true;
            }
        } else {
            self.stall_count = 0;
        }
        self.prev_violation = v;
    }

    fn should_stop(&self) -> bool {
        self.stop
    }
}
