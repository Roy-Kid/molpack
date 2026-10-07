//! The [`Callback`] trait the packing lifecycle calls at its start, steps,
//! phases, stages and finish, and the built-in callbacks (progress, LAMMPS-style
//! log, early stop, XYZ trajectory).

use molrs::core::SimBox;
use molrs::op::F;

use std::time::Instant;

use crate::context::PackContext;
use crate::outcome::StageOutcome;

// ── Info structs ─────────────────────────────────────────────────────────────

/// Identifies the stage a callback comes from.
///
/// The field shape follows [`PhaseInfo`] on purpose — a stage is to a run
/// what a phase is to the GENCAN loop, so the two identities read the same
/// way — and, like `PhaseInfo`, this is a plain `Copy` record a caller may
/// build by literal (a callback test drives the two stage hooks with one).
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
/// Passed to [`Callback::on_phase_end`]. Fields mirror the counters
/// reported by `ProgressCallback::on_phase_start`; when no summary is
/// needed (e.g. early termination), callers may omit the callback hook.
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
/// [`State`](crate::State).
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
    /// longest chain's steps + 1)`; see the growth driver (`grow/driver.rs`).
    pub max_loops: usize,
    /// Current phase info.
    pub phase: PhaseInfo,
    /// Max inter-molecular overlap violation (0.0 = no overlap).
    pub fdist: F,
    /// Max constraint violation (0.0 = all constraints satisfied).
    pub frest: F,
    /// GENCAN: the objective at the user's radii — Packmol's `fx` after a
    /// loop (packmol.f90:833-841), the value `fimp` and `bestf` are measured
    /// on. Growth: a constructive `0.0`, like `fdist` / `frest`.
    pub f: F,
    /// Improvement from last iteration, as percentage (positive = improving).
    pub improvement_pct: F,
    /// GENCAN: current radius scaling factor (starts at discale, decays to
    /// 1.0). Growth: the dimensionless hard-core scale — the factor multiplying
    /// the pair contact distance a placement must clear (`1.0` = full declared
    /// contact) — which starts at 1.0 and steps down by
    /// [`GrowConfig::SOFTEN_RUNG`](crate::grow::GrowConfig::SOFTEN_RUNG) to the
    /// floor set by
    /// [`GrowConfig::with_min_hard_scale`](crate::grow::GrowConfig::with_min_hard_scale).
    pub radscale: F,
    /// Convergence precision target.
    pub precision: F,
}

// ── Trait ─────────────────────────────────────────────────────────────────────

/// The hooks called by the [`PackEngine`](crate::PackEngine) lifecycle during packing.
pub trait Callback: Send {
    /// Called immediately at the start of [`run`][crate::PackEngine::run],
    /// before any computation. Use this for immediate user feedback.
    fn on_start(&mut self, _ntotat: usize, _ntotmol: usize) {}

    /// Called once after initialization completes, with valid `xcart` positions.
    /// Use this to write the initial conformation (e.g. `XyzTrajectoryCallback`).
    fn on_initialized(&mut self, _sys: &PackContext) {}

    /// Called after each outer optimization loop iteration.
    fn on_step(&mut self, info: &StepInfo, sys: &PackContext);

    /// Called at the start of each packing phase (per-type and all-types).
    /// Allows stateful callbacks to reset between phases.
    fn on_phase_start(&mut self, _info: &PhaseInfo) {}

    /// Called once after the packing loop finishes (convergence or max loops).
    fn on_finish(&mut self, _sys: &PackContext) {}

    /// Return `true` to request early termination of the packing loop.
    fn should_stop(&self) -> bool {
        false
    }

    /// Called at the end of each packing phase, with a summary report.
    ///
    /// Default: no-op. Paired with [`on_phase_start`] for symmetric
    /// setup / teardown hooks.
    ///
    /// [`on_phase_start`]: Callback::on_phase_start
    fn on_phase_end(&mut self, _info: &PhaseInfo, _report: &PhaseReport) {}

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
    /// [`StageOutcome`] deliberately carries no verdict. A callback that wants
    /// the run's violation maxima reads them off `sys` — `sys.fdist` and
    /// `sys.frest`, the shared objective's numbers on the post-run state,
    /// exactly as [`on_finish`] does. Same shape, same authority.
    ///
    /// [`on_finish`]: Callback::on_finish
    fn on_stage_end(&mut self, _info: &StageInfo, _outcome: &StageOutcome, _sys: &PackContext) {}
}

// ── XyzTrajectoryCallback ────────────────────────────────────────────────────────────────

/// Writes packing snapshots as a multi-frame extended XYZ trajectory.
///
/// Writes a frame on every `every`-th step (starting from step 0).
/// No automatic initial or final writes. Each snapshot is the
/// coordinates-only frame the final result falls back to (`atoms`: `id`,
/// 1-based `mol_id`, `x`/`y`/`z`, `element`) with the step in the frame's
/// `step` meta key, written by molrs's extended XYZ writer
/// (`molrs::io::xyz::XyzWriter`); molpack keeps no XYZ format code of its
/// own. Needs the `io` feature.
#[cfg(feature = "io")]
pub struct XyzTrajectoryCallback {
    path: std::path::PathBuf,
    /// Write every `n` steps (must be >= 1).
    every: usize,
    file: Option<std::io::BufWriter<std::fs::File>>,
}

#[cfg(feature = "io")]
impl XyzTrajectoryCallback {
    /// Create a new XYZ trajectory callback.
    ///
    /// `every` controls writing frequency: a frame is written on every step
    /// where `loop_idx % every == 0` (step 0 is always included).
    /// `every` must be >= 1.
    pub fn new(path: impl Into<std::path::PathBuf>, every: usize) -> Self {
        assert!(every >= 1, "every must be >= 1");
        Self {
            path: path.into(),
            every,
            file: None,
        }
    }

    fn open(&mut self) {
        if self.file.is_some() {
            return;
        }
        match std::fs::File::create(&self.path) {
            Ok(f) => self.file = Some(std::io::BufWriter::new(f)),
            Err(e) => log::warn!(
                "XyzTrajectoryCallback: cannot open {}: {e}",
                self.path.display()
            ),
        }
    }

    /// The snapshot as a frame: the coordinates-only frame of `xcart`, one
    /// `atoms` row per entry, with the step in its meta.
    fn snapshot(step: usize, sys: &PackContext) -> molrs::core::Frame {
        let elements = (0..sys.xcart.len())
            .map(|i| {
                sys.elements
                    .get(i)
                    .and_then(|e| *e)
                    .map_or("X", |e| e.symbol())
                    .to_string()
            })
            .collect();
        let groups = (0..sys.ntype_with_fixed).map(|t| (sys.natoms[t], sys.nmols[t]));
        let mut frame =
            crate::assemble::coords_frame(&sys.xcart, elements, crate::assemble::mol_ids(groups));
        frame.meta.insert("step", step as u64);
        frame
    }

    fn write_snapshot(&mut self, step: usize, sys: &PackContext) {
        use molrs::io::writer::{FrameWriter, Writer};

        self.open();
        let frame = Self::snapshot(step, sys);
        let Some(ref mut w) = self.file else { return };
        let written = molrs::io::xyz::XyzWriter::new(&mut *w)
            .write(&frame)
            .and_then(|()| std::io::Write::flush(w));
        if let Err(e) = written {
            log::warn!(
                "XyzTrajectoryCallback: writing {}: {e}",
                self.path.display()
            );
        }
    }
}

#[cfg(feature = "io")]
impl Callback for XyzTrajectoryCallback {
    fn on_step(&mut self, info: &StepInfo, sys: &PackContext) {
        if info.loop_idx.is_multiple_of(self.every) {
            self.write_snapshot(info.loop_idx, sys);
        }
    }
}

// ── ProgressCallback ───────────────────────────────────────────────────────────

/// Prints human-readable progress lines to `stderr`.
///
/// Attach it with [`PackEngine::with_callback`](crate::PackEngine::with_callback);
/// the screen log the engine installs from its log level is
/// [`LammpsLogCallback`].
pub struct ProgressCallback {
    start: Option<Instant>,
}

impl ProgressCallback {
    pub fn new() -> Self {
        Self { start: None }
    }
}

impl Default for ProgressCallback {
    fn default() -> Self {
        Self::new()
    }
}

impl Callback for ProgressCallback {
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

// ── LammpsLogCallback ─────────────────────────────────────────────────────────

/// LAMMPS-style screen log for packing runs.
///
/// Most users should enable this through
/// [`PackEngine::with_log_level`][crate::PackEngine::with_log_level]
/// or [`with_log_frequency`][crate::PackEngine::with_log_frequency]
/// instead of attaching the callback manually.
pub struct LammpsLogCallback {
    level: LogLevel,
    every: usize,
    tolerance: F,
    precision: F,
    seed: u64,
    max_loops: usize,
    ntypes: usize,
    cell: Option<SimBox>,
    start: Option<Instant>,
    phase_start: Option<Instant>,
}

impl LammpsLogCallback {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        level: LogLevel,
        every: usize,
        tolerance: F,
        precision: F,
        seed: u64,
        max_loops: usize,
        ntypes: usize,
        cell: Option<SimBox>,
    ) -> Self {
        Self {
            level,
            every: every.max(1),
            tolerance,
            precision,
            seed,
            max_loops,
            ntypes,
            cell,
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

impl Callback for LammpsLogCallback {
    fn on_start(&mut self, ntotat: usize, ntotmol: usize) {
        if !self.level.is_enabled() {
            return;
        }

        let now = Instant::now();
        self.start = Some(now);
        self.phase_start = Some(now);

        eprintln!("molpack screen log");
        eprintln!("System information:");
        eprintln!("  molecule types = {}", self.ntypes);
        eprintln!("  molecules      = {ntotmol}");
        eprintln!("  atoms          = {ntotat}");
        eprintln!("Settings:");
        eprintln!("  tolerance      = {:.6} A", self.tolerance);
        eprintln!("  precision      = {:.6}", self.precision);
        eprintln!("  seed           = {}", self.seed);
        eprintln!("  nloop          = {}", self.max_loops);
        match &self.cell {
            Some(cell) => {
                let (o, h) = (cell.origin_view(), cell.h_view());
                eprintln!(
                    "  cell origin    = [{:.6}, {:.6}, {:.6}]  pbc={:?}",
                    o[0],
                    o[1],
                    o[2],
                    cell.pbc()
                );
                for (k, axis) in ["a", "b", "c"].iter().enumerate() {
                    eprintln!(
                        "  cell {axis}         = [{:.6}, {:.6}, {:.6}]",
                        h[[0, k]],
                        h[[1, k]],
                        h[[2, k]]
                    );
                }
            }
            None => eprintln!("  cell           = none"),
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

// ── EarlyStopCallback ──────────────────────────────────────────────────────────

/// Ends a GENCAN phase whose best objective has stopped improving.
///
/// Packmol has no early stop — a phase runs until it converges or reaches
/// `nloop` — so this callback is molpack's, but it is phrased entirely in
/// Packmol's own quantities (`app/packmol.f90`):
///
/// * the objective is `fx`, the function value at the user's radii after each
///   loop ([`StepInfo::f`]), and `bestf` is its per-phase minimum
///   (packmol.f90:808, 869, 894);
/// * improvement is Packmol's `fimprov`, `-100 * (fx - bestf) / bestf`, in
///   percent, with `bestf == 0` counting as 100 % (packmol.f90:844-845) —
///   here measured across a window: `bestf` now against `bestf` `patience`
///   loops ago;
/// * a phase is judged only once `radscale` has come down to 1.0. That is
///   where Packmol itself starts treating poor progress as a stall (movebad
///   fires only at `radscale == 1` and `fimp <= 10`, packmol.f90:815); above
///   1.0 its own schedule is still shrinking the radii (packmol.f90:940-947);
/// * the default threshold is that same 10 %.
///
/// So the rule is: at the user's radii, if `bestf` has improved by less than
/// `threshold_pct` over the last `patience` loops, end the phase. `bestf`
/// rather than `fx` because movebad kicks `fx` up on purpose.
///
/// The stop ends the current phase only; later phases and later pipeline
/// stages still run. [`GencanPack`](crate::GencanPack) installs one by
/// default ([`Default`] values); `with_early_stop` replaces or removes it. A
/// run whose final phase was stopped is not converged:
/// [`State::converged`](crate::State) reports `false`, as Packmol's
/// "maximum number of GENCAN loops achieved" would.
#[derive(Debug, Clone)]
pub struct EarlyStopCallback {
    /// Minimum improvement of `bestf` across `patience` loops, in percent
    /// (Packmol's `fimprov` units). Default: `10.0`, Packmol's movebad
    /// threshold.
    pub threshold_pct: F,
    /// Window length, in loops at `radscale == 1`. Default: `10`.
    pub patience: usize,
    /// Per-phase best `fx` (Packmol's `bestf`).
    bestf: F,
    /// `bestf` after each loop judged so far in this phase.
    window: Vec<F>,
    stop: bool,
}

impl EarlyStopCallback {
    pub fn new(threshold_pct: F) -> Self {
        Self {
            threshold_pct,
            patience: 10,
            bestf: F::INFINITY,
            window: Vec::new(),
            stop: false,
        }
    }

    pub fn with_patience(mut self, patience: usize) -> Self {
        self.patience = patience.max(1);
        self
    }

    fn reset(&mut self) {
        self.bestf = F::INFINITY;
        self.window.clear();
        self.stop = false;
    }

    /// One outer loop: `fx` at the user's radii, and the radius scale the
    /// loop ran with.
    fn observe(&mut self, fx: F, radscale: F) {
        self.bestf = self.bestf.min(fx);
        if radscale != 1.0 {
            return;
        }
        self.window.push(self.bestf);
        let n = self.window.len();
        if n <= self.patience {
            return;
        }
        let before = self.window[n - 1 - self.patience];
        let fimprov = if before > 0.0 {
            -100.0 * (self.bestf - before) / before
        } else {
            100.0
        };
        if fimprov < self.threshold_pct {
            log::debug!(
                "EarlyStop: bestf {before:.3e} -> {:.3e} ({fimprov:.2} %) over {} loops",
                self.bestf,
                self.patience
            );
            self.stop = true;
        }
    }
}

impl Default for EarlyStopCallback {
    fn default() -> Self {
        Self::new(10.0)
    }
}

impl Callback for EarlyStopCallback {
    fn on_initialized(&mut self, _sys: &PackContext) {
        self.reset();
    }

    fn on_phase_start(&mut self, _info: &PhaseInfo) {
        self.reset();
    }

    /// The stop is phase-scoped: GENCAN polls `should_stop` right after
    /// `on_step` and ends the phase, then calls this. Clearing here keeps the
    /// flag from reaching the pipeline's between-stage check, where it would
    /// also cancel every later stage.
    fn on_phase_end(&mut self, _info: &PhaseInfo, _report: &PhaseReport) {
        self.reset();
    }

    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        self.observe(info.f, info.radscale);
    }

    fn should_stop(&self) -> bool {
        self.stop
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(h: &mut EarlyStopCallback, fx: &[F], radscale: F) -> Option<usize> {
        for (i, &f) in fx.iter().enumerate() {
            h.observe(f, radscale);
            if h.should_stop() {
                return Some(i);
            }
        }
        None
    }

    #[test]
    fn early_stop_fires_after_patience_loops_of_a_plateau() {
        // The first judged loop only seeds the window: a flat bestf is
        // stopped on loop `patience`.
        assert_eq!(
            run(&mut EarlyStopCallback::default(), &[3.0; 30], 1.0),
            Some(10)
        );
    }

    #[test]
    fn early_stop_never_judges_above_the_user_radii() {
        // Packmol's own radscale schedule is still running: never stop there.
        assert_eq!(
            run(&mut EarlyStopCallback::default(), &[3.0; 60], 1.1),
            None
        );
    }

    #[test]
    fn early_stop_keeps_a_phase_whose_bestf_falls_ten_percent_per_window() {
        // fx is kicked up every 5 loops (movebad) but bestf still falls
        // 2 % a loop, ~18 % per 10-loop window.
        let fx: Vec<F> = (0..60)
            .map(|i| {
                if i % 5 == 4 {
                    1e3
                } else {
                    100.0 * (0.98 as F).powi(i)
                }
            })
            .collect();
        assert_eq!(run(&mut EarlyStopCallback::default(), &fx, 1.0), None);
    }

    #[test]
    fn early_stop_carries_bestf_from_the_scaled_loops() {
        // bestf is per phase, as in Packmol, so a good point found at
        // radscale > 1 is what the first judged loop is measured against.
        let mut h = EarlyStopCallback::default();
        assert_eq!(run(&mut h, &[1.0], 1.1), None);
        assert_eq!(run(&mut h, &[5.0; 20], 1.0), Some(10));
    }

    #[test]
    fn early_stop_resets_at_phase_end_and_never_fires_on_zero() {
        let mut h = EarlyStopCallback::default();
        assert_eq!(run(&mut h, &[3.0; 11], 1.0), Some(10));
        let info = PhaseInfo {
            phase: 0,
            total_phases: 1,
            molecule_type: None,
        };
        h.on_phase_end(&info, &PhaseReport::default());
        assert!(!h.should_stop(), "the stop must not outlive its phase");
        assert_eq!(run(&mut h, &[0.0; 40], 1.0), None);
    }
}
