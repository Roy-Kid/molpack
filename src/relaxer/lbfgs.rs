//! Force-field geometry relaxer (`ff` feature).
//!
//! [`LBFGSRelaxer`] is a [`Relaxer`] that relaxes a flexible molecule's
//! **internal geometry** during packing by energy-minimizing it under a molrs
//! force-field [`Potential`] with the limited-memory BFGS optimizer from
//! `molrs::ff::optimize`. It complements [`TorsionMcRelaxer`] (stochastic
//! torsion sampling) with a deterministic gradient-based relaxation.
//!
//! Like every relaxer it operates on the **reference geometry** axis (it changes
//! the molecule's shape), not the objective axis. The relaxed conformer is only
//! swapped in when it does **not** worsen the packer objective, so the relaxer
//! is a strict non-harm operation with respect to GENCAN convergence.
//!
//! The force field itself is supplied by the caller as an `Arc<dyn Potential>`:
//! build it from a real molecule via `molrs::ff::typifier::mmff::MMFFTypifier`,
//! or hand-assemble one from the `molrs::ff::potential` kernels. Keeping the
//! potential caller-supplied means molpack stays agnostic to the chemistry and
//! the user controls the force field.
//!
//! [`TorsionMcRelaxer`]: crate::relaxer::TorsionMcRelaxer

use std::sync::Arc;

use molrs::Frame;
use molrs::ff::ForceField;
use molrs::ff::optimize::{LBFGS, LbfgsConfig};
use molrs::ff::potential::{Potential, intramolecular_pairs};
use molrs::types::F;
use rand::RngCore;

use super::{Relaxer, RelaxerRunner, recenter};

// ── LBFGSRelaxer ────────────────────────────────────────────────────────

/// Force-field energy-minimization relaxer. Stored on `Target` (immutable config).
///
/// Analogous to [`TorsionMcRelaxer`][crate::relaxer::TorsionMcRelaxer] — a
/// concrete [`Relaxer`] — but driven by a molrs force-field [`Potential`] and the
/// L-BFGS minimizer instead of Monte Carlo torsion moves.
///
/// The potential source is cloned into every spawned runner, so a
/// `LBFGSRelaxer` is itself `Clone` despite wrapping a trait object.
#[derive(Clone)]
pub struct LBFGSRelaxer {
    /// Where the molecule potential comes from (pre-built, or a force field
    /// compiled lazily against the molecule frame at `spawn`).
    source: RelaxerSource,
    /// L-BFGS convergence + trust-region controls (defaults from molrs).
    cfg: LbfgsConfig,
}

/// The two ways a [`LBFGSRelaxer`] obtains its [`Potential`].
#[derive(Clone)]
enum RelaxerSource {
    /// A pre-built potential, used directly — chemistry-agnostic, the original
    /// path (handy for tests and callers who assemble their own kernels).
    Potential(Arc<dyn Potential>),
    /// A molrs [`ForceField`] compiled to a potential **lazily** at `spawn`,
    /// once the molecule's frame is known: build the neighbour list
    /// (`intramolecular_pairs`) on the frame, then `to_potentials`. This is the
    /// `molpack.LBFGSRelaxer(ff)` path.
    ForceField(ForceField),
}

impl std::fmt::Debug for LBFGSRelaxer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // `dyn Potential` is not `Debug`; surface only the optimizer config.
        f.debug_struct("LBFGSRelaxer")
            .field("cfg", &self.cfg)
            .finish_non_exhaustive()
    }
}

impl LBFGSRelaxer {
    /// Create a relaxer bound to `potential`, using molrs's default L-BFGS
    /// configuration (`fmax = 0.05`, `max_steps = 500`, `max_step = 0.2`,
    /// `memory = 8`).
    ///
    /// The potential must compute energy + forces for exactly the atoms of the
    /// target molecule, in the same order as the target's reference coords.
    pub fn new(potential: Arc<dyn Potential>) -> Self {
        Self {
            source: RelaxerSource::Potential(potential),
            cfg: LbfgsConfig::default(),
        }
    }

    /// Create a relaxer from a molrs [`ForceField`]. The potential is compiled
    /// **lazily** — the first time the relaxer is spawned with the molecule's
    /// frame during packing — by building the neighbour list
    /// (`intramolecular_pairs`) and calling `to_potentials`. This is the
    /// binding-facing constructor behind `molpack.LBFGSRelaxer(ff)`: the relaxer
    /// need not know the molecule (or carry a `pairs` block) until packing starts.
    pub fn from_forcefield(ff: ForceField) -> Self {
        Self {
            source: RelaxerSource::ForceField(ff),
            cfg: LbfgsConfig::default(),
        }
    }

    /// Set the convergence tolerance: stop when the max per-atom force magnitude
    /// drops below `fmax` (kcal/mol/Å).
    pub fn with_fmax(mut self, fmax: F) -> Self {
        self.cfg.fmax = fmax;
        self
    }

    /// Set the outer-iteration cap per relaxation call.
    pub fn with_max_steps(mut self, max_steps: usize) -> Self {
        self.cfg.max_steps = max_steps;
        self
    }

    /// Set the per-step trust region in Å. `f64::INFINITY` disables it.
    pub fn with_max_step(mut self, max_step: F) -> Self {
        self.cfg.max_step = max_step;
        self
    }

    /// Set the L-BFGS correction-pair history size.
    pub fn with_memory(mut self, memory: usize) -> Self {
        self.cfg.memory = memory;
        self
    }
}

impl Relaxer for LBFGSRelaxer {
    fn spawn(&self, frame: Option<&Frame>, _ref_coords: &[[F; 3]]) -> Box<dyn RelaxerRunner> {
        let potential = match &self.source {
            RelaxerSource::Potential(p) => Some(Arc::clone(p)),
            RelaxerSource::ForceField(ff) => compile_forcefield(ff, frame),
        };
        Box::new(LBFGSRelaxerRunner {
            potential,
            cfg: self.cfg,
            attempts: 0,
            accepts: 0,
        })
    }
}

/// Compile a [`ForceField`] into a molecule-bound [`Potential`].
///
/// The consumer owns the neighbour list: build it on the frame with
/// `intramolecular_pairs`, insert the `pairs` block, then `to_potentials`.
/// Returns `None` — and the runner becomes a no-op — when there is no frame
/// (coordinate-only target) or compilation fails (logged), so a misconfigured
/// relaxer degrades to "leave the geometry unchanged" rather than aborting the pack.
fn compile_forcefield(ff: &ForceField, frame: Option<&Frame>) -> Option<Arc<dyn Potential>> {
    let mut frame = frame?.clone();
    let pairs = intramolecular_pairs(&frame);
    frame.insert("pairs", pairs);
    match ff.to_potentials(&frame) {
        Ok(pots) => Some(Arc::new(pots)),
        Err(e) => {
            log::warn!(
                "LBFGSRelaxer: compiling the force field to a potential failed, \
                 relaxer disabled for this molecule: {e}"
            );
            None
        }
    }
}

// ── LBFGSRelaxerRunner ──────────────────────────────────────────────────

/// Runtime state for a [`LBFGSRelaxer`]. Created by `spawn`, used in `pack()`.
///
/// `potential` is `None` when the force field could not be compiled (no frame or
/// a `to_potentials` error); such a runner is a no-op.
struct LBFGSRelaxerRunner {
    potential: Option<Arc<dyn Potential>>,
    cfg: LbfgsConfig,
    attempts: usize,
    accepts: usize,
}

impl RelaxerRunner for LBFGSRelaxerRunner {
    fn on_iter(
        &mut self,
        coords: &[[F; 3]],
        f_current: F,
        evaluate: &mut dyn FnMut(&[[F; 3]]) -> F,
        _rng: &mut dyn RngCore,
    ) -> Option<Vec<[F; 3]>> {
        // No compiled potential (no frame, or force-field compilation failed) →
        // nothing to relax against; leave the geometry unchanged.
        let potential = self.potential.as_ref()?;

        // Flatten to molrs's `[x0,y0,z0, x1,y1,z1, …]` 3N layout and relax the
        // internal geometry under the force field.
        let mut flat: Vec<F> = coords.iter().flat_map(|p| *p).collect();
        // `run` only errors on a non-3N buffer, which cannot happen here, but
        // surface it rather than silently no-op'ing in case molrs grows new
        // error conditions (e.g. a non-finite energy).
        let report = match LBFGS::new(&**potential, self.cfg).run(&mut flat) {
            Ok(report) => report,
            Err(e) => {
                log::debug!(
                    "LBFGSRelaxer: L-BFGS minimization failed, leaving geometry unchanged: {e}"
                );
                return None;
            }
        };

        // Once a conformer has reached the force-field minimum, later packing
        // iterations only move its placement (COM + Euler), never its reference
        // geometry — so L-BFGS converges in ≤1 step with nothing to apply. Skip
        // the rebuild and the full-objective packing gate in that common case
        // (it is not a relaxation attempt, so it does not count toward the rate).
        if report.converged && report.n_steps <= 1 {
            return None;
        }
        self.attempts += 1;

        let mut relaxed: Vec<[F; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();

        // Keep the centroid at the origin: the packer places each molecule by
        // COM + Euler rotation applied to the reference coords, so the reference
        // conformer is centered by convention (matches `TorsionMcRelaxer`).
        recenter(&mut relaxed);

        // Packing gate: only accept the relaxed conformer when it does not
        // worsen the packer objective (overlap + restraints). This guarantees
        // the relaxer never destabilizes GENCAN's convergence.
        if evaluate(&relaxed) <= f_current {
            self.accepts += 1;
            Some(relaxed)
        } else {
            None
        }
    }

    fn acceptance_rate(&self) -> F {
        if self.attempts == 0 {
            0.0
        } else {
            self.accepts as F / self.attempts as F
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Two atoms joined by a single harmonic bond `E = ½k(r − r0)²` along x.
    /// Mirrors the molrs optimizer test potential so we can relax without
    /// depending on MMFF typification.
    struct HarmonicBond {
        k: F,
        r0: F,
    }

    impl Potential for HarmonicBond {
        fn calc_energy_forces(&self, coords: &[F]) -> (F, Vec<F>) {
            let d = [
                coords[3] - coords[0],
                coords[4] - coords[1],
                coords[5] - coords[2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let e = 0.5 * self.k * (r - self.r0) * (r - self.r0);
            let mut f = vec![0.0; 6];
            if r > 1e-12 {
                let coeff = self.k * (r - self.r0) / r; // dE/dr / r
                for i in 0..3 {
                    let fi = coeff * d[i];
                    f[i] = fi; // force on atom 0
                    f[3 + i] = -fi; // force on atom 1
                }
            }
            (e, f)
        }
    }

    /// A 1.5 Å bond along x, centered at the origin (r0 = 1.0 → stretched).
    fn stretched() -> Vec<[F; 3]> {
        vec![[-0.75, 0.0, 0.0], [0.75, 0.0, 0.0]]
    }

    fn bond_len(c: &[[F; 3]]) -> F {
        let d = [c[1][0] - c[0][0], c[1][1] - c[0][1], c[1][2] - c[0][2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }

    fn relaxer() -> LBFGSRelaxer {
        let pot: Arc<dyn Potential> = Arc::new(HarmonicBond { k: 100.0, r0: 1.0 });
        LBFGSRelaxer::new(pot)
    }

    #[test]
    fn relaxes_stretched_bond_to_equilibrium() {
        let coords = stretched();
        let mut runner = relaxer().spawn(None, &coords);
        let mut rng = rand::rng();

        // `evaluate` returns 0 → packing never worsens → accept the relaxation.
        let out = runner.on_iter(&coords, 0.0, &mut |_| 0.0, &mut rng);
        let relaxed = out.expect("relaxed conformer should be accepted");

        let r = bond_len(&relaxed);
        assert!(
            (r - 1.0).abs() < 1e-3,
            "bond should relax to r0 = 1.0, got {r}"
        );
    }

    #[test]
    fn rejects_when_packing_worsens() {
        let coords = stretched();
        let mut runner = relaxer().spawn(None, &coords);
        let mut rng = rand::rng();

        // The relaxed conformer scores far worse for packing than f_current →
        // the gate must reject it (return None), leaving the geometry untouched.
        let out = runner.on_iter(&coords, 1.0, &mut |_| 1.0e9, &mut rng);
        assert!(out.is_none(), "must reject when packing objective worsens");
    }

    #[test]
    fn accepts_when_packing_neutral_or_better() {
        let coords = stretched();
        let mut runner = relaxer().spawn(None, &coords);
        let mut rng = rand::rng();

        // f_trial (0.5) <= f_current (0.5) → accept.
        let out = runner.on_iter(&coords, 0.5, &mut |_| 0.5, &mut rng);
        assert!(out.is_some(), "must accept when packing does not worsen");
    }

    #[test]
    fn returned_coords_are_centered() {
        let coords = stretched();
        let mut runner = relaxer().spawn(None, &coords);
        let mut rng = rand::rng();

        let relaxed = runner
            .on_iter(&coords, 0.0, &mut |_| 0.0, &mut rng)
            .expect("accepted");
        let n = relaxed.len() as F;
        let cx: F = relaxed.iter().map(|p| p[0]).sum::<F>() / n;
        let cy: F = relaxed.iter().map(|p| p[1]).sum::<F>() / n;
        let cz: F = relaxed.iter().map(|p| p[2]).sum::<F>() / n;
        assert!(cx.abs() < 1e-9 && cy.abs() < 1e-9 && cz.abs() < 1e-9);
    }

    #[test]
    fn acceptance_rate_tracks_accepts() {
        let coords = stretched();
        let mut runner = relaxer().spawn(None, &coords);
        let mut rng = rand::rng();

        assert_eq!(runner.acceptance_rate(), 0.0);
        let _ = runner.on_iter(&coords, 0.0, &mut |_| 0.0, &mut rng);
        assert!((runner.acceptance_rate() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn already_minimized_geometry_is_skipped() {
        // A bond already at r0 → L-BFGS converges in ≤1 step with nothing to
        // move, so on_iter short-circuits to None *without* running the packing
        // gate (the `evaluate` closure panics if reached) and does not count the
        // call as a relaxation attempt.
        let pot: Arc<dyn Potential> = Arc::new(HarmonicBond { k: 100.0, r0: 1.0 });
        let mut runner = LBFGSRelaxer::new(pot).spawn(None, &[]);
        let mut rng = rand::rng();
        let at_min: Vec<[F; 3]> = vec![[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]];

        let out = runner.on_iter(
            &at_min,
            0.0,
            &mut |_| panic!("packing gate must be skipped for a no-op relaxation"),
            &mut rng,
        );
        assert!(out.is_none(), "already-minimized geometry yields no change");
        assert_eq!(
            runner.acceptance_rate(),
            0.0,
            "a skipped no-op must not count as a relaxation attempt"
        );
    }

    #[test]
    fn config_builders_apply() {
        let r = relaxer()
            .with_fmax(0.01)
            .with_max_steps(50)
            .with_max_step(0.1)
            .with_memory(4);
        assert_eq!(r.cfg.fmax, 0.01);
        assert_eq!(r.cfg.max_steps, 50);
        assert_eq!(r.cfg.max_step, 0.1);
        assert_eq!(r.cfg.memory, 4);
    }

    #[test]
    fn from_forcefield_compiles_lazily_and_relaxes() {
        use molrs::store::block::Block;
        use molrs::types::U;
        use ndarray::Array1;

        // A single harmonic bond A-A with equilibrium length r0 = 1.0.
        let mut ff = ForceField::new("test");
        ff.def_bondstyle("harmonic")
            .def_type("A-A", &[("k", 100.0), ("r0", 1.0)]);

        // A 2-atom molecule carrying that bond, stretched to 1.5 Å.
        let mut frame = Frame::new();
        let mut atoms = Block::new();
        atoms
            .insert("x", Array1::from_vec(vec![-0.75 as F, 0.75]).into_dyn())
            .unwrap();
        for col in ["y", "z"] {
            atoms
                .insert(col, Array1::from_vec(vec![0.0 as F, 0.0]).into_dyn())
                .unwrap();
        }
        frame.insert("atoms", atoms);
        let mut bonds = Block::new();
        bonds
            .insert("atomi", Array1::from_vec(vec![0 as U]).into_dyn())
            .unwrap();
        bonds
            .insert("atomj", Array1::from_vec(vec![1 as U]).into_dyn())
            .unwrap();
        bonds
            .insert("type", Array1::from_vec(vec!["A-A".to_string()]).into_dyn())
            .unwrap();
        frame.insert("bonds", bonds);

        // `from_forcefield` holds only the FF; the potential is compiled at spawn
        // against the molecule frame.
        let coords = stretched();
        let mut runner = LBFGSRelaxer::from_forcefield(ff).spawn(Some(&frame), &coords);
        let mut rng = rand::rng();
        let relaxed = runner
            .on_iter(&coords, 0.0, &mut |_| 0.0, &mut rng)
            .expect("force-field relaxer should compile and relax");
        assert!(
            (bond_len(&relaxed) - 1.0).abs() < 1e-3,
            "bond should relax to r0 = 1.0, got {}",
            bond_len(&relaxed)
        );
    }

    #[test]
    fn from_forcefield_with_nonbonded_compiles_and_relaxes() {
        use molrs::store::block::Block;
        use molrs::types::U;
        use ndarray::Array1;

        // A force field with bonded *and* non-bonded styles: harmonic bond +
        // per-atom lj/cut + coul/cut. Exercises the full compile path through
        // `from_forcefield` (neighbour list + per-atom pair kernels).
        let mut ff = ForceField::new("test");
        ff.def_bondstyle("harmonic")
            .def_type("c3-c3", &[("k", 300.0), ("r0", 1.5)]);
        ff.def_pairstyle("lj/cut", &[])
            .def_type("c3", &[("epsilon", 0.05), ("sigma", 2.0)]);
        ff.def_pairstyle("coul/cut", &[]);

        // A 4-atom c3 chain with charges, bonds stretched to 1.6 (r0 = 1.5).
        let mut frame = Frame::new();
        let mut atoms = Block::new();
        atoms
            .insert(
                "x",
                Array1::from_vec(vec![0.0 as F, 1.6, 3.2, 4.8]).into_dyn(),
            )
            .unwrap();
        for col in ["y", "z"] {
            atoms
                .insert(col, Array1::from_vec(vec![0.0 as F; 4]).into_dyn())
                .unwrap();
        }
        atoms
            .insert(
                "type",
                Array1::from_vec(vec!["c3".to_string(); 4]).into_dyn(),
            )
            .unwrap();
        atoms
            .insert(
                "charge",
                Array1::from_vec(vec![0.1 as F, -0.1, 0.1, -0.1]).into_dyn(),
            )
            .unwrap();
        frame.insert("atoms", atoms);
        let mut bonds = Block::new();
        bonds
            .insert("atomi", Array1::from_vec(vec![0 as U, 1, 2]).into_dyn())
            .unwrap();
        bonds
            .insert("atomj", Array1::from_vec(vec![1 as U, 2, 3]).into_dyn())
            .unwrap();
        bonds
            .insert(
                "type",
                Array1::from_vec(vec!["c3-c3".to_string(); 3]).into_dyn(),
            )
            .unwrap();
        frame.insert("bonds", bonds);

        let coords: Vec<[F; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.6, 0.0, 0.0],
            [3.2, 0.0, 0.0],
            [4.8, 0.0, 0.0],
        ];
        let mut runner = LBFGSRelaxer::from_forcefield(ff).spawn(Some(&frame), &coords);
        let mut rng = rand::rng();
        // evaluate returns 0 → packing never worsens → accept the relaxation.
        let relaxed = runner
            .on_iter(&coords, 0.0, &mut |_| 0.0, &mut rng)
            .expect("a force field with non-bonded styles should compile and relax");
        let d = [
            relaxed[1][0] - relaxed[0][0],
            relaxed[1][1] - relaxed[0][1],
            relaxed[1][2] - relaxed[0][2],
        ];
        let b01 = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        assert!(
            b01 < 1.59,
            "stretched 1.6 bonds should relax toward r0 = 1.5, got {b01}"
        );
    }

    #[test]
    fn from_forcefield_without_frame_is_noop() {
        // No frame at spawn → the potential cannot be compiled → no-op runner
        // (the evaluate closure must never be reached).
        let mut ff = ForceField::new("test");
        ff.def_bondstyle("harmonic")
            .def_type("A-A", &[("k", 100.0), ("r0", 1.0)]);
        let coords = stretched();
        let mut runner = LBFGSRelaxer::from_forcefield(ff).spawn(None, &coords);
        let mut rng = rand::rng();
        let out = runner.on_iter(&coords, 0.0, &mut |_| panic!("must not evaluate"), &mut rng);
        assert!(out.is_none(), "no frame → the relaxer is a no-op");
    }
}
