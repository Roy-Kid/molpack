//! The growth driver: synchronous configurational-bias chain growth.
//!
//! By default every pending chain advances one step per **round**. Within a
//! round, every chain's candidate torsions are proposed and scored against the
//! **round-start snapshot** of the overlap field; commits then run serially in
//! the round's (shuffled) order, re-validating each selected candidate against
//! the atoms committed earlier in the same round. This snapshot semantics is
//! part of the algorithm's definition, not an implementation detail: a future
//! parallel version (parallel proposals + serial commits) is bit-identical
//! to this serial one by construction.
//!
//! [`GrowConfig::with_serial`] changes *which* chains a round touches, not what
//! a round does: only the first still-pending chain advances, so each chain is
//! finished before the next one starts. The round cap below still covers all
//! chains together, so a serial run needs `max_loops ≥ n_chains` to finish.
//!
//! Randomness is drawn from counter-based streams hashed per
//! `(seed, molecule, stage, visit)` — no global stream. Removing a chain
//! from the system therefore cannot change another chain's proposal
//! sequence, which is both a debugging property and the other half of the
//! parallel-equivalence guarantee.
//!
//! Dead ends (a proposal or commit that cannot place the pending atoms)
//! retract already-grown steps and may soften the hard core down to a floor.
//! Three clocks feed three readers: `deadend_streak` (cleared on a successful
//! commit) drives `retract_depth` and floor-level `force_due`;
//! `deadends_total` with `rungs_earned` is the cumulative watermark that
//! `rung_due` reads to earn a softening rung (one rung multiplies **this
//! chain's** dimensionless hard-core scale by [`GrowConfig::SOFTEN_RUNG`]).
//! Another chain that did not earn the rung keeps scale `1.0`.
//! [`StageOutcome::degraded`](crate::StageOutcome::degraded) counts
//! each such shrink *and* each forced
//! placement, and a structure is only `converged` when that counter is zero, so
//! the constructive no-overlap guarantee is asserted, never hoped for.
//!
//! The round loop is bounded. Growth reads `Budget::max_loops` as an
//! allowance of *passes over a chain*: growing a chain once costs
//! `n_steps + 1` rounds (the seed plus every step), so the cap is
//! `max(max_loops, 1) × (max n_steps + 1)` rounds — a zero budget still buys
//! one pass. Reaching it is a surrender, not a result: the pending chains are
//! force-completed exactly as on a handler abort, each forced placement counted
//! in `degraded`, and the outcome is `converged == false`.

use molrs::spatial::simbox::SimBox;
use molrs::types::F;

use crate::context::pack_state::evaluate_unscaled;
use crate::context::{PackState, Placed};
use crate::error::PackError;
use crate::grow::GrowError;
use crate::grow::config::GrowConfig;
use crate::grow::field::{BlockKind, OverlapField};
use crate::grow::internal::InternalTree;
use crate::grow::moves::{
    Chain, DeadEnd, Proposal, RestraintTable, SALT_SHUFFLE, Species, commit, force_place, propose,
    relax, retract, stream, uniform,
};
use crate::handler::{Handler, PhaseInfo, StageInfo, StepInfo};
use crate::stage::{Budget, Guarantees, Requires, Stage, StageOutcome};
use crate::target::Target;

/// The chain-growth stage. Built from the Grow targets before the first
/// [`run`](Stage::run); the targets handed to `run` must be the same objects.
///
/// The per-species trees are the stage's own configuration and stay on it
/// across runs, as the seam's re-entrancy contract requires.
///
/// **Rust-only:** not mirrored in the Python wheel — Python reaches growth
/// through the `CbmcGrow` entry, which constructs this stage itself.
pub struct GrowStage {
    seed: u64,
    species: Vec<Species>,
    /// The box this stage grows into, and the radius up-scaling its cell
    /// grid is sized from. `None` only for a stage built directly from
    /// templates and never handed a cell — it then grows in whatever box the
    /// state already carries; the entry always supplies one.
    cell: Option<(SimBox, F)>,
}

impl GrowStage {
    /// The name this stage reports, in one place: [`Stage::name`] returns it
    /// and `StepInfo.stage.name` is filled from it, so the two cannot drift
    /// apart.
    pub(crate) const NAME: &'static str = "growth";

    /// Build the per-species trees from the targets' templates. The `i`-th
    /// error names the offending target.
    pub fn from_targets(
        targets: &[Target],
        config: &GrowConfig,
        seed: u64,
    ) -> Result<Self, (usize, GrowError)> {
        let mut species = Vec::with_capacity(targets.len());
        for (i, t) in targets.iter().enumerate() {
            let tree = super::tree_from_target(t).map_err(|e| (i, e))?;
            species.push(Species {
                tree,
                prior: config.torsion_prior.clone(),
                cfg: config.clone(),
            });
        }
        Ok(Self {
            seed,
            species,
            cell: None,
        })
    }

    /// The compiled internal-coordinate tree of species `i` (the `i`-th
    /// target passed to [`from_targets`](Self::from_targets)).
    pub fn tree(&self, i: usize) -> &InternalTree {
        &self.species[i].tree
    }

    /// The box to grow into, with the `discale` its cell grid is sized from.
    ///
    /// Installed at the top of [`run`](Stage::run) rather than here: a stage
    /// is re-entrant, and the box is state the run owns, not configuration
    /// the stage consumes.
    pub(crate) fn with_resolved_cell(mut self, cell: SimBox, discale: F) -> Self {
        self.cell = Some((cell, discale));
        self
    }

    fn max_soft_shell(&self) -> F {
        self.species
            .iter()
            .map(|s| s.cfg.soft_shell)
            .fold(0.0, F::max)
    }
}

impl Stage for GrowStage {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    /// Nothing: growth constructs every placement itself, atom by atom.
    fn requires(&self) -> Requires {
        Requires::new(Placed::None)
    }

    /// Every free molecule placed.
    fn guarantees(&self) -> Guarantees {
        Guarantees::new(Placed::All)
    }

    fn run(
        &mut self,
        state: &mut PackState,
        _targets: &[Target],
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> Result<StageOutcome, PackError> {
        // ── The box and its cell grid ──────────────────────────────────────
        // See `install_resolved_cell` for why `radmax` reads `radius_ini`.
        if let Some((cell, discale)) = &self.cell {
            let sys = state.ctx_mut();
            crate::context::grid::install_resolved_cell(sys, cell, *discale);
        }

        let (sys, x) = state.rigid_split_mut();
        // ── Field over the final box ───────────────────────────────────────
        let origin: [F; 3] = {
            let v = sys.simbox.origin_view();
            [v[0], v[1], v[2]]
        };
        let lengths: [F; 3] = {
            let v = sys.simbox.lengths();
            [v[0], v[1], v[2]]
        };
        let pbc = sys.simbox.pbc();

        let mut mol_of = vec![0u32; sys.ntotat];
        let mut atom_of = vec![0u32; sys.ntotat];
        let mut chains: Vec<Chain> = Vec::with_capacity(sys.ntotmol);
        let mut mol = 0usize;
        for itype in 0..sys.ntype {
            let na = sys.natoms[itype];
            let n_steps = self.species[itype].tree.n_steps();
            for imol in 0..sys.nmols[itype] {
                let base = sys.idfirst[itype] + imol * na;
                for a in 0..na {
                    mol_of[base + a] = mol as u32;
                    atom_of[base + a] = a as u32;
                }
                chains.push(Chain {
                    itype,
                    mol,
                    base,
                    stage: 0,
                    coords: vec![[0.0; 3]; na],
                    vars: vec![0.0; self.species[itype].tree.n_vars()],
                    visits: vec![0; n_steps + 1],
                    deadends_total: 0,
                    deadend_streak: 0,
                    rungs_earned: 0,
                    hard_scale: 1.0,
                    relax_epoch: 0,
                });
                mol += 1;
            }
        }
        let n_chains = chains.len();

        let radmax = sys.radius.iter().cloned().fold(0.0, F::max);
        let cutoff = 2.0 * radmax + self.max_soft_shell();
        let mut field = OverlapField::new(
            origin,
            lengths,
            pbc,
            sys.radius.clone(),
            mol_of,
            atom_of,
            cutoff.max(1.0),
        );

        // ── Round loop ─────────────────────────────────────────────────────
        let mut degraded = 0usize;
        // Upper bound on the round loop: `max_loops` passes over the longest
        // chain, one round per stage (module docs). Without it a density the
        // hard core cannot satisfy spins forever instead of returning an
        // honest unconverged result (debt D-01 (ii)).
        let max_rounds = max_rounds(
            budget.max_loops,
            self.species
                .iter()
                .map(|s| s.tree.n_steps())
                .max()
                .unwrap_or(0),
        );
        let restraint_table = RestraintTable::from_context(sys);
        let mut aborted = false;
        let mut self_blocked = 0usize;
        let mut inter_chain = 0usize;
        let mut round: u64 = 0;
        loop {
            let done = chains
                .iter()
                .all(|c| c.stage > self.species[c.itype].tree.n_steps());
            if done {
                break;
            }
            if round >= max_rounds {
                // Out of rounds. Same contract as a handler abort: the
                // pending chains are completed below and nothing is claimed.
                aborted = true;
                break;
            }
            round += 1;

            // Chain order for this round: shuffled round-robin by default;
            // serial scheduling advances only the first pending chain, so
            // each chain completes into a finished matrix before the next
            // starts.
            // One config per run (the entry is the method), so the first
            // species speaks for all.
            let serial = self.species.first().is_some_and(|s| s.cfg.serial);
            let mut order: Vec<usize> = if serial {
                chains
                    .iter()
                    .position(|c| c.stage <= self.species[c.itype].tree.n_steps())
                    .into_iter()
                    .collect()
            } else {
                (0..n_chains).collect()
            };
            if !serial {
                let mut rng = stream(self.seed, u64::MAX, round, 0, SALT_SHUFFLE);
                for i in (1..order.len()).rev() {
                    let j = (uniform(&mut rng) * (i + 1) as F) as usize;
                    order.swap(i, j.min(i));
                }
            }

            // Proposal phase: every pending chain scores against the
            // round-start snapshot (no commits have happened yet).
            // Outer `None` = this chain was skipped this round.
            let mut proposals: Vec<Option<Result<Proposal, DeadEnd>>> =
                (0..n_chains).map(|_| None).collect();
            for &c in &order {
                let sp = &self.species[chains[c].itype];
                if chains[c].stage > sp.tree.n_steps() {
                    continue;
                }
                let p = propose(
                    &chains[c],
                    sp,
                    &field,
                    &restraint_table,
                    chains[c].hard_scale,
                    origin,
                    lengths,
                    self.seed,
                );
                proposals[c] = Some(p);
            }

            // Commit phase: serial, in round order, re-validating against
            // the live field (which accumulates this round's commits).
            for &c in &order {
                let Some(prop) = proposals[c].take() else {
                    continue;
                };
                let sp_idx = chains[c].itype;
                // Force only a chain that is truly wedged: the core is at its
                // floor AND the streak from *before* this round kept failing
                // through the escalating retractions. Anything less keeps
                // retrying — forced placements are the last resort that
                // breaks the constructive guarantee.
                let scale = chains[c].hard_scale;
                let min_scale = self.species[sp_idx].cfg.min_hard_scale;
                let force = force_due(
                    scale,
                    min_scale,
                    chains[c].deadend_streak,
                    self.species[sp_idx].cfg.soften_after,
                );
                let committed = match prop {
                    Ok(p) => {
                        match commit(&mut chains[c], &self.species[sp_idx], &mut field, p, scale) {
                            Ok(()) => true,
                            Err(kind) => {
                                match kind {
                                    BlockKind::SelfBlocked => self_blocked += 1,
                                    BlockKind::InterChain => inter_chain += 1,
                                }
                                false
                            }
                        }
                    }
                    Err(DeadEnd::Overlap(kind)) => {
                        match kind {
                            BlockKind::SelfBlocked => self_blocked += 1,
                            BlockKind::InterChain => inter_chain += 1,
                        }
                        false
                    }
                    Err(DeadEnd::Restraint) => false,
                };
                if committed {
                    chains[c].clear_streak();
                    continue;
                }
                // Every failed attempt — proposal dead end or commit-time
                // conflict — consumes its visit, so the next try at this
                // stage draws a fresh stream instead of replaying the same
                // rejected candidates (livelock otherwise).
                {
                    let stage = chains[c].stage;
                    let idx = stage.min(chains[c].visits.len() - 1);
                    chains[c].visits[idx] += 1;
                }
                // A chain that exhausted the escape ladder at the softening
                // floor is force-placed so it stops blocking the round loop —
                // each forced placement is counted as a softening event and
                // the result is not converged. Force does not increment clocks.
                if force {
                    force_place(
                        &mut chains[c],
                        &self.species[sp_idx],
                        &mut field,
                        origin,
                        lengths,
                        self.seed,
                    );
                    degraded += 1;
                    continue;
                }
                chains[c].record_dead_end();
                let cfg = &self.species[sp_idx].cfg;
                // Streak (cleared on commit) feeds retract; total never
                // resets, so a missed rung stays due. The ladder is the ONLY
                // path to a smaller hard core — a per-dead-end shrink is the
                // free-fall of debt D-01 (i).
                if rung_due(
                    chains[c].deadends_total,
                    chains[c].rungs_earned,
                    cfg.soften_after,
                ) && chains[c].hard_scale > cfg.min_hard_scale
                {
                    chains[c].hard_scale =
                        (chains[c].hard_scale * GrowConfig::SOFTEN_RUNG).max(cfg.min_hard_scale);
                    chains[c].rungs_earned += 1;
                    degraded += 1;
                }
                let depth = retract_depth(cfg.retract, chains[c].deadend_streak);
                retract(&mut chains[c], &self.species[sp_idx], &mut field, depth);
            }

            // Periodic tail regrowth against the now-denser field, guarded
            // so it never degrades (U_new must not exceed U_old). Cadence and
            // window are per species — each chain follows its own target's
            // config.
            for chain in chains.iter_mut() {
                let sp = &self.species[chain.itype];
                let every = sp.cfg.relax_every;
                if every > 0 && round.is_multiple_of(every as u64) {
                    relax(
                        chain,
                        sp,
                        &mut field,
                        &restraint_table,
                        chain.hard_scale,
                        sp.cfg.relax_window,
                        self.seed,
                    );
                }
            }

            // Handler visibility: one StepInfo per round. `radscale` carries
            // the current hard-core scale (softening is visible live);
            // fdist/frest are constructively 0 while growth holds its
            // guarantees, and `loop_idx` is the 1-based round number.
            for chain in &chains {
                for (a, p) in chain.coords.iter().enumerate() {
                    sys.xcart[chain.base + a] = *p;
                }
            }
            let info = StepInfo {
                // One stage per run until the pipeline lands; the name comes
                // from the stage type so the two cannot drift.
                stage: StageInfo {
                    index: 0,
                    total: 1,
                    name: Self::NAME,
                },
                loop_idx: round as usize,
                max_loops: budget.max_loops,
                phase: PhaseInfo {
                    phase: 0,
                    total_phases: 1,
                    molecule_type: None,
                },
                fdist: 0.0,
                frest: 0.0,
                f: 0.0,
                improvement_pct: 0.0,
                // Softest core this round. One chain's rung no longer moves
                // the others, so the live report is the minimum.
                radscale: chains.iter().map(|c| c.hard_scale).fold(1.0, F::min),
                precision: budget.precision,
            };
            for h in handlers.iter_mut() {
                h.on_step(&info, sys);
            }
            if handlers.iter().any(|h| h.should_stop()) {
                aborted = true;
                break;
            }
        }

        // An abort — a handler stop or the exhausted round cap — leaves later
        // stages at the origin sentinel in `Chain.coords`. Completing them
        // with force_place (hard core ignored) keeps bonded geometry chemical
        // so assemble_frame does not emit 0-length or box-scale bonds. Each
        // forced placement breaks the constructive guarantee exactly like a
        // softening rung and is counted as one; the result stays unconverged.
        if aborted {
            let denom = self_blocked + inter_chain;
            if denom > 0 {
                log::warn!(
                    "growth aborted: {self_blocked}/{denom} overlap dead ends were self-blocked. \
                     If intramolecular contacts dominate, shrink radii with Target::with_atom_radius \
                     or widen intramolecular exclusions / special bonds"
                );
            } else {
                log::warn!(
                    "growth aborted. If intramolecular contacts dominate, shrink radii with \
                     Target::with_atom_radius or widen intramolecular exclusions / special bonds"
                );
            }
            for chain in &mut chains {
                let sp = &self.species[chain.itype];
                while chain.stage <= sp.tree.n_steps() {
                    force_place(chain, sp, &mut field, origin, lengths, self.seed);
                    degraded += 1;
                }
                // `sys.xcart` is the single lab-frame home of the placed
                // atoms; the round loop syncs it at every round end, so the
                // forced completion has to sync it as well or the stages
                // placed here would live only in `Chain.coords`. Same shape
                // as the lattice driver's abort path (`grow/lattice/mod.rs`).
                for (a, p) in chain.coords.iter().enumerate() {
                    sys.xcart[chain.base + a] = *p;
                }
            }
        }

        // ── Writeback: `RigidView::capture_from_xcart` ─────────────────────
        // `sys.xcart` is the lab-frame home every chain has synced into — at
        // each round end and, for the forced completion above, in the abort
        // loop — so the view captures the placements from there. One
        // derivation, shared with the lattice driver.
        x.capture_from_xcart(sys);

        // ── Final verdict from the shared objective, never self-reported ──
        // One unscaled evaluation, the crate's single primitive for it: it
        // sets the unscaled `scale` / `scale2` pair this site used to write
        // inline and gives the caller's values back afterwards.
        let (_, fdist, frest) = evaluate_unscaled(sys, x.as_slice());
        let converged = !aborted && degraded == 0 && fdist == 0.0 && frest < budget.precision;
        Ok(StageOutcome::new(converged, degraded))
    }
}

/// Upper bound on the round loop.
///
/// Growth reads `max_loops` as an allowance of *passes over a chain*, and one
/// pass costs `n_steps + 1` rounds (the seed plus every step), so the cap is
/// `max(max_loops, 1) × (max_n_steps + 1)`: a zero budget still buys one pass.
/// Saturating, because the product of two user-supplied counts must not wrap
/// into a small cap — that would turn a generous budget into an early
/// surrender. Without any cap a density the hard core cannot satisfy spins
/// forever instead of returning an honest unconverged result (debt D-01 (ii)).
fn max_rounds(max_loops: usize, max_n_steps: usize) -> u64 {
    (max_loops.max(1) as u64).saturating_mul((max_n_steps as u64).saturating_add(1))
}

/// Retract depth from the consecutive dead-end streak.
///
/// `base.saturating_mul(1 << (streak / 4).min(12))`. Reads `deadend_streak`
/// only — never `deadends_total`. Feeding the cumulative clock would restore
/// the exponential whole-chain retract of debt D-01.
fn retract_depth(base: usize, streak: usize) -> usize {
    base.saturating_mul(1 << (streak / 4).min(12))
}

/// Whether this chain has earned another softening rung.
///
/// `total >= (rungs_earned + 1).saturating_mul(soften_after)`. Reads
/// `deadends_total` against the `rungs_earned` watermark; never
/// `deadend_streak`. A missed take stays due, and a successful commit does
/// not reset the cumulative counter.
fn rung_due(total: usize, rungs_earned: usize, soften_after: usize) -> bool {
    total >= (rungs_earned + 1).saturating_mul(soften_after)
}

/// Whether this chain is wedged at the softening floor and must force-place.
///
/// `hard_scale <= min && streak >= 2 * soften_after` (saturating). Reads
/// `deadend_streak` at the floor; never `deadends_total` or `rungs_earned`.
/// `hard_scale` is dimensionless (`1.0` = full declared contact).
fn force_due(hard_scale: F, min: F, streak: usize, soften_after: usize) -> bool {
    hard_scale <= min && streak >= 2usize.saturating_mul(soften_after)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cap is an allowance of passes over a chain, not of rounds: a zero
    /// budget still buys one pass, and the product never wraps.
    #[test]
    fn max_rounds_is_passes_times_chain_length() {
        assert_eq!(max_rounds(1, 4), 5, "one pass = seed + 4 steps");
        assert_eq!(max_rounds(3, 4), 15);
        assert_eq!(max_rounds(0, 4), 5, "a zero budget still buys one pass");
        assert_eq!(
            max_rounds(2, 0),
            2,
            "a rigid template is one round per pass"
        );
        assert_eq!(
            max_rounds(usize::MAX, usize::MAX),
            u64::MAX,
            "the cap saturates — wrapping would turn a generous budget into an \
             early surrender"
        );
    }

    #[test]
    fn retract_depth_reads_deadend_streak() {
        assert_eq!(retract_depth(10, 0), 10);
        assert_eq!(retract_depth(10, 3), 10);
        assert_eq!(retract_depth(10, 4), 20);
        assert_eq!(retract_depth(10, 8), 40);
        // Shift cap: `(streak / 4).min(12)` — 48 consecutive is 12 bins.
        assert_eq!(retract_depth(10, 48), 10 * (1usize << 12));
    }

    #[test]
    fn force_due_reads_streak_at_floor() {
        // Spec edge: streak 0 at the floor must not force-place.
        assert!(!force_due(0.8 as F, 0.8 as F, 0, 50));
        assert!(!force_due(0.8 as F, 0.8 as F, 99, 50));
        assert!(force_due(0.8 as F, 0.8 as F, 100, 50));
        assert!(!force_due(0.81 as F, 0.8 as F, 100, 50));
    }

    #[test]
    fn rung_due_uses_watermark() {
        assert!(!rung_due(0, 0, 50));
        assert!(!rung_due(49, 0, 50));
        assert!(rung_due(50, 0, 50));
        // First rung already taken; next is due at 100, not at 50 again.
        assert!(!rung_due(50, 1, 50));
        assert!(rung_due(100, 1, 50));
        assert!(rung_due(usize::MAX, 0, 2));
    }

    #[test]
    fn missed_rung_still_due() {
        // The old `is_multiple_of` gate skipped a missed take (3 % 2 != 0)
        // until the next multiple. A watermark stays due once crossed.
        assert!(rung_due(3, 0, 2));
        assert!(rung_due(2, 0, 2));
        assert!(!rung_due(3, 1, 2));
    }
}
