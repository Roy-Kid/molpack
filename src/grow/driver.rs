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
//! Dead ends retract (recoil at feeler depth 1 ≡ CBMC); repeated dead ends
//! soften the hard core down to a floor. Softening is a **ladder**: a rung is
//! earned when a single chain reaches another `soften_after` dead ends, and the
//! driver takes at most one rung per round even when two chains earn one in the
//! same round, because the rung is the driver's global step, not a per-chain
//! one. [`StageOutcome::softened`](crate::StageOutcome::softened) counts
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
//! in `softened`, and the outcome is `converged == false`.

use molrs::spatial::simbox::SimBox;
use molrs::types::F;

use crate::context::pack_state::evaluate_unscaled;
use crate::context::{PackState, Placed};
use crate::error::PackError;
use crate::grow::config::{GrowConfig, GrowError};
use crate::grow::field::OverlapField;
use crate::grow::internal::InternalTree;
use crate::grow::moves::{
    Chain, Proposal, RestraintTable, SALT_SHUFFLE, Species, commit, force_place, propose, relax,
    retract, stream, uniform,
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
            let frame = t.template.as_ref().ok_or((i, GrowError::MissingTemplate))?;
            let tree = InternalTree::from_frame_with_depth(frame, config.exclusion_depth)
                .map_err(|e| (i, e))?;
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
            crate::initial::install_resolved_cell(sys, cell, *discale);
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
                    deadends: 0,
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
        let mut hard_scale: F = 1.0;
        let mut softened = 0usize;
        // Upper bound on the round loop: `max_loops` passes over the longest
        // chain, one round per stage (module docs). Without it a density the
        // hard core cannot satisfy spins forever instead of returning an
        // honest unconverged result (debt D-01 (ii)).
        let max_rounds = {
            let max_n_steps = self
                .species
                .iter()
                .map(|s| s.tree.n_steps())
                .max()
                .unwrap_or(0);
            (budget.max_loops.max(1) as u64).saturating_mul(max_n_steps as u64 + 1)
        };
        let min_hard_scale = self
            .species
            .iter()
            .map(|s| s.cfg.min_hard_scale)
            .fold(1.0, F::min);

        let restraint_table = RestraintTable::from_context(sys);
        let mut aborted = false;
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
            // At most one softening rung per round (see the ladder in the
            // module docs). Two chains can legitimately earn a rung in the
            // same round; the core still shrinks once, because the rung is
            // the driver's global step, not a per-chain one.
            let mut rung_this_round = false;

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
            let mut proposals: Vec<Option<Option<Proposal>>> =
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
                    hard_scale,
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
                // floor AND the chain kept failing through the escalating
                // retractions (including full restarts). Anything less keeps
                // retrying — forced placements are the last resort that
                // breaks the constructive guarantee.
                // `2 * soften_after` = two full ladder rungs of dead ends past
                // the floor; one rung's worth is normal softening pressure, not
                // a wedged chain.
                let force = hard_scale <= min_hard_scale
                    && chains[c].deadends >= 2 * self.species[sp_idx].cfg.soften_after;
                let committed = match prop {
                    Some(p) => commit(
                        &mut chains[c],
                        &self.species[sp_idx],
                        &mut field,
                        p,
                        hard_scale,
                    ),
                    None => false,
                };
                if committed {
                    chains[c].deadends = 0;
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
                // Dead end. A chain that exhausted the escape ladder at the
                // softening floor is force-placed so it stops blocking the
                // round loop — each forced placement is counted as a
                // softening event and the result is not converged.
                if force {
                    force_place(
                        &mut chains[c],
                        &self.species[sp_idx],
                        &mut field,
                        origin,
                        lengths,
                        self.seed,
                    );
                    softened += 1;
                    continue;
                }
                chains[c].deadends += 1;
                let cfg = &self.species[sp_idx].cfg;
                // The failure count resets only on a successful commit —
                // softening must NOT reset it, or the retraction depth
                // never escalates and a chain whose growth front is wedged
                // in a dense pocket replays shallow retractions forever
                // (a bonded atom is not a free insertion: it must sit on
                // its parent's bond sphere, so escaping a pocket needs the
                // deep retractions).
                // The ladder is the ONLY path to a smaller hard core: a rung
                // costs `soften_after` dead ends on THIS chain. An exhausted
                // budget of any kind may not shrink the core per dead end —
                // that turns the ladder into a free-fall (1.0 → the floor
                // inside one round, debt D-01 (i)); a run that cannot make
                // progress terminates through the round cap instead.
                if chains[c].deadends.is_multiple_of(cfg.soften_after)
                    && !rung_this_round
                    && hard_scale > min_hard_scale
                {
                    hard_scale = (hard_scale * 0.97).max(min_hard_scale);
                    rung_this_round = true;
                    softened += 1;
                }
                let depth = cfg
                    .retract
                    .saturating_mul(1 << (chains[c].deadends / 4).min(12));
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
                        hard_scale,
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
                improvement_pct: 0.0,
                radscale: hard_scale,
                precision: budget.precision,
                relaxer_acceptance: Vec::new(),
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
            for chain in &mut chains {
                let sp = &self.species[chain.itype];
                while chain.stage <= sp.tree.n_steps() {
                    force_place(chain, sp, &mut field, origin, lengths, self.seed);
                    softened += 1;
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
        let converged = !aborted && softened == 0 && fdist == 0.0 && frest < budget.precision;
        Ok(StageOutcome::new(converged, softened))
    }
}
