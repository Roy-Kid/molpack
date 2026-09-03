//! Growth-move machinery: proposal, scoring, commit, retraction, forced
//! placement, and the W-guarded tail regrowth — plus the per-chain state
//! ([`Chain`], [`Species`]) and the hashed RNG streams they draw from.
//!
//! The [`driver`](super::driver) owns the round loop and the `Solver`
//! contract; this module owns everything a single chain does within a round.
//! The split follows the file-size budget, not a semantic boundary shift:
//! the round-snapshot semantics documented on the parent module bind both.

use std::sync::Arc;

use molrs::types::F;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::context::PackContext;
use crate::euler::eulerrmat;
use crate::grow::config::GrowConfig;
use crate::grow::field::{OverlapField, Probe};
use crate::grow::internal::InternalTree;
use crate::grow::prior::{AnglePrior, TorsionPrior};
use crate::random::uniform01;
use crate::restraint::AtomRestraint;

pub(super) const TWO_PI: F = std::f64::consts::TAU as F;

/// Stream salts: keep the (mol, stage, visit) keying collision-free across
/// the driver's independent random purposes.
pub(super) const SALT_PROPOSE: u64 = 0x6772_6f77_0000_0001;
pub(super) const SALT_RELAX: u64 = 0x6772_6f77_0000_0002;
pub(super) const SALT_SHUFFLE: u64 = 0x6772_6f77_0000_0003;

/// One growable species: its internal-coordinate tree and resolved config.
pub(super) struct Species {
    pub(super) tree: InternalTree,
    pub(super) prior: TorsionPrior,
    pub(super) cfg: GrowConfig,
}

/// Mutable growth state of one molecule copy.
pub(super) struct Chain {
    pub(super) itype: usize,
    /// Global molecule index (type-major, copy-major — the `x` order).
    pub(super) mol: usize,
    /// First `icart` of this copy.
    pub(super) base: usize,
    /// 0 = seed pending; `1 + k` = tree step `k` pending; done at
    /// `1 + n_steps`.
    pub(super) stage: usize,
    /// Lab-frame coordinates, valid for placed atoms only (unwrapped).
    pub(super) coords: Vec<[F; 3]>,
    /// Free-torsion values chosen so far.
    pub(super) vars: Vec<F>,
    /// Per-stage visit counters — retraction revisits draw fresh streams.
    pub(super) visits: Vec<u32>,
    /// Consecutive dead ends at the current stage.
    pub(super) deadends: usize,
    /// Relax epochs completed (keys the relax streams).
    pub(super) relax_epoch: u32,
}

/// `(icart, template_atom, lab position)` of one placed atom.
pub(super) type PlacedAtom = (usize, usize, [F; 3]);

/// A selected candidate for one chain's pending stage, plus the surviving
/// alternatives for commit-time conflict resolution.
pub(super) struct Proposal {
    pub(super) atoms: Vec<PlacedAtom>,
    pub(super) var: Option<F>,
    /// Alternatives, best-first (excluding the selected one).
    pub(super) alternatives: Vec<Trial>,
}

pub(super) struct Trial {
    pub(super) atoms: Vec<PlacedAtom>,
    pub(super) var: Option<F>,
    pub(super) penalty: F,
}

/// Per-atom restraint lookup, cloned out of the context once so the round
/// loop holds no borrow on `sys`. Restraints are **hard** during growth: a
/// candidate violating any of its atom's restraints is rejected outright,
/// exactly like a hard-core overlap — which is what makes `frest == 0` a
/// constructive guarantee rather than a convergence hope (spec Design §3).
pub(super) struct RestraintTable {
    pub(super) offsets: Vec<usize>,
    pub(super) data: Vec<usize>,
    pub(super) restraints: Vec<Arc<dyn AtomRestraint>>,
}

impl RestraintTable {
    pub(super) fn from_context(sys: &PackContext) -> Self {
        Self {
            offsets: sys.iratom_offsets.clone(),
            data: sys.iratom_data.clone(),
            restraints: sys.restraints.clone(),
        }
    }

    /// `true` when any restraint on atom `icart` is violated at `p`.
    /// Scales mirror the shared objective's final-evaluation settings.
    pub(super) fn violated(&self, icart: usize, p: &[F; 3]) -> bool {
        self.data[self.offsets[icart]..self.offsets[icart + 1]]
            .iter()
            .any(|&r| self.restraints[r].f(p, 1.0, 0.01) > 0.0)
    }
}

/// splitmix64-style hash of the stream key into an RNG seed.
pub(super) fn stream(seed: u64, mol: u64, stage: u64, visit: u64, salt: u64) -> SmallRng {
    let mut z = seed
        ^ mol.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ stage.wrapping_mul(0xBF58_476D_1CE4_E5B9)
        ^ visit.wrapping_mul(0x94D0_49BB_1331_11EB)
        ^ salt.wrapping_mul(0xD6E8_FEB8_6659_FD93);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    SmallRng::seed_from_u64(z ^ (z >> 31))
}

pub(super) fn uniform(rng: &mut SmallRng) -> F {
    uniform01(rng)
}

/// Seed-anchor draw. Uniform over the box by default; with
/// [`GrowConfig::with_void_bias`] the anchor comes from a uniformly chosen
/// *empty* field cell (cavity seeding), falling back to uniform when no cell
/// is empty. Consumes exactly 3 uniforms per call (4 with void bias), so the
/// stream layout is fixed within a mode.
fn draw_anchor(
    cfg: &GrowConfig,
    field: &OverlapField,
    origin: [F; 3],
    lengths: [F; 3],
    rng: &mut SmallRng,
) -> [F; 3] {
    if cfg.void_bias {
        let u = [uniform(rng), uniform(rng), uniform(rng), uniform(rng)];
        return field.empty_cell_point(u).unwrap_or([
            origin[0] + u[1] * lengths[0],
            origin[1] + u[2] * lengths[1],
            origin[2] + u[3] * lengths[2],
        ]);
    }
    [
        origin[0] + uniform(rng) * lengths[0],
        origin[1] + uniform(rng) * lengths[1],
        origin[2] + uniform(rng) * lengths[2],
    ]
}

/// Place step `k` into `scratch`, sampling placement angles when the
/// species' angle prior asks for it ([`AnglePrior::Template`] consumes no
/// randomness, so the all-atom stream layout is unchanged).
pub(super) fn place_step_sampled(
    sp: &Species,
    k: usize,
    vars: &[F],
    scratch: &mut [[F; 3]],
    rng: &mut SmallRng,
) {
    match &sp.cfg.angle_prior {
        AnglePrior::Template => sp.tree.place_step(k, vars, scratch),
        prior => {
            let angles: Vec<F> = sp
                .tree
                .step_angles(k)
                .map(|t| prior.sample_interior(t, rng))
                .collect();
            sp.tree.place_step_with_angles(k, vars, &angles, scratch);
        }
    }
}

/// Propose one stage for `chain` against the snapshot field. `None` = dead
/// end (every trial hard-blocked or restraint-rejected).
#[allow(clippy::too_many_arguments)]
pub(super) fn propose(
    chain: &Chain,
    sp: &Species,
    field: &OverlapField,
    restraints: &RestraintTable,
    hard_scale: F,
    origin: [F; 3],
    lengths: [F; 3],
    seed: u64,
) -> Option<Proposal> {
    let visit = chain.visits[chain.stage.min(chain.visits.len() - 1)];
    let mut rng = stream(
        seed,
        chain.mol as u64,
        chain.stage as u64,
        visit as u64,
        SALT_PROPOSE,
    );
    let cfg = &sp.cfg;
    let cap = 60.0 / cfg.selectivity.max(0.1);

    let mut trials: Vec<Trial> = Vec::with_capacity(cfg.trials);
    if chain.stage == 0 {
        // Seed placement: anchor position + orientation trials.
        for _ in 0..cfg.trials {
            let anchor = draw_anchor(cfg, field, origin, lengths, &mut rng);
            let (v1, v2, v3) = eulerrmat(
                uniform(&mut rng) * TWO_PI,
                uniform(&mut rng) * TWO_PI,
                uniform(&mut rng) * TWO_PI,
            );
            let rot = [
                [v1[0], v2[0], v3[0]],
                [v1[1], v2[1], v3[1]],
                [v1[2], v2[2], v3[2]],
            ];
            let mut scratch = chain.coords.clone();
            sp.tree.place_seed(anchor, &rot, &mut scratch);
            if let Some((atoms, penalty)) = score_atoms(
                sp.tree.seed_atoms().iter().copied(),
                &scratch,
                chain,
                sp,
                field,
                restraints,
                hard_scale,
                cap,
            ) {
                trials.push(Trial {
                    atoms,
                    var: None,
                    penalty,
                });
            }
        }
    } else {
        let k = chain.stage - 1;
        let n_trials = if sp.tree.step_var(k).is_some() {
            cfg.trials
        } else {
            1
        };
        for _ in 0..n_trials {
            let mut vars = chain.vars.clone();
            let var = sp.tree.step_var(k).map(|v| {
                let phi = sp.prior.sample(sp.tree.template_var(k), &mut rng);
                vars[v] = phi;
                phi
            });
            let mut scratch = chain.coords.clone();
            place_step_sampled(sp, k, &vars, &mut scratch, &mut rng);
            if let Some((atoms, penalty)) = score_atoms(
                sp.tree.step_atoms(k),
                &scratch,
                chain,
                sp,
                field,
                restraints,
                hard_scale,
                cap,
            ) {
                trials.push(Trial {
                    atoms,
                    var,
                    penalty,
                });
            }
        }
    }
    if trials.is_empty() {
        return None;
    }

    // Rosenbluth selection among the surviving trials (log-w safe: weights
    // are exp(-β(U - U_min))).
    let beta = cfg.selectivity;
    let u_min = trials.iter().map(|t| t.penalty).fold(F::INFINITY, F::min);
    let weights: Vec<F> = trials
        .iter()
        .map(|t| (-beta * (t.penalty - u_min)).exp())
        .collect();
    let total: F = weights.iter().sum();
    let mut ticket = uniform(&mut rng) * total;
    let mut chosen = trials.len() - 1;
    for (i, w) in weights.iter().enumerate() {
        if ticket < *w {
            chosen = i;
            break;
        }
        ticket -= *w;
    }
    let selected = trials.swap_remove(chosen);
    trials.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));
    Some(Proposal {
        atoms: selected.atoms,
        var: selected.var,
        alternatives: trials,
    })
}

/// Score one trial's atoms against the field and the restraints.
/// `Some((atoms, ΣU))`, or `None` when any atom is hard-blocked or violates
/// a restraint (restraints are hard during growth — spec Design §3).
#[allow(clippy::too_many_arguments)]
pub(super) fn score_atoms(
    atoms: impl Iterator<Item = usize>,
    scratch: &[[F; 3]],
    chain: &Chain,
    sp: &Species,
    field: &OverlapField,
    restraints: &RestraintTable,
    hard_scale: F,
    cap: F,
) -> Option<(Vec<PlacedAtom>, F)> {
    let mut placed = Vec::new();
    let mut penalty = 0.0;
    for a in atoms {
        let slot = chain.base + a;
        if restraints.violated(slot, &scratch[a]) {
            return None;
        }
        match field.probe(
            slot,
            scratch[a],
            sp.tree.exclusions(a),
            hard_scale,
            sp.cfg.soft_shell,
            cap,
        ) {
            Probe::Blocked => return None,
            Probe::Room(p) => penalty += p,
        }
        placed.push((slot, a, scratch[a]));
    }
    Some((placed, penalty))
}

/// Commit a proposal against the live field: probe-then-insert per atom so
/// this round's earlier commits (and the step's own atoms) are seen. Falls
/// back to the proposal's alternatives; `false` = nothing survived.
pub(super) fn commit(
    chain: &mut Chain,
    sp: &Species,
    field: &mut OverlapField,
    prop: Proposal,
    hard_scale: F,
) -> bool {
    let cap = 60.0 / sp.cfg.selectivity.max(0.1);
    let mut candidates = Vec::with_capacity(1 + prop.alternatives.len());
    candidates.push((prop.atoms, prop.var));
    for t in prop.alternatives {
        candidates.push((t.atoms, t.var));
    }
    'cand: for (atoms, var) in candidates {
        let mut inserted: Vec<usize> = Vec::with_capacity(atoms.len());
        for &(slot, a, p) in &atoms {
            match field.probe(slot, p, sp.tree.exclusions(a), hard_scale, 0.0, cap) {
                Probe::Blocked => {
                    for &s in &inserted {
                        field.remove(s);
                    }
                    continue 'cand;
                }
                Probe::Room(_) => {
                    field.insert(slot, p);
                    inserted.push(slot);
                }
            }
        }
        // Accepted: record coordinates and advance.
        for &(_, a, p) in &atoms {
            chain.coords[a] = p;
        }
        if let (Some(phi), Some(v)) = (var, step_var_of(chain, sp)) {
            chain.vars[v] = phi;
        }
        let idx = chain.stage.min(chain.visits.len() - 1);
        chain.visits[idx] += 1;
        chain.stage += 1;
        return true;
    }
    // Nothing survived commit-time validation. The caller counts the visit
    // (uniformly for commit failures and proposal dead ends) so the next
    // proposal draws a fresh stream.
    false
}

pub(super) fn step_var_of(chain: &Chain, sp: &Species) -> Option<usize> {
    if chain.stage == 0 {
        None
    } else {
        sp.tree.step_var(chain.stage - 1)
    }
}

/// Remove the last `depth` stages' atoms and step back. Retraction may reach
/// stage 0 (the seed is re-placed from a fresh stream) — a variable-less
/// molecule has no other escape.
pub(super) fn retract(chain: &mut Chain, sp: &Species, field: &mut OverlapField, depth: usize) {
    let target_stage = chain.stage.saturating_sub(depth);
    while chain.stage > target_stage {
        chain.stage -= 1;
        if chain.stage == 0 {
            for a in sp.tree.seed_atoms() {
                field.remove(chain.base + a);
            }
        } else {
            for a in sp.tree.step_atoms(chain.stage - 1) {
                field.remove(chain.base + a);
            }
        }
    }
}

/// Last resort: place the pending stage at the least-bad candidate, ignoring
/// the hard core, so the driver can terminate.
///
/// The driver calls this in exactly two situations — a chain still wedged at
/// the softening floor, and an abort (a handler stop, or the exhausted round
/// cap) that leaves later stages unplaced. Either way the constructive
/// no-overlap guarantee is broken, so the caller counts every call as a
/// softening event and the run's result is never `converged`.
pub(super) fn force_place(
    chain: &mut Chain,
    sp: &Species,
    field: &mut OverlapField,
    origin: [F; 3],
    lengths: [F; 3],
    seed: u64,
) {
    let visit = chain.visits[chain.stage.min(chain.visits.len() - 1)];
    let mut rng = stream(
        seed,
        chain.mol as u64,
        chain.stage as u64,
        visit as u64,
        SALT_PROPOSE,
    );
    // Least-bad among the usual number of trials: rank by the smallest
    // non-excluded neighbour distance a candidate's worst atom would have,
    // and take the roomiest. Never a blind single draw — a forced placement
    // is already a counted failure, and a near-coincident pair would poison
    // the downstream relaxation.
    let mut best: Option<(F, Vec<PlacedAtom>, Option<F>)> = None;
    for _ in 0..sp.cfg.trials.max(1) {
        let (atoms, var): (Vec<PlacedAtom>, Option<F>) = if chain.stage == 0 {
            let anchor = draw_anchor(&sp.cfg, field, origin, lengths, &mut rng);
            let (v1, v2, v3) = eulerrmat(
                uniform(&mut rng) * TWO_PI,
                uniform(&mut rng) * TWO_PI,
                uniform(&mut rng) * TWO_PI,
            );
            let rot = [
                [v1[0], v2[0], v3[0]],
                [v1[1], v2[1], v3[1]],
                [v1[2], v2[2], v3[2]],
            ];
            let mut scratch = chain.coords.clone();
            sp.tree.place_seed(anchor, &rot, &mut scratch);
            (
                sp.tree
                    .seed_atoms()
                    .iter()
                    .map(|&a| (chain.base + a, a, scratch[a]))
                    .collect(),
                None,
            )
        } else {
            let k = chain.stage - 1;
            let mut vars = chain.vars.clone();
            let var = sp.tree.step_var(k).map(|v| {
                let phi = sp.prior.sample(sp.tree.template_var(k), &mut rng);
                vars[v] = phi;
                phi
            });
            let mut scratch = chain.coords.clone();
            place_step_sampled(sp, k, &vars, &mut scratch, &mut rng);
            (
                sp.tree
                    .step_atoms(k)
                    .map(|a| (chain.base + a, a, scratch[a]))
                    .collect(),
                var,
            )
        };
        let worst = atoms
            .iter()
            .map(|&(slot, a, p)| field.nearest(slot, p, sp.tree.exclusions(a)))
            .fold(F::INFINITY, F::min);
        if best.as_ref().is_none_or(|(b, _, _)| worst > *b) {
            best = Some((worst, atoms, var));
        }
    }
    let (_, atoms, var) = best.expect("at least one forced trial");
    for &(slot, a, p) in &atoms {
        field.insert(slot, p);
        chain.coords[a] = p;
    }
    if let (Some(phi), Some(v)) = (var, step_var_of(chain, sp)) {
        chain.vars[v] = phi;
    }
    let idx = chain.stage.min(chain.visits.len() - 1);
    chain.visits[idx] += 1;
    chain.stage += 1;
}

/// W-guarded tail regrowth: retract `window` steps and regrow them against
/// the current field; keep the new tail only when its total crowding penalty
/// does not exceed the old tail's (re-scored the same way). Restores the old
/// tail verbatim otherwise.
#[allow(clippy::too_many_arguments)]
pub(super) fn relax(
    chain: &mut Chain,
    sp: &Species,
    field: &mut OverlapField,
    restraints: &RestraintTable,
    hard_scale: F,
    window: usize,
    seed: u64,
) {
    let n_steps = sp.tree.n_steps();
    if chain.stage <= window + 1 || chain.stage <= n_steps {
        // Only relax completed chains (simplest sound policy: partial chains
        // are still growing and will be revisited anyway).
        return;
    }
    let cap = F::INFINITY;
    let first = n_steps - window;

    // Save the old tail and remove it from the field.
    let old_coords = chain.coords.clone();
    let old_vars = chain.vars.clone();
    for k in first..n_steps {
        for a in sp.tree.step_atoms(k) {
            field.remove(chain.base + a);
        }
    }

    // Score the old tail exactly as a regrowth would see it.
    let mut u_old = 0.0;
    let mut old_ok = true;
    'old: for k in first..n_steps {
        for a in sp.tree.step_atoms(k) {
            let slot = chain.base + a;
            match field.probe(
                slot,
                old_coords[a],
                sp.tree.exclusions(a),
                hard_scale,
                sp.cfg.soft_shell,
                cap,
            ) {
                Probe::Blocked => {
                    old_ok = false;
                    break 'old;
                }
                Probe::Room(p) => u_old += p,
            }
            field.insert(slot, old_coords[a]);
        }
    }
    // Remove whatever the scoring pass re-inserted.
    for k in first..n_steps {
        for a in sp.tree.step_atoms(k) {
            field.remove(chain.base + a);
        }
    }
    if !old_ok {
        u_old = F::INFINITY;
    }

    // Regrow the tail greedily against the live field.
    chain.relax_epoch += 1;
    let mut rng = stream(
        seed,
        chain.mol as u64,
        chain.relax_epoch as u64,
        0,
        SALT_RELAX,
    );
    let mut u_new = 0.0;
    let mut new_ok = true;
    'grow: for k in first..n_steps {
        type RelaxBest = (F, Vec<(usize, [F; 3])>, Option<F>);
        let mut best: Option<RelaxBest> = None;
        let n_trials = if sp.tree.step_var(k).is_some() {
            sp.cfg.trials
        } else {
            1
        };
        for _ in 0..n_trials {
            let mut vars = chain.vars.clone();
            let var = sp.tree.step_var(k).map(|v| {
                let phi = sp.prior.sample(sp.tree.template_var(k), &mut rng);
                vars[v] = phi;
                phi
            });
            let mut scratch = chain.coords.clone();
            place_step_sampled(sp, k, &vars, &mut scratch, &mut rng);
            let mut u = 0.0;
            let mut ok = true;
            let mut atoms = Vec::new();
            for a in sp.tree.step_atoms(k) {
                if restraints.violated(chain.base + a, &scratch[a]) {
                    ok = false;
                    break;
                }
                match field.probe(
                    chain.base + a,
                    scratch[a],
                    sp.tree.exclusions(a),
                    hard_scale,
                    sp.cfg.soft_shell,
                    cap,
                ) {
                    Probe::Blocked => {
                        ok = false;
                        break;
                    }
                    Probe::Room(p) => u += p,
                }
                atoms.push((a, scratch[a]));
            }
            if ok && best.as_ref().is_none_or(|(bu, _, _)| u < *bu) {
                best = Some((u, atoms, var));
            }
        }
        match best {
            Some((u, atoms, var)) => {
                for &(a, p) in &atoms {
                    chain.coords[a] = p;
                    field.insert(chain.base + a, p);
                }
                if let (Some(phi), Some(v)) = (var, sp.tree.step_var(k)) {
                    chain.vars[v] = phi;
                }
                u_new += u;
            }
            None => {
                new_ok = false;
                break 'grow;
            }
        }
    }

    if !new_ok || u_new > u_old {
        // Restore the old tail verbatim.
        for k in first..n_steps {
            for a in sp.tree.step_atoms(k) {
                let slot = chain.base + a;
                if field.is_placed(slot) {
                    field.remove(slot);
                }
            }
        }
        chain.coords = old_coords;
        chain.vars = old_vars;
        for k in first..n_steps {
            for a in sp.tree.step_atoms(k) {
                field.insert(chain.base + a, chain.coords[a]);
            }
        }
    }
}
