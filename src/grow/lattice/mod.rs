//! Lattice growth: diamond-lattice SAW pre-generation for dense melts
//! (lattice-growth-phase spec; tetrahedral trees in lattice-branch-saw).
//!
//! The continuum growth solver hits a density ceiling long before melt
//! density: its torsion candidates are proposed in continuous space and
//! ground down by retraction. On the diamond lattice the same RIS states
//! are EXACT lattice moves — trans/gauche± are the three non-backtracking
//! continuations — and excluded volume is an O(1) site-occupancy check, so
//! a melt-density walk completes in milliseconds where continuum growth
//! grinds.
//!
//! **A template supplies its topology, not its geometry.** The walk chooses
//! the route and `decorate` seats each backbone atom on the site the walk
//! chose, so the route's guarantees — self-avoidance, the occupancy guard,
//! and any attached region — are the molecule's rather than merely the
//! walk's, and the backbone's bond lengths and angles are the lattice's.
//! Only hydrogens and side atoms keep the template's local geometry. See
//! `decorate` for why, and `saw::DiamondLattice::fit` for what the bond
//! length ends up being.
//!
//! The shared objective judges the decorated result at full tolerance —
//! residual contacts are reported honestly and belong to the seeded GENCAN
//! push-off (`GenCanPack::with_restart`), never hidden.

pub mod config;
pub(crate) mod decorate;
pub mod entry;
pub(crate) mod saw;

pub use config::LatticeConfig;
pub use entry::LatticeGrow;

use std::collections::HashSet;
use std::sync::Arc;

use molrs::spatial::simbox::SimBox;
use molrs::types::F;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::context::pack_state::evaluate_unscaled;
use crate::context::{PackState, Placed};
use crate::error::PackError;
use crate::grow::GrowError;
use crate::grow::internal::InternalTree;
use crate::grow::prior::TorsionPrior;
use crate::handler::{Handler, PhaseInfo, StageInfo, StepInfo};
use crate::restraint::AtomRestraint;
use crate::stage::{Budget, Guarantees, Requires, Stage, StageOutcome};
use crate::target::Target;

use decorate::{Backbone, analyze_backbone, decorate_chain};
use saw::{DiamondLattice, RisWeights, SawField, forced_zigzag, grow_walk};

struct LatticeSpecies {
    tree: InternalTree,
    backbone: Backbone,
}

/// Sites whose continuum position violates any molecule-level restraint
/// (`f(x, 1, 0.01) > 0`, same scale as CBMC). Empty input → empty mask.
fn blocked_sites(lat: &DiamondLattice, restraints: &[Arc<dyn AtomRestraint>]) -> HashSet<[i64; 3]> {
    if restraints.is_empty() {
        return HashSet::new();
    }
    lat.iter_sites()
        .filter(|&p| {
            let x = lat.to_continuum(p);
            restraints.iter().any(|r| r.f(&x, 1.0, 0.01) > 0.0)
        })
        .collect()
}

/// The diamond-lattice growth algorithm on the [`Stage`] seam.
///
/// The per-species trees and backbones are the stage's own configuration and
/// stay on it across runs, as the seam's re-entrancy contract requires.
pub struct LatticeStage {
    seed: u64,
    cfg: LatticeConfig,
    species: Vec<LatticeSpecies>,
    /// The box this stage tiles, and the radius up-scaling its cell grid is
    /// sized from. `None` only for a stage built directly from templates and
    /// never handed a cell — it then tiles whatever box the state already
    /// carries; the entry always supplies one.
    cell: Option<(SimBox, F)>,
}

impl LatticeStage {
    /// The name this stage reports, in one place: [`Stage::name`] returns it
    /// and `StepInfo.stage.name` is filled from it, so the two cannot drift
    /// apart.
    pub(crate) const NAME: &'static str = "lattice";

    /// Build the per-species trees and backbones. The `usize` names the
    /// offending target.
    ///
    /// Each tree is compiled from that target's [`Target::special_bonds`].
    pub fn from_targets(
        targets: &[Target],
        cfg: &LatticeConfig,
        seed: u64,
    ) -> Result<Self, (usize, GrowError)> {
        let mut species = Vec::with_capacity(targets.len());
        for (i, t) in targets.iter().enumerate() {
            let tree = crate::grow::tree_from_target(t).map_err(|e| (i, e))?;
            let frame = t
                .template
                .as_ref()
                .expect("tree_from_target requires a template");
            let backbone =
                analyze_backbone(frame, &tree, &t.hydrogen_mask()).map_err(|e| (i, e))?;
            species.push(LatticeSpecies { tree, backbone });
        }
        Ok(Self {
            seed,
            cfg: cfg.clone(),
            species,
            cell: None,
        })
    }

    /// The box to tile, with the `discale` its cell grid is sized from.
    ///
    /// Installed at the top of [`run`](Stage::run) rather than here: a stage
    /// is re-entrant, and the box is state the run owns, not configuration
    /// the stage consumes.
    pub(crate) fn with_resolved_cell(mut self, cell: SimBox, discale: F) -> Self {
        self.cell = Some((cell, discale));
        self
    }

    /// trans/gauche weights from the torsion prior: `States` sums the
    /// weight of states in the trans basin (|φ| > π/2, matching the
    /// calibration convention trans = π, gauche = ±π/3); the continuous
    /// priors degrade to uniform thirds.
    fn ris_weights(&self) -> RisWeights {
        match &self.cfg.torsion_prior {
            TorsionPrior::States(states) => {
                let total: F = states.iter().map(|&(_, w)| w.max(0.0)).sum();
                if total <= 0.0 {
                    return RisWeights { p_t: 1.0, p_g: 1.0 };
                }
                let trans: F = states
                    .iter()
                    .filter(|&&(a, _)| a.abs() > std::f64::consts::FRAC_PI_2 as F)
                    .map(|&(_, w)| w.max(0.0))
                    .sum();
                let p_t = trans / total;
                RisWeights {
                    p_t,
                    p_g: ((1.0 - p_t) / 2.0).max(1e-12),
                }
            }
            TorsionPrior::Uniform | TorsionPrior::Template { .. } => RisWeights {
                p_t: 1.0 / 3.0,
                p_g: 1.0 / 3.0,
            },
        }
    }
}

impl Stage for LatticeStage {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    /// Nothing: the walk builds its own placements on the lattice.
    fn requires(&self) -> Requires {
        Requires::new(Placed::None)
    }

    /// Every free molecule placed.
    fn guarantees(&self) -> Guarantees {
        Guarantees::new(Placed::All)
    }

    #[allow(clippy::needless_range_loop)]
    fn run(
        &mut self,
        state: &mut PackState,
        targets: &[Target],
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
        let origin: [F; 3] = {
            let v = sys.simbox.origin_view();
            [v[0], v[1], v[2]]
        };
        let lengths: [F; 3] = {
            let v = sys.simbox.lengths();
            [v[0], v[1], v[2]]
        };

        // One lattice for the whole box, its constant set by the mean
        // backbone bond across species (copy-count weighted).
        let (mut bond_sum, mut bond_n) = (0.0 as F, 0usize);
        for (itype, sp) in self.species.iter().enumerate() {
            let copies = sys.nmols[itype];
            bond_sum += sp.backbone.mean_bond * copies as F;
            bond_n += copies;
        }
        let mean_bond = if bond_n > 0 {
            bond_sum / bond_n as F
        } else {
            1.53
        };
        let lat = DiamondLattice::fit(origin, lengths, mean_bond);

        let weights = self.ris_weights();
        let mut field = SawField::new();
        let mut rng = SmallRng::seed_from_u64(self.seed);

        let mut degraded = 0usize;
        let mut aborted = false;

        // Bases of the chains already walked + decorated, in xcart order:
        // their coordinates live in `sys.xcart` (the lab-frame home the
        // writeback reads), so the driver only tracks which copies are done.
        let mut done: Vec<usize> = Vec::with_capacity(sys.ntotmol);
        let mut mol = 0usize;
        'outer: for itype in 0..sys.ntype {
            let sp = &self.species[itype];
            let na = sys.natoms[itype];
            let mask = blocked_sites(&lat, &targets[itype].molecule_restraints);
            if !lat.has_allowed_a_site(&mask) {
                return Err(PackError::Grow {
                    target: itype,
                    source: GrowError::LatticeRegionEmpty,
                });
            }
            field.set_blocked(mask);
            for imol in 0..sys.nmols[itype] {
                let base = sys.idfirst[itype] + imol * na;
                // Escape ladder: guarded walk → unguarded walk (site
                // self-avoidance only) → forced zigzag. Every escape below
                // the guard is a relaxation of the constructive guarantee
                // and is counted.
                let mut relaxed = 0usize;
                let walk = if let Some(w) = grow_walk(
                    &lat,
                    &mut field,
                    mol as u32,
                    &sp.backbone.parent,
                    &sp.backbone.children,
                    &sp.backbone.follows,
                    &weights,
                    self.cfg.occupancy_guard,
                    config::MAX_BACKTRACK,
                    config::MAX_RESEED,
                    &mut rng,
                ) {
                    w
                } else if let Some(w) = grow_walk(
                    &lat,
                    &mut field,
                    mol as u32,
                    &sp.backbone.parent,
                    &sp.backbone.children,
                    &sp.backbone.follows,
                    &weights,
                    false,
                    config::MAX_BACKTRACK,
                    config::MAX_RESEED,
                    &mut rng,
                ) {
                    relaxed = 1;
                    w
                } else {
                    relaxed = 1;
                    forced_zigzag(
                        &lat,
                        &mut field,
                        mol as u32,
                        &sp.backbone.parent,
                        &sp.backbone.children,
                        &sp.backbone.follows,
                        &mut rng,
                    )
                    .ok_or(PackError::Grow {
                        target: itype,
                        source: GrowError::LatticeRegionEmpty,
                    })?
                };
                degraded += relaxed;

                let mut coords = vec![[0.0 as F; 3]; na];
                decorate_chain(&sp.tree, &sp.backbone, &lat, &walk.sites, &mut coords);
                for (a, p) in coords.iter().enumerate() {
                    sys.xcart[base + a] = *p;
                }
                done.push(base);
                mol += 1;

                // Handler visibility: one StepInfo per finished chain.
                let info = StepInfo {
                    // One stage per run until the pipeline lands; the name
                    // comes from the stage type so the two cannot drift.
                    stage: StageInfo {
                        index: 0,
                        total: 1,
                        name: Self::NAME,
                    },
                    loop_idx: mol,
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
                    radscale: 1.0,
                    precision: budget.precision,
                };
                for h in handlers.iter_mut() {
                    h.on_step(&info, sys);
                }
                if handlers.iter().any(|h| h.should_stop()) {
                    aborted = true;
                    break 'outer;
                }
            }
        }

        // An abort leaves later chains unplaced: complete them as forced
        // zigzags so assembly stays chemical (same contract as the
        // continuum driver's abort path).
        if aborted {
            let mut m = done.len();
            for itype in 0..sys.ntype {
                let sp = &self.species[itype];
                let na = sys.natoms[itype];
                field.set_blocked(blocked_sites(&lat, &targets[itype].molecule_restraints));
                for imol in 0..sys.nmols[itype] {
                    let base = sys.idfirst[itype] + imol * na;
                    if done.contains(&base) {
                        continue;
                    }
                    let walk = forced_zigzag(
                        &lat,
                        &mut field,
                        m as u32,
                        &sp.backbone.parent,
                        &sp.backbone.children,
                        &sp.backbone.follows,
                        &mut rng,
                    )
                    .ok_or(PackError::Grow {
                        target: itype,
                        source: GrowError::LatticeRegionEmpty,
                    })?;
                    let mut coords = vec![[0.0 as F; 3]; na];
                    decorate_chain(&sp.tree, &sp.backbone, &lat, &walk.sites, &mut coords);
                    for (a, p) in coords.iter().enumerate() {
                        sys.xcart[base + a] = *p;
                    }
                    done.push(base);
                    m += 1;
                }
            }
        }

        // Writeback contract: per-copy COM into x, centered conformer into
        // this copy's own `coor` block — one derivation, shared with the
        // continuum driver. Every chain wrote its coordinates into
        // `sys.xcart` as it finished, including the forced completions above.
        x.capture_from_xcart(sys);

        // Final verdict from the shared objective, never self-reported: the
        // same unscaled primitive the continuum driver calls, which also owns
        // the `scale` / `scale2` handling this site used to spell out.
        let (_, fdist, frest) = evaluate_unscaled(sys, x.as_slice());
        let converged = !aborted && degraded == 0 && fdist == 0.0 && frest < budget.precision;
        Ok(StageOutcome::new(converged, degraded))
    }
}
