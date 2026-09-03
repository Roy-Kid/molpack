//! Lattice growth: diamond-lattice SAW pre-generation for dense melts
//! (lattice-growth-phase spec).
//!
//! The continuum growth solver hits a density ceiling long before melt
//! density: its torsion candidates are proposed in continuous space and
//! ground down by retraction. On the diamond lattice the same RIS states
//! are EXACT lattice moves — trans/gauche± are the three non-backtracking
//! continuations — and excluded volume is an O(1) site-occupancy check, so
//! a melt-density walk completes in milliseconds where continuum growth
//! grinds. The walk decides only the torsion sequence; decoration rebuilds
//! every atom from the template's true internal coordinates
//! ([`decorate`]), and the shared objective judges the decorated result at
//! full tolerance — residual contacts are reported honestly and belong to
//! the seeded GENCAN push-off (`GenCanPack::seeded_from`), never hidden.

pub mod config;
pub(crate) mod decorate;
pub(crate) mod saw;

pub use config::LatticeConfig;

use molrs::types::F;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::constraints::EvalMode;
use crate::context::{PackContext, RigidView};
use crate::entry::{EngineSetup, PackEngine, PackSettings};
use crate::error::PackError;
use crate::grow::internal::InternalTree;
use crate::grow::prior::TorsionPrior;
use crate::grow::{GrowError, validate_grow_cell};
use crate::handler::{Handler, PhaseInfo, StepInfo};
use crate::solver::{Budget, SolveOutcome, Solver};
use crate::target::Target;

use decorate::{Backbone, analyze_backbone, decorate_chain};
use saw::{DiamondLattice, RisWeights, SawField, forced_zigzag, grow_walk};

struct LatticeSpecies {
    tree: InternalTree,
    backbone: Backbone,
}

/// The diamond-lattice growth algorithm on the [`Solver`] seam.
pub struct LatticeSolver {
    seed: u64,
    cfg: LatticeConfig,
    species: Vec<LatticeSpecies>,
}

impl LatticeSolver {
    /// Build the per-species trees and backbones. The `usize` names the
    /// offending target.
    pub fn from_targets(
        targets: &[Target],
        cfg: &LatticeConfig,
        seed: u64,
    ) -> Result<Self, (usize, GrowError)> {
        let mut species = Vec::with_capacity(targets.len());
        for (i, t) in targets.iter().enumerate() {
            let frame = t.template.as_ref().ok_or((i, GrowError::MissingTemplate))?;
            let tree = InternalTree::from_frame(frame).map_err(|e| (i, e))?;
            let backbone = analyze_backbone(frame, &tree).map_err(|e| (i, e))?;
            species.push(LatticeSpecies { tree, backbone });
        }
        Ok(Self {
            seed,
            cfg: cfg.clone(),
            species,
        })
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

impl Solver for LatticeSolver {
    fn name(&self) -> &'static str {
        "lattice"
    }

    fn solve(
        &mut self,
        sys: &mut PackContext,
        _targets: &[Target],
        x: &mut RigidView,
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> SolveOutcome {
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

        let mut softened = 0usize;
        let mut aborted = false;

        // Bases of the chains already walked + decorated, in xcart order:
        // their coordinates live in `sys.xcart` (the lab-frame home the
        // writeback reads), so the driver only tracks which copies are done.
        let mut done: Vec<usize> = Vec::with_capacity(sys.ntotmol);
        let mut mol = 0usize;
        'outer: for itype in 0..sys.ntype {
            let sp = &self.species[itype];
            let na = sys.natoms[itype];
            let n_bb = sp.backbone.atoms.len();
            for imol in 0..sys.nmols[itype] {
                let base = sys.idfirst[itype] + imol * na;
                // Escape ladder: guarded walk → unguarded walk (site
                // self-avoidance only) → forced zigzag. Every escape below
                // the guard is a relaxation of the constructive guarantee
                // and is counted.
                let walk = if let Some(w) = grow_walk(
                    &lat,
                    &mut field,
                    mol as u32,
                    n_bb,
                    &weights,
                    self.cfg.occupancy_guard,
                    self.cfg.max_backtrack,
                    self.cfg.max_reseed,
                    &mut rng,
                ) {
                    w
                } else if let Some(w) = grow_walk(
                    &lat,
                    &mut field,
                    mol as u32,
                    n_bb,
                    &weights,
                    false,
                    self.cfg.max_backtrack,
                    self.cfg.max_reseed,
                    &mut rng,
                ) {
                    softened += 1;
                    w
                } else {
                    softened += 1;
                    forced_zigzag(&lat, &mut field, mol as u32, n_bb, &mut rng)
                };

                let mut coords = vec![[0.0 as F; 3]; na];
                decorate_chain(
                    &sp.tree,
                    &sp.backbone,
                    &lat,
                    &walk.sites,
                    self.cfg.track_tweak,
                    &mut coords,
                );
                for (a, p) in coords.iter().enumerate() {
                    sys.xcart[base + a] = *p;
                }
                done.push(base);
                mol += 1;

                // Handler visibility: one StepInfo per finished chain.
                let info = StepInfo {
                    loop_idx: mol,
                    max_loops: budget.max_loops,
                    phase: PhaseInfo {
                        phase: 0,
                        total_phases: 1,
                        molecule_type: None,
                    },
                    fdist: 0.0,
                    frest: 0.0,
                    improvement_pct: 0.0,
                    radscale: 1.0,
                    precision: budget.precision,
                    relaxer_acceptance: Vec::new(),
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
                let n_bb = sp.backbone.atoms.len();
                for imol in 0..sys.nmols[itype] {
                    let base = sys.idfirst[itype] + imol * na;
                    if done.contains(&base) {
                        continue;
                    }
                    let walk = forced_zigzag(&lat, &mut field, m as u32, n_bb, &mut rng);
                    let mut coords = vec![[0.0 as F; 3]; na];
                    decorate_chain(
                        &sp.tree,
                        &sp.backbone,
                        &lat,
                        &walk.sites,
                        self.cfg.track_tweak,
                        &mut coords,
                    );
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

        // Final verdict from the shared objective, never self-reported.
        sys.scale = 1.0;
        sys.scale2 = 0.01;
        let _ = sys.evaluate(x.as_slice(), EvalMode::FOnly, None);
        let converged =
            !aborted && softened == 0 && sys.fdist == 0.0 && sys.frest < budget.precision;
        SolveOutcome::new(converged, sys.fdist, sys.frest, softened)
    }
}

/// Diamond-lattice growth as its own entry (lattice-growth-phase spec).
///
/// The torsion prior is mandatory — it decides the walk's trans/gauche±
/// weights — so it is the one constructor argument. The entry reports its
/// outcome honestly: decoration drift and hydrogen crowding leave real
/// contacts at melt density, `fdist` says so, and the remedy is the
/// explicit seeded push-off chain
/// ([`GenCanPack::seeded_from`](crate::GenCanPack::seeded_from)).
pub struct LatticeGrow {
    settings: PackSettings,
    handlers: Vec<Box<dyn Handler>>,
    config: LatticeConfig,
}

impl LatticeGrow {
    pub fn new(torsion_prior: TorsionPrior) -> Self {
        Self::from_config(LatticeConfig::new(torsion_prior))
    }

    /// Build the entry around an existing [`LatticeConfig`].
    pub fn from_config(config: LatticeConfig) -> Self {
        Self {
            settings: PackSettings::default(),
            handlers: Vec::new(),
            config,
        }
    }

    /// Nearest-neighbour site exclusion (default on).
    pub fn with_occupancy_guard(mut self, on: bool) -> Self {
        self.config = self.config.with_occupancy_guard(on);
        self
    }

    /// Path-tracking torsion tweak in radians (default 0.35 ≈ 20°); `0.0`
    /// disables tracking. See [`LatticeConfig::with_track_tweak`].
    pub fn with_track_tweak(mut self, radians: F) -> Self {
        self.config = self.config.with_track_tweak(radians);
        self
    }
}

impl PackEngine for LatticeGrow {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }
    fn settings_mut(&mut self) -> &mut PackSettings {
        &mut self.settings
    }
    fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>> {
        &mut self.handlers
    }

    fn validate(&self, targets: &[Target]) -> Result<(), PackError> {
        for (i, t) in targets.iter().enumerate() {
            crate::grow::validate_template(t.template.as_ref())
                .map_err(|source| PackError::Grow { target: i, source })?;
            if t.fixed_at.is_some() {
                return Err(PackError::Grow {
                    target: i,
                    source: GrowError::FixedTarget,
                });
            }
            // Backbone analyzability is part of validation: named rejection
            // before any context is built.
            let frame = t.template.as_ref().expect("validate_template checked");
            let tree = InternalTree::from_frame(frame)
                .map_err(|source| PackError::Grow { target: i, source })?;
            analyze_backbone(frame, &tree)
                .map_err(|source| PackError::Grow { target: i, source })?;
        }
        Ok(())
    }

    fn prepare(
        &self,
        sys: &mut PackContext,
        _x: &mut [F],
        setup: &EngineSetup<'_>,
    ) -> Result<(), PackError> {
        let simbox = validate_grow_cell(setup.cell.clone(), 0)?;
        let radmax = sys.radius.iter().cloned().fold(0.0 as F, F::max);
        crate::initial::install_simbox_and_grid(
            sys,
            simbox,
            radmax,
            self.settings.discale(),
            setup.ntotat_free,
        );
        Ok(())
    }

    fn solver(&mut self, setup: &EngineSetup<'_>) -> Result<Box<dyn Solver>, PackError> {
        let solver = LatticeSolver::from_targets(setup.targets, &self.config, self.settings.seed())
            .map_err(|(target, source)| PackError::Grow { target, source })?;
        Ok(Box::new(solver))
    }
}
