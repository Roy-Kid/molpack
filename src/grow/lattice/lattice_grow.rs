//! `LatticeGrow` — the diamond-lattice chain-growth entry.
//!
//! Split out of `lattice/mod.rs` so the lattice module has the same shape as
//! its two peers (`gencan/gencan_pack.rs`, `grow/cbmc_grow.rs`): the algorithm in one
//! file, the entry that selects it in another.

use crate::Callback;
use crate::PackError;
use crate::Stage;
use crate::Target;
use crate::grow::prior::TorsionPrior;
use crate::grow::{GrowError, validate_grow_cell};
use crate::pipeline::{EngineSetup, PackEngine, Pipeline, StageFactory};
use crate::{PackSettings, State};

use super::LatticeStage;
use super::config::LatticeConfig;
use super::decorate::analyze_backbone;

/// Diamond-lattice growth as its own entry (lattice-growth-phase spec).
///
/// A tetrahedral heavy-atom tree of degree ≤ 4 is grown as a diamond-lattice
/// self-avoiding walk; a linear chain is the `d = 2` degeneracy of that
/// walk, not a second algorithm. Degree `> 4` is
/// [`GrowError::NonTetrahedralTemplate`]. Rings stay
/// [`GrowError::RingTemplate`].
///
/// The torsion prior is mandatory — it decides the walk's trans/gauche±
/// weights — so it is the one constructor argument. The entry reports its
/// outcome honestly: decoration drift and hydrogen crowding leave real
/// contacts at melt density, `fdist` says so, and the remedy is the
/// explicit seeded push-off chain
/// ([`GencanPack::with_restart`](crate::GencanPack::with_restart)).
///
/// A molecule-level geometric restraint (including
/// any molrs region lifted by [`RegionRestraint`](crate::RegionRestraint)) masks the diamond lattice: sites whose
/// continuum position lies outside the region are blocked, and the SAW
/// only walks Region ∩ lattice. An empty intersection is
/// [`GrowError::LatticeRegionEmpty`].
pub struct LatticeGrow {
    settings: PackSettings,
    callbacks: Vec<Box<dyn Callback>>,
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
            callbacks: Vec::new(),
            config,
        }
    }

    /// Nearest-neighbour site exclusion (default on).
    pub fn with_occupancy_guard(mut self, on: bool) -> Self {
        self.config = self.config.with_occupancy_guard(on);
        self
    }
}

impl StageFactory for LatticeGrow {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }

    /// Each tree is compiled from that target's [`Target::special_bonds`].
    fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> {
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
            // before any system is built.
            let frame = t.template.as_ref().expect("validate_template checked");
            let tree = crate::grow::tree_from_target(t)
                .map_err(|source| PackError::Grow { target: i, source })?;
            analyze_backbone(frame, &tree, &t.hydrogen_mask())
                .map_err(|source| PackError::Grow { target: i, source })?;
        }
        Ok(())
    }

    fn take_callbacks(&mut self) -> Vec<Box<dyn Callback>> {
        std::mem::take(self.callbacks_mut())
    }

    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        // The lattice needs a box to tile — a named error, resolved here so
        // it is reported before any stage runs. The stage installs it (and
        // the matching cell grid) at the top of its own run.
        let stage = LatticeStage::from_targets(setup.targets, &self.config, setup.settings.seed())
            .map_err(|(target, source)| PackError::Grow { target, source })?;
        let cell = validate_grow_cell(setup.cell.clone())
            .map_err(|source| PackError::Grow { target: 0, source })?;
        Ok(vec![Box::new(
            stage.with_resolved_cell(cell, setup.settings.discale()),
        )])
    }
}

impl PackEngine for LatticeGrow {
    fn settings_mut(&mut self) -> &mut PackSettings {
        &mut self.settings
    }
    fn callbacks_mut(&mut self) -> &mut Vec<Box<dyn Callback>> {
        &mut self.callbacks
    }

    fn run(self, targets: &[Target], max_loops: usize) -> Result<State, PackError> {
        Pipeline::single(self).run(targets, max_loops)
    }
}
