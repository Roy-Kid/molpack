//! `LatticeGrow` — the diamond-lattice chain-growth entry.
//!
//! Split out of `lattice/mod.rs` so the lattice module has the same shape as
//! its two peers (`gencan/entry.rs`, `grow/entry.rs`): the algorithm in one
//! file, the entry that selects it in another.

use molrs::types::F;

use crate::context::PackContext;
use crate::entry::{EngineSetup, PackEngine, PackSettings};
use crate::error::PackError;
use crate::grow::internal::InternalTree;
use crate::grow::prior::TorsionPrior;
use crate::grow::{GrowError, validate_grow_cell};
use crate::handler::Handler;
use crate::stage::Stage;
use crate::target::Target;

use super::LatticeStage;
use super::config::LatticeConfig;
use super::decorate::analyze_backbone;

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

    fn solver(&mut self, setup: &EngineSetup<'_>) -> Result<Box<dyn Stage>, PackError> {
        let stage = LatticeStage::from_targets(setup.targets, &self.config, self.settings.seed())
            .map_err(|(target, source)| PackError::Grow { target, source })?;
        Ok(Box::new(stage))
    }
}
