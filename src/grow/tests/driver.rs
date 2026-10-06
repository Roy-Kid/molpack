//! The growth driver's two constructive guarantees, at the level where they
//! are observable: a bounded round loop, and an honest verdict when that
//! bound is what ended the run.
//!
//! The ladder predicates themselves (`rung_due`, `retract_depth`,
//! `force_due`, `max_rounds`) are pinned next to their definitions in
//! `src/grow/driver.rs`. What cannot be read off a predicate is the
//! *surrender path*: on reaching the cap the pending chains are force-
//! completed, every forced placement is counted in `degraded`, and the
//! outcome is `converged == false`. The fixture below is the smallest one
//! that reaches the cap — a hard core that no placement can satisfy, and a
//! one-pass budget — so it costs milliseconds, not a melt.

use super::*;

use crate::AtomRestraint;
use crate::context::PackContext;
use crate::handler::{Handler, StepInfo};

/// A restraint that refuses every point, so every growth attempt is a dead end.
#[derive(Debug)]
struct RefuseEverywhere;

impl AtomRestraint for RefuseEverywhere {
    fn f(&self, _x: &[F; 3], _scale: F, _scale2: F) -> F {
        1.0
    }
    fn fg(&self, _x: &[F; 3], _scale: F, _scale2: F, _g: &mut [F; 3]) -> F {
        1.0
    }
}

/// Stops the driver after the first round, before a second chain could
/// spend a rung that the old global ladder had queued.
#[derive(Default)]
struct StopAfterOne {
    seen: bool,
}

impl Handler for StopAfterOne {
    fn on_step(&mut self, _info: &StepInfo, _sys: &PackContext) {
        self.seen = true;
    }
    fn should_stop(&self) -> bool {
        self.seen
    }
}

/// A density the strict hard core cannot satisfy must still *return*, with a
/// verdict that says so. Without the cap the round loop spins forever
/// (debt D-01 (ii)); with it, the run is a surrender, never a silent success.
#[test]
fn an_unsatisfiable_hard_core_terminates_and_says_so() {
    let (copies, n_beads, bond, l) = (4usize, 5usize, 1.53 as F, 6.0 as F);
    // Strict core (no softening rung is ever earned, none is ever taken), so
    // the only way out of the round loop is the cap.
    let cfg = GrowConfig::new(TorsionPrior::Uniform)
        .with_min_hard_scale(1.0)
        .with_soften_after(usize::MAX / 4);
    let target = Target::new(chain_frame(n_beads, bond, true), copies);

    let state = CbmcGrow::from_config(cfg)
        .with_seed(11)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [l; 3], [true; 3])
        .run(&[target], 1)
        .expect(
            "an impossible density is a result, not an error: the run must \
             return Ok with converged == false",
        );

    assert!(
        !state.converged,
        "a run that forced its way out of an unsatisfiable hard core must \
         report converged == false — the cap is a surrender, not a proof",
    );
    assert!(
        state.degraded > 0,
        "degraded = 0 after a capped run — every forced completion breaks the \
         constructive guarantee and must be counted, or the caller cannot \
         tell a capped result from a clean one",
    );

    // A surrender still owes the caller a complete, finite structure.
    let pos = state.positions();
    assert_eq!(pos.len(), copies * n_beads, "every bead is written back");
    for (i, p) in pos.iter().enumerate() {
        assert!(
            p.iter().all(|v| v.is_finite()),
            "atom {i} at {p:?} — a capped run must still write real coordinates",
        );
    }
}

/// Two chains that both earn a rung in the same round each shrink their own
/// core. The old ladder took one global rung per round, so the second chain
/// waited and `degraded` counted one shrink.
#[test]
fn a_rung_shrinks_only_the_chain_that_earned_it() {
    let cfg = GrowConfig::new(TorsionPrior::Uniform)
        .with_soften_after(1)
        .with_min_hard_scale(GrowConfig::SOFTEN_RUNG);
    let frame = chain_frame(3, 1.53, true);
    let wedged = |frame| Target::new(frame, 1).with_restraint(RefuseEverywhere);
    let state = CbmcGrow::from_config(cfg)
        .with_handler(Box::new(StopAfterOne::default()))
        .with_seed(3)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [20.0; 3], [true; 3])
        .run(&[wedged(frame.clone()), wedged(frame)], 1)
        .expect("both chains are refused, and the run still returns");

    // One round, then the handler stops. A 3-bead template is a seed plus
    // one step, so the abort force-places both for each chain (4) and each
    // chain has taken its own rung (2). A shared core would count 5.
    assert_eq!(
        state.degraded, 6,
        "each chain's rung is its own — a shared core would count 5 \
         (one rung, four forced placements)",
    );
    assert!(!state.converged);
}
