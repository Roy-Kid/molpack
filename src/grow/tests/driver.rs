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
