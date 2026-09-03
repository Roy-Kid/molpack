//! Chain growth: a configurational-bias packing algorithm ranked alongside
//! GENCAN.
//!
//! Rigid-body placement + GENCAN descent solves packing problems whose
//! molecular shapes are fixed before packing begins (water, urea, lipids,
//! proteins). A dense polymer melt is not such a problem: at melt density the
//! box is covered by many interpenetrating chains, and no downhill search
//! over rigid placements reaches an interdigitated state. A melt must be
//! **grown** into place, one torsion at a time, inside the final box.
//!
//! Growth is therefore a *peer* of the GENCAN path, not a component of it:
//! it consumes the same [`PackContext`](crate::PackContext) (same radii, same
//! restraints, same cell), is judged by the same objective, and is selected
//! by the [`CbmcGrow`](crate::CbmcGrow) entry. It never
//! calls the GENCAN internals, and the GENCAN path never calls it — the
//! [`Solver`](crate::solver::Solver) seam is the only shared contract.
//!
//! The module is purely geometric: no force-field dependency, no `ff`
//! feature. Conformer statistics come from user-supplied geometric priors
//! ([`prior::TorsionPrior`] / [`prior::AnglePrior`]).
//!
//! # What growth delivers — and what it does not
//!
//! Growth delivers **geometry**: no contact below the declared tolerance
//! (constructive — a placement violating the hard core or a restraint is
//! rejected, never penalized; every relaxation of that guarantee is counted
//! in `softened`), chain statistics governed by the priors, and homogeneous
//! density. It does **not** deliver an equilibrium Boltzmann ensemble: greedy
//! Rosenbluth selection has a known, characterizable bias (Consta et al.
//! 1999), and equilibration is the downstream MD's job — the classic
//! generate → push-off → equilibrate pipeline (Auhl et al. 2003). When the
//! hard core had to soften at melt density, the packer itself chains into
//! the GENCAN phases as the push-off stage, on the same state, with zero
//! coordinate conversion.
//!
//! # Round-snapshot semantics (part of the algorithm, not an implementation
//! detail)
//!
//! All chains advance one step per round. Within a round every chain's
//! candidates are proposed and scored against the round-start snapshot of
//! the overlap field; commits run serially in the round's shuffled order,
//! re-validating against this round's earlier commits. Randomness comes from
//! counter-based streams hashed per `(seed, molecule, stage, visit)`. Both
//! choices exist so that a future parallel driver (parallel proposals +
//! serial commits) is bit-identical to this serial one — see
//! [`driver`] for the full contract.

pub mod config;
pub mod driver;
pub mod entry;
pub mod field;
pub mod internal;
pub mod lattice;
pub(crate) mod moves;
pub mod prior;

pub use config::{GrowConfig, GrowError};
pub use driver::GrowthSolver;
pub use prior::{AnglePrior, TorsionPrior};

use molrs::store::frame::Frame;

use crate::topology::Topology;

/// Validate that a template frame can seed growth.
///
/// Named rejections only — a target that cannot grow is refused, never
/// silently degraded to rigid-body packing (the caller picks the method per
/// target, and molpack does not second-guess it).
pub(crate) fn validate_template(frame: Option<&Frame>) -> Result<(), GrowError> {
    let frame = frame.ok_or(GrowError::MissingTemplate)?;
    let topo = Topology::from_frame(frame).map_err(GrowError::Topology)?;
    if topo.natoms() < 3 {
        return Err(GrowError::TemplateTooSmall(topo.natoms()));
    }
    // A connected tree has exactly n-1 edges; e >= n implies a cycle. (An
    // e == n-1 graph with a cycle is disconnected and named separately.)
    if topo.bonds().len() >= topo.natoms() {
        return Err(GrowError::RingTemplate);
    }
    Ok(())
}

/// Growth needs a resolved, orthorhombic box: the v1 overlap field supports
/// nothing else, and the box must exist from the first atom. Named errors
/// per spec principle 3.
pub(crate) fn validate_grow_cell(
    cell: Option<molrs::spatial::simbox::SimBox>,
    target: usize,
) -> Result<molrs::spatial::simbox::SimBox, crate::error::PackError> {
    let Some(simbox) = cell else {
        return Err(crate::error::PackError::Grow {
            target,
            source: crate::grow::GrowError::NoBox,
        });
    };
    let h = simbox.h_view();
    let ortho = (0..3).all(|i| (0..3).all(|j| i == j || h[(i, j)].abs() < 1e-9));
    if !ortho {
        return Err(crate::error::PackError::Grow {
            target,
            source: crate::grow::GrowError::TriclinicCell,
        });
    }
    Ok(simbox)
}
