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
//! [`Stage`](crate::stage::Stage) seam is the only shared contract.
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
//! in `degraded`), chain statistics governed by the priors, and homogeneous
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
#[cfg(test)]
mod tests;

pub use config::{GrowConfig, GrowError};
pub use driver::GrowStage;
pub use prior::{AnglePrior, TorsionPrior};

use molrs::store::frame::Frame;
use molrs::types::F;

use crate::grow::internal::InternalTree;
use crate::target::Target;

/// Read the template's bond graph and coordinates for growth.
///
/// Coordinates are in Å. This is the only grow-module caller of
/// `crate::template`, and the unique implementer of the refusal order
/// (including `RingTemplate`):
/// `NoAtomsBlock → BondOutOfRange → NoBonds → TemplateTooSmall →
/// Disconnected → RingTemplate`. Bond graphs come from
/// `molrs::Topology::from_frame`; this function does not rebuild CSR.
pub(crate) fn topology_for_growth(
    frame: &Frame,
) -> Result<(molrs::Topology, Vec<[F; 3]>), GrowError> {
    let xyz = crate::template::coord_rows(&frame.coords().map_err(|_| GrowError::NoAtomsBlock)?);
    let n = xyz.len();
    if let Some(err) = first_bond_out_of_range(frame, n) {
        return Err(err);
    }
    // molrs treats a missing or empty bonds block as Ok with zero edges.
    // A still-failing `from_frame` (non-empty bonds missing atomi/atomj)
    // is named `NoBonds`; MolRsError is never wrapped.
    let topo = molrs::Topology::from_frame(frame).map_err(|_| GrowError::NoBonds)?;
    if topo.n_bonds() == 0 {
        return Err(GrowError::NoBonds);
    }
    if n < 3 {
        return Err(GrowError::TemplateTooSmall(n));
    }
    if topo.n_components() != 1 {
        return Err(GrowError::Disconnected);
    }
    // A tree has n−1 edges; e ≥ n is a cycle.
    if topo.n_bonds() >= n {
        return Err(GrowError::RingTemplate);
    }
    Ok((topo, xyz))
}

fn first_bond_out_of_range(frame: &Frame, n: usize) -> Option<GrowError> {
    let block = frame.get("bonds")?;
    let atomi = block.get_uint("atomi")?;
    let atomj = block.get_uint("atomj")?;
    for (&a, &b) in atomi.iter().zip(atomj.iter()) {
        let a = a as usize;
        let b = b as usize;
        if a >= n || b >= n {
            return Some(GrowError::BondOutOfRange { a, b, n });
        }
    }
    None
}

/// Validate that a template frame can seed growth.
///
/// Named rejections only — a target that cannot grow is refused, never
/// silently degraded to rigid-body packing (the caller picks the method per
/// target, and molpack does not second-guess it).
pub(crate) fn validate_template(frame: Option<&Frame>) -> Result<(), GrowError> {
    let frame = frame.ok_or(GrowError::MissingTemplate)?;
    topology_for_growth(frame).map(|_| ())
}

/// Compile the target's intramolecular skip table into an internal-coordinate
/// tree. Missing template is [`GrowError::MissingTemplate`]; a non-binary
/// table is [`GrowError::NonBinarySpecialBond`].
pub(crate) fn tree_from_target(t: &Target) -> Result<InternalTree, GrowError> {
    let frame = t.template.as_ref().ok_or(GrowError::MissingTemplate)?;
    InternalTree::from_frame(frame, &t.special_bonds)
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
    // molrs's own classification — the one minimum image and the cell grid
    // use — so growth never accepts a box the rest of the run treats as tilted.
    if !matches!(simbox.kind(), molrs::spatial::simbox::BoxKind::Ortho { .. }) {
        return Err(crate::error::PackError::Grow {
            target,
            source: crate::grow::GrowError::TriclinicCell,
        });
    }
    Ok(simbox)
}
