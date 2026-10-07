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
//! [`Stage`](crate::Stage) seam is the only shared contract.
//!
//! The module is purely geometric: no force-field dependency, no `ff`
//! feature. Conformer statistics come from user-supplied geometric priors
//! ([`TorsionPrior`] / [`AnglePrior`]).
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
//! generate → push-off → equilibrate pipeline (Auhl et al. 2003). Push-off,
//! when the caller wants it, is [`GencanPack::with_restart`](crate::GencanPack::with_restart)
//! on this run's [`State`](crate::State).
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
//! the growth driver (`grow/driver.rs`) for the full contract.

pub(crate) mod cbmc_grow;
mod config;
pub(crate) mod driver;
mod error;
pub(crate) mod field;
pub(crate) mod internal;
pub(crate) mod lattice;
pub(crate) mod moves;
mod prior;
#[cfg(test)]
mod tests;

// The entries (`CbmcGrow`, `LatticeGrow`) are published at the crate root;
// this module publishes their configuration vocabulary. Leaves are private:
// one path per item.
pub use config::GrowConfig;
pub use error::GrowError;
pub use lattice::LatticeConfig;
pub use prior::{AnglePrior, TorsionPrior};

pub(crate) use driver::GrowStage;

use molrs::core::Frame;
use molrs::op::F;

use crate::Target;
use crate::grow::internal::InternalTree;

/// Read the template's bond graph and coordinates for growth.
///
/// Coordinates are in Å. This is the only grow-module caller of
/// `crate::template`, and the unique implementer of the refusal order
/// (including `RingTemplate`):
/// `NoAtomsBlock → MissingBondEndpoint → BondOutOfRange → NoBonds →
/// TemplateTooSmall → Disconnected → RingTemplate`. Bond graphs come from
/// `molrs::core::Topology::from_frame`; this function does not rebuild CSR.
pub(crate) fn topology_for_growth(
    frame: &Frame,
) -> Result<(molrs::core::Topology, Vec<[F; 3]>), GrowError> {
    let xyz = crate::template::coord_rows(&frame.coords().map_err(|_| GrowError::NoAtomsBlock)?);
    let n = xyz.len();
    // molrs treats a missing or empty bonds block as Ok with zero edges.
    // A still-failing `from_frame` names the column or the row; MolRsError
    // is never wrapped.
    let topo = match molrs::core::Topology::from_frame(frame) {
        Ok(topo) => topo,
        Err(molrs::core::TopologyError::EndpointOutOfRange {
            row, atom, n_atoms, ..
        }) => {
            return Err(GrowError::BondOutOfRange {
                row,
                atom,
                n: n_atoms,
            });
        }
        Err(molrs::core::TopologyError::MissingEndpoint { column, .. }) => {
            return Err(GrowError::MissingBondEndpoint { column });
        }
        Err(
            molrs::core::TopologyError::MissingBlock { .. }
            | molrs::core::TopologyError::NoRows { .. },
        ) => return Err(GrowError::NoAtomsBlock),
    };
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
    cell: Option<molrs::core::SimBox>,
) -> Result<molrs::core::SimBox, GrowError> {
    let Some(simbox) = cell else {
        return Err(GrowError::NoBox);
    };
    // molrs's own classification — the one minimum image and the cell grid
    // use — so growth never accepts a box the rest of the run treats as tilted.
    if !matches!(simbox.kind(), molrs::core::BoxKind::Ortho { .. }) {
        return Err(GrowError::TriclinicCell);
    }
    Ok(simbox)
}
