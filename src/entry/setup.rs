//! Shared pack-space setup: density / periodic-box / cell resolution and
//! global-restraint broadcast.
//!
//! Every entry needs the same answers before any algorithm runs: what box
//! does the system live in, and which restraints apply to every target.
//! This machinery moved out of the packer (engine-entry-split) so the entry
//! lifecycle owns it once — never re-declared per algorithm.

use std::sync::Arc;

use molrs::Element;
use molrs::spatial::simbox::SimBox;
use molrs::types::F;
use ndarray::array;

use crate::error::PackError;
use crate::restraint::AtomRestraint;
use crate::target::Target;

pub(crate) type PeriodicSpec = ([F; 3], [F; 3], [bool; 3]);

/// A packing cell as the caller declared it, resolved to a [`SimBox`]
/// inside the engine lifecycle so the builders can stay infallible.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum CellDecl {
    LengthsAngles {
        lengths: [F; 3],
        angles_deg: [F; 3],
        pbc: [bool; 3],
    },
    Matrix {
        h: [[F; 3]; 3],
        origin: [F; 3],
        pbc: [bool; 3],
    },
}

impl CellDecl {
    pub(crate) fn resolve(self) -> Result<SimBox, PackError> {
        let (h, origin, pbc) = match self {
            CellDecl::LengthsAngles {
                lengths,
                angles_deg,
                pbc,
            } => (
                SimBox::matrix_from_lengths_angles(lengths, angles_deg).map_err(|_| {
                    PackError::InvalidCell {
                        detail: format!(
                            "lengths {lengths:?} and angles {angles_deg:?} do not describe a cell"
                        ),
                    }
                })?,
                [0.0; 3],
                pbc,
            ),
            CellDecl::Matrix { h, origin, pbc } => (
                array![
                    [h[0][0], h[0][1], h[0][2]],
                    [h[1][0], h[1][1], h[1][2]],
                    [h[2][0], h[2][1], h[2][2]]
                ],
                origin,
                pbc,
            ),
        };
        SimBox::new(h, array![origin[0], origin[1], origin[2]], pbc).map_err(|_| {
            PackError::InvalidCell {
                detail: "lattice matrix is singular".to_string(),
            }
        })
    }
}

/// Reject a half-space restraint declared across a periodic lattice direction.
///
/// A plane splits space in two. Under periodicity along lattice vector `a`,
/// translating a point by `a` must leave it on the same side, which holds only
/// when the normal is orthogonal to `a`. Packmol evaluates such a constraint
/// anyway, in the frame of whichever cell the origin sits in, so whether it is
/// satisfied depends on where the user put the origin. There is no
/// interpretation to salvage — the declaration is refused, naming the axis.
pub(crate) fn reject_planes_across_periodic_axes(
    targets: &[Target],
    bx: &SimBox,
) -> Result<(), PackError> {
    let pbc = bx.pbc();
    if !pbc.iter().any(|&p| p) {
        return Ok(());
    }
    let h = bx.h_view();
    // Lattice vectors are the columns of H.
    let lattice = |k: usize| [h[[0, k]], h[[1, k]], h[[2, k]]];

    for target in targets {
        let restraints = target
            .molecule_restraints
            .iter()
            .chain(target.atom_restraints.iter().map(|(_, r)| r));
        for r in restraints {
            let Some(n) = r.plane_normal() else { continue };
            let n_norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if n_norm <= 0.0 {
                continue;
            }
            for (k, &periodic) in pbc.iter().enumerate() {
                if !periodic {
                    continue;
                }
                let a = lattice(k);
                let a_norm = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
                let cos = (n[0] * a[0] + n[1] * a[1] + n[2] * a[2]) / (n_norm * a_norm);
                if cos.abs() > 1e-9 {
                    return Err(PackError::PlaneAcrossPeriodicAxis { axis: k, normal: n });
                }
            }
        }
    }
    Ok(())
}

/// Do two lattices describe the same cell?
pub(crate) fn cells_agree(a: &SimBox, b: &SimBox) -> bool {
    let (ha, hb) = (a.h_view(), b.h_view());
    let (oa, ob) = (a.origin_view(), b.origin_view());
    let tol: F = 1e-9;
    (0..3).all(|i| (0..3).all(|j| (ha[[i, j]] - hb[[i, j]]).abs() <= tol))
        && (0..3).all(|k| (oa[k] - ob[k]).abs() <= tol)
        && a.pbc() == b.pbc()
}

/// Scan every restraint on every target for an `AtomRestraint::declared_cell`.
///
/// Mirrors [`derive_periodic_box`]: at most one distinct lattice may be
/// declared across the whole system, otherwise the packing has no
/// well-defined cell.
pub(crate) fn derive_cell(targets: &[Target]) -> Result<Option<CellDecl>, PackError> {
    let mut found: Option<CellDecl> = None;
    for target in targets {
        let restraints = target
            .molecule_restraints
            .iter()
            .chain(target.atom_restraints.iter().map(|(_, r)| r));
        for r in restraints {
            if let Some((h, origin, pbc)) = r.declared_cell() {
                let candidate = CellDecl::Matrix { h, origin, pbc };
                match found {
                    None => found = Some(candidate),
                    Some(existing) if existing == candidate => {}
                    Some(_) => {
                        return Err(PackError::InvalidCell {
                            detail: "restraints declare more than one lattice".to_string(),
                        });
                    }
                }
            }
        }
    }
    Ok(found)
}

/// Scan every restraint on every target for a `AtomRestraint::periodic_box`
/// override. Returns `Ok(None)` if no restraint declares one, `Ok(Some(...))`
/// if exactly one unique declaration exists (duplicates with identical
/// bounds + flags are allowed — they come from `with_global_restraint`
/// broadcast and from two targets sharing the same restraint object).
/// Returns `Err(ConflictingPeriodicBoxes)` when two declarations disagree
/// and `Err(InvalidPBCBox)` if the declared box has a non-positive extent
/// on any axis.
pub(crate) fn derive_periodic_box(targets: &[Target]) -> Result<Option<PeriodicSpec>, PackError> {
    let mut found: Option<PeriodicSpec> = None;
    for target in targets {
        let restraints = target
            .molecule_restraints
            .iter()
            .chain(target.atom_restraints.iter().map(|(_, r)| r));
        for r in restraints {
            if let Some(candidate) = r.periodic_box() {
                let (min, max, _periodic) = candidate;
                let length = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
                if length.iter().any(|&v| v <= 0.0) {
                    return Err(PackError::InvalidPBCBox { min, max });
                }
                match found {
                    None => found = Some(candidate),
                    Some(existing) if existing == candidate => {}
                    Some(existing) => {
                        return Err(PackError::ConflictingPeriodicBoxes {
                            first: existing,
                            second: candidate,
                        });
                    }
                }
            }
        }
    }
    Ok(found)
}

/// The resolved packing space: one cell for the whole run — a declared
/// lattice if there is one, else the periodic box, else nothing.
#[derive(Debug, Clone)]
pub(crate) struct ResolvedSpace {
    pub cell: Option<SimBox>,
    pub pbc: Option<PeriodicSpec>,
}

/// Resolve the caller's space declarations against the targets' own.
///
/// Verbatim extraction of the packer's stage-① block (engine-entry-split;
/// behavior-preserving, guarded by `examples_batch`): density → cubic box
/// (mass over ALL targets), box/cell exclusivity, restraint-derived
/// PBC/cell agreement, and the plane-across-periodic-axis rejection.
pub(crate) fn resolve_pack_space(
    density: Option<F>,
    periodic_box: Option<PeriodicSpec>,
    declared_cell: Option<CellDecl>,
    targets: &[Target],
) -> Result<ResolvedSpace, PackError> {
    let density_box: Option<PeriodicSpec> = if let Some(rho) = density {
        if periodic_box.is_some() || declared_cell.is_some() {
            return Err(PackError::DensityConflictsWithBox);
        }
        let mut total_amu: F = 0.0;
        for (i, t) in targets.iter().enumerate() {
            let per_copy = match t.mass {
                Some(m) => m,
                None => {
                    let mut m: F = 0.0;
                    for sym in &t.elements {
                        match Element::by_symbol(sym) {
                            Some(el) => m += el.atomic_mass() as F,
                            None => return Err(PackError::UnknownMass { target: i }),
                        }
                    }
                    m
                }
            };
            total_amu += per_copy * t.count.max(1) as F;
        }
        let l = (total_amu / (rho * 6.022_140_76e23) * 1e24).cbrt();
        Some(([0.0; 3], [l, l, l], [true; 3]))
    } else {
        None
    };
    let declared_box = periodic_box.or(density_box);

    if let Some((min, max, _)) = declared_box {
        let length = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
        if length.iter().any(|&v| v <= 0.0) {
            return Err(PackError::InvalidPBCBox { min, max });
        }
    }
    let derived = derive_periodic_box(targets)?;
    let derived_cell = derive_cell(targets)?;
    let declared_cell = match (declared_cell, derived_cell) {
        (None, other) | (other, None) => other,
        (Some(builder), Some(restraint)) => {
            let from_builder = builder.resolve()?;
            let from_restraint = restraint.resolve()?;
            if cells_agree(&from_builder, &from_restraint) {
                Some(builder)
            } else {
                return Err(PackError::InvalidCell {
                    detail: "the cell declared on the packer and the cell declared by a \
                             restraint describe different lattices"
                        .to_string(),
                });
            }
        }
    };
    if declared_cell.is_some() && (declared_box.is_some() || derived.is_some()) {
        return Err(PackError::InvalidCell {
            detail: "a declared cell and a periodic box are mutually exclusive; \
                     drop the `pbc` declaration or express it as the cell"
                .to_string(),
        });
    }
    let pbc = match (declared_box, derived) {
        (None, derived) => derived,
        (Some(global), None) => Some(global),
        (Some(global), Some(derived)) if global == derived => Some(global),
        (Some(global), Some(derived)) => {
            return Err(PackError::ConflictingPeriodicBoxes {
                first: global,
                second: derived,
            });
        }
    };

    let cell = match (declared_cell, pbc) {
        (Some(decl), _) => Some(decl.resolve()?),
        (None, Some((min, max, periodic))) => Some(
            SimBox::ortho(
                array![max[0] - min[0], max[1] - min[1], max[2] - min[2]],
                array![min[0], min[1], min[2]],
                periodic,
            )
            .map_err(|_| PackError::InvalidPBCBox { min, max })?,
        ),
        (None, None) => None,
    };
    if let Some(bx) = cell.as_ref() {
        reject_planes_across_periodic_axes(targets, bx)?;
    }
    Ok(ResolvedSpace { cell, pbc })
}

/// Broadcast global restraints to each target (scope equivalence, spec §4).
/// A no-op clone-free borrow when there are none.
pub(crate) fn broadcast_global_restraints<'a>(
    targets: &'a [Target],
    global: &[Arc<dyn AtomRestraint>],
) -> std::borrow::Cow<'a, [Target]> {
    if global.is_empty() {
        std::borrow::Cow::Borrowed(targets)
    } else {
        std::borrow::Cow::Owned(
            targets
                .iter()
                .map(|t| {
                    let mut t = t.clone();
                    for r in global {
                        t.molecule_restraints.push(Arc::clone(r));
                    }
                    t
                })
                .collect(),
        )
    }
}
