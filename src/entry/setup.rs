//! Shared pack-space setup: density / periodic-box / cell resolution and
//! global-restraint broadcast.
//!
//! Every entry needs the same answers before any algorithm runs: what box
//! does the system live in, and which restraints apply to every target.
//! This machinery moved out of the packer (engine-entry-split) so the entry
//! lifecycle owns it once — never re-declared per algorithm.

use std::sync::Arc;

use molrs::core::Element;
use molrs::core::SimBox;
use molrs::core::constants::{ANGSTROM3_PER_CM3, AVOGADRO};
use molrs::op::F;
use ndarray::array;

use crate::AtomRestraint;
use crate::PackError;
use crate::Target;
use crate::restraint::cell::{simbox_from_lengths_angles, simbox_from_matrix};

pub(crate) type PeriodicSpec = ([F; 3], [F; 3], [bool; 3]);

/// A packing cell as the caller declared it, resolved to a [`SimBox`]
/// inside the engine lifecycle so the builders can stay infallible.
/// `Resolved` is a cell that is already a [`SimBox`] (a seed's inherited
/// cell), carried as-is rather than taken apart and rebuilt (boxed: a
/// `SimBox` is several times the size of the other variants).
#[derive(Clone, Debug)]
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
    Resolved(Box<SimBox>),
}

impl CellDecl {
    pub(crate) fn resolve(self) -> Result<SimBox, PackError> {
        match self {
            CellDecl::LengthsAngles {
                lengths,
                angles_deg,
                pbc,
            } => simbox_from_lengths_angles(lengths, angles_deg, pbc),
            CellDecl::Matrix { h, origin, pbc } => simbox_from_matrix(h, origin, pbc),
            CellDecl::Resolved(bx) => Ok(*bx),
        }
    }
}

/// Reject a restraint that is open along a periodic lattice direction.
///
/// Restraints are evaluated at lab coordinates, never at wrapped ones, so a
/// restraint used with periodic boundaries has to answer the same way in every
/// image. [`AtomRestraint::holds_along`] is where each restraint says
/// whether it does — by confining atoms to one image, or by repeating with the
/// lattice. One that does neither is satisfied or not depending on which image
/// an atom drifted into, and on where the user happened to put the origin.
/// There is no interpretation to salvage, so the declaration is refused,
/// naming the axis and the restraint.
///
/// This covers every restraint the same way: the `.inp` plane kernels, a
/// molrs region lifted by [`RegionRestraint`](crate::RegionRestraint) — a
/// half-space across a periodic axis is refused whichever spelling built it —
/// and a user-supplied `f` / `fg`, which is taken at its word by the trait's
/// default.
pub(crate) fn reject_restraints_across_periodic_axes(
    targets: &[Target],
    bx: &SimBox,
) -> Result<(), PackError> {
    let pbc = bx.pbc();
    if !pbc.iter().any(|&p| p) {
        return Ok(());
    }
    let lattice = |k: usize| {
        let v = bx.lattice(k);
        [v[0], v[1], v[2]]
    };

    for target in targets {
        let restraints = target
            .molecule_restraints
            .iter()
            .chain(target.atom_restraints.iter().map(|(_, r)| r));
        for r in restraints {
            for (k, &periodic) in pbc.iter().enumerate() {
                if periodic && !r.holds_along(lattice(k)) {
                    return Err(PackError::RestraintAcrossPeriodicAxis {
                        axis: k,
                        restraint: r.name(),
                    });
                }
            }
        }
    }
    Ok(())
}

/// Absolute tolerance (Å) under which two declared lattices are the same cell,
/// per [`SimBox::approx_eq`]: every matrix and origin entry within it, identical
/// PBC flags, and the same `is_cell_defined`. The last one is stricter than the
/// hand-rolled check it replaced, on purpose: a no-cell box never agrees with
/// a real lattice.
const CELL_AGREE_TOL: F = 1e-9;

/// Scan every restraint on every target for an `AtomRestraint::declared_cell`.
///
/// At most one distinct lattice may be declared across the whole system,
/// otherwise the packing has no well-defined cell.
pub(crate) fn derive_cell(targets: &[Target]) -> Result<Option<SimBox>, PackError> {
    let mut found: Option<SimBox> = None;
    for target in targets {
        let restraints = target
            .molecule_restraints
            .iter()
            .chain(target.atom_restraints.iter().map(|(_, r)| r));
        for r in restraints {
            if let Some(candidate) = r.declared_cell() {
                match &found {
                    None => found = Some(candidate),
                    Some(existing) if existing.approx_eq(&candidate, CELL_AGREE_TOL) => {}
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

/// The resolved packing space: one cell for the whole run — a declared
/// lattice if there is one, else the periodic box, else nothing.
#[derive(Debug, Clone)]
pub(crate) struct ResolvedSpace {
    pub cell: Option<SimBox>,
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
        // rho in g/cm³, the mass in g/mol: the volume in cm³, then in Å³.
        let l = (total_amu / (rho * AVOGADRO) * ANGSTROM3_PER_CM3).cbrt();
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
    let derived_cell = derive_cell(targets)?;
    let declared_cell = match (
        declared_cell.map(CellDecl::resolve).transpose()?,
        derived_cell,
    ) {
        (None, other) | (other, None) => other,
        (Some(builder), Some(restraint)) => {
            if builder.approx_eq(&restraint, CELL_AGREE_TOL) {
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
    if declared_cell.is_some() && declared_box.is_some() {
        return Err(PackError::InvalidCell {
            detail: "a declared cell and a periodic box are mutually exclusive; \
                     drop the `pbc` declaration or express it as the cell"
                .to_string(),
        });
    }
    let cell = match (declared_cell, declared_box) {
        (Some(cell), _) => Some(cell),
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
        reject_restraints_across_periodic_axes(targets, bx)?;
    }
    Ok(ResolvedSpace { cell })
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

#[cfg(test)]
mod periodic_declaration_tests {
    //! Periodicity is declared on the engine entry and nowhere else. What a
    //! restraint still answers for is whether it survives the wrap, and
    //! [`reject_restraints_across_periodic_axes`] refuses the ones that do
    //! not — for the `.inp` plane kernels here, and for a lifted molrs region
    //! in `region_under_wrap_tests` below.

    use crate::restraint::geometric::{AbovePlaneRestraint, BelowPlaneRestraint};
    use crate::{GenCanPack, PackEngine, PackError, Target};

    fn one_atom(n: usize) -> Target {
        Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], n)
    }

    /// The `.inp` plane kernels declare a normal, and a normal with a
    /// component along a periodic lattice vector has no well-defined side:
    /// [`reject_restraints_across_periodic_axes`] names the first such axis.
    #[test]
    fn a_plane_across_a_periodic_axis_is_rejected_naming_the_axis() {
        let target = one_atom(4).with_restraint(AbovePlaneRestraint::new([0.0, 0.0, 1.0], 5.0));
        let err = GenCanPack::new()
            .with_seed(1)
            .with_periodic_box([0.0; 3], [30.0; 3], [true; 3])
            .run(std::slice::from_ref(&target), 2)
            .expect_err("plane across a periodic axis must be rejected");
        assert!(
            matches!(err, PackError::RestraintAcrossPeriodicAxis { axis: 2, .. }),
            "unexpected error: {err:?}"
        );
        assert!(
            format!("{err}").contains("periodic lattice direction 2"),
            "the error must name the offending axis, got: {err}"
        );

        let target = one_atom(4).with_restraint(BelowPlaneRestraint::new([1.0, 0.0, 0.0], 5.0));
        let err = GenCanPack::new()
            .with_seed(1)
            .with_periodic_box([0.0; 3], [30.0; 3], [true, true, false])
            .run(std::slice::from_ref(&target), 2)
            .expect_err("plane across a periodic axis must be rejected");
        assert!(
            matches!(err, PackError::RestraintAcrossPeriodicAxis { axis: 0, .. }),
            "unexpected error: {err:?}"
        );
    }

    /// The same plane along a confined axis is the slab an interface needs.
    #[test]
    fn a_plane_along_a_confined_axis_is_accepted() {
        let target = one_atom(4).with_restraint(AbovePlaneRestraint::new([0.0, 0.0, 1.0], 5.0));
        let result = GenCanPack::new()
            .with_seed(1)
            .with_periodic_box([0.0; 3], [30.0; 3], [true, true, false])
            .run(std::slice::from_ref(&target), 5);
        assert!(
            result.is_ok(),
            "expected the slab plane to be accepted: {result:?}"
        );
    }

    /// A periodic box with no extent on an axis is named, not silently
    /// collapsed into a two-dimensional cell.
    #[test]
    fn zero_extent_declaration_is_rejected() {
        let result = GenCanPack::new()
            .with_seed(7)
            .with_periodic_box([0.0; 3], [10.0, 0.0, 10.0], [true; 3])
            .run(&[one_atom(1)], 5);
        assert!(
            matches!(result, Err(PackError::InvalidPBCBox { .. })),
            "expected InvalidPBCBox, got: {result:?}"
        );
    }
}

#[cfg(test)]
mod region_under_wrap_tests {
    //! The same rule, reached through the public spelling: a molrs region
    //! lifted by `RegionRestraint`. This is what used to go unchecked — a
    //! half-space handed in as a region was accepted whatever the cell did,
    //! and its "inside" then depended on where the origin sat.

    use std::sync::Arc;

    use crate::AtomRestraint;
    use crate::{GenCanPack, PackEngine, PackError, RegionRestraint, Target};
    use molrs::core::{AndRegion, Cuboid, HalfSpace, NotRegion, Region, Sphere};
    use ndarray::array;

    fn one_atom(n: usize) -> Target {
        Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], n)
    }

    fn lift(r: impl Region + 'static) -> RegionRestraint {
        RegionRestraint(Arc::new(r))
    }

    fn half_space_above_z() -> RegionRestraint {
        // Inside is `z <= 5`; open downwards, so it repeats along x and y
        // and not along z.
        lift(HalfSpace::new([0.0, 0.0, 1.0], [0.0, 0.0, 5.0]).expect("plane"))
    }

    #[test]
    fn a_lifted_half_space_across_a_periodic_axis_is_rejected() {
        let target = one_atom(4).with_restraint(half_space_above_z());
        let err = GenCanPack::new()
            .with_seed(1)
            .with_periodic_box([0.0; 3], [30.0; 3], [true; 3])
            .run(&[target], 2)
            .expect_err("a region open along z cannot be used with periodic z");

        match err {
            PackError::RestraintAcrossPeriodicAxis { axis, restraint } => {
                assert_eq!(axis, 2, "the error names the offending axis");
                assert_eq!(restraint, "RegionRestraint", "and the restraint");
            }
            other => panic!("expected RestraintAcrossPeriodicAxis, got {other:?}"),
        }
    }

    #[test]
    fn the_same_region_along_confined_axes_is_the_slab_it_should_be() {
        let target = one_atom(4).with_restraint(half_space_above_z());
        let result = GenCanPack::new()
            .with_seed(1)
            .with_periodic_box([0.0; 3], [30.0; 3], [true, true, false])
            .run(&[target], 5);
        assert!(
            result.is_ok(),
            "the region repeats along x and y, which is exactly the slab an \
             interface system asks for: {result:?}"
        );
    }

    /// The decision is per restraint, not per shape family: a bounded region
    /// keeps its atoms inside one image, so every axis may wrap.
    #[test]
    fn a_bounded_region_is_usable_along_every_axis() {
        let cuboid = lift(Cuboid::new(array![0.0, 0.0, 0.0], array![30.0, 30.0, 30.0]));
        let ball = lift(Sphere::new(array![15.0, 15.0, 15.0], 10.0));
        for shift in [[30.0, 0.0, 0.0], [0.0, 30.0, 0.0], [0.0, 0.0, 30.0]] {
            assert!(cuboid.holds_along(shift), "cuboid, {shift:?}");
            assert!(ball.holds_along(shift), "sphere, {shift:?}");
        }
    }

    /// A confining member rescues an open one: the intersection is bounded,
    /// so the composition is usable even though the half-space alone is not.
    #[test]
    fn a_confined_composition_is_usable_where_its_open_member_is_not() {
        let open = HalfSpace::new([0.0, 0.0, 1.0], [0.0, 0.0, 5.0]).expect("plane");
        let cell = Cuboid::new(array![0.0, 0.0, 0.0], array![30.0, 30.0, 30.0]);
        let both = lift(AndRegion::new(Arc::new(cell), Arc::new(open)));
        assert!(
            both.holds_along([0.0, 0.0, 30.0]),
            "the box confines z, so the plane inside it has one meaning",
        );
    }

    /// A void — pack everything *outside* a ball — keeps working: the region
    /// molpack must not silently mishandle is the open one, not this.
    #[test]
    fn a_void_region_stays_usable() {
        let void = lift(NotRegion::new(Arc::new(Sphere::new(
            array![15.0, 15.0, 15.0],
            5.0,
        ))));
        assert!(void.holds_along([0.0, 0.0, 30.0]));
    }
}

#[cfg(test)]
mod broadcast_tests {
    //! The scope equivalence a global restraint promises: declaring it once on
    //! the engine is the same thing as attaching it to every target by hand.

    use std::sync::Arc;

    use super::broadcast_global_restraints;
    use crate::AtomRestraint;
    use crate::Target;
    use crate::restraint::geometric::InsideCubeRestraint;

    fn one_atom(name: &str) -> Target {
        Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 2).with_name(name)
    }

    #[test]
    fn no_global_restraint_borrows_the_targets_untouched() {
        let targets = vec![one_atom("a"), one_atom("b")];
        let out = broadcast_global_restraints(&targets, &[]);
        assert!(matches!(out, std::borrow::Cow::Borrowed(_)));
        assert!(out.iter().all(|t| t.molecule_restraints.is_empty()));
    }

    #[test]
    fn a_global_restraint_reaches_every_target() {
        let targets = vec![one_atom("a"), one_atom("b")];
        let global: Vec<Arc<dyn AtomRestraint>> =
            vec![Arc::new(InsideCubeRestraint::new([0.0; 3], 30.0))];

        let out = broadcast_global_restraints(&targets, &global);

        assert_eq!(out.len(), 2, "broadcasting must not add or drop a target");
        for t in out.iter() {
            assert_eq!(
                t.molecule_restraints.len(),
                1,
                "every target carries the global restraint"
            );
        }
    }

    #[test]
    fn a_global_restraint_is_appended_to_the_target_s_own() {
        let own = InsideCubeRestraint::new([0.0; 3], 10.0);
        let targets = vec![one_atom("a").with_restraint(own)];
        let global: Vec<Arc<dyn AtomRestraint>> =
            vec![Arc::new(InsideCubeRestraint::new([0.0; 3], 30.0))];

        let out = broadcast_global_restraints(&targets, &global);

        assert_eq!(
            out[0].molecule_restraints.len(),
            2,
            "the global restraint is added to the target's own, never replacing it"
        );
    }
}
