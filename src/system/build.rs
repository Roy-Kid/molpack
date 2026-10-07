//! Stage ② of the engine lifecycle: lower targets into a fully built
//! [`PackSystem`] (counts, per-copy conformers, radii, restraints, fixed
//! placements, AoS sync, frame constants).
//!
//! Every engine shares this one system construction; it is callable more
//! than once — chained engines build one system per stage.

use molrs::core::Element;
use molrs::op::F;

use crate::PackError;
use crate::euler::{compcart, eulerfixed};
use crate::system::PackSystem;
use crate::target::{CenteringMode, Target};

/// The three shared knobs system construction actually consumes.
#[derive(Debug, Clone)]
pub(crate) struct SystemKnobs {
    pub tolerance: F,
    pub short_tolerance: Option<(F, F)>,
    pub parallel_eval: bool,
}

/// Everything stage ② produces: the lowered system plus the counts the
/// later stages need. See [`build_system`].
pub(crate) struct BuiltSystem {
    pub(crate) sys: PackSystem,
    pub(crate) maxmove_per_type: Vec<usize>,
    pub(crate) ntype: usize,
    pub(crate) ntype_with_fixed: usize,
    pub(crate) ntotmol_free: usize,
    pub(crate) ntotat: usize,
    pub(crate) ntotat_free: usize,
}

/// One target's resolved per-atom properties — the per-type template Packmol
/// broadcasts to every copy (`app/packmol.f90` lines 503-515).
struct AtomPropsTemplate {
    radii: Vec<F>,
    fscale: Vec<F>,
    short_radii: Vec<F>,
    short_scale: Vec<F>,
    use_short: Vec<bool>,
}

impl AtomPropsTemplate {
    fn of(
        target: &Target,
        default_radius: F,
        global_short_radius: F,
        global_short_scale: F,
        short_on_globally: bool,
    ) -> Self {
        let opted_in = target.uses_short_radius();
        Self {
            radii: target.resolved_radii(default_radius),
            fscale: target.resolved_fscale(),
            short_radii: target.resolved_short_radii(global_short_radius),
            short_scale: target.resolved_short_radius_scale(global_short_scale),
            // A target opts in per atom; the packer's global switch turns it on
            // for everything else.
            use_short: opted_in
                .into_iter()
                .map(|opted| opted || short_on_globally)
                .collect(),
        }
    }

    /// The short penalty is meaningless unless it is the tighter of the two.
    fn validate(&self, itype: usize) -> Result<(), PackError> {
        for (iatom, &use_short) in self.use_short.iter().enumerate() {
            if use_short && self.short_radii[iatom] >= self.radii[iatom] {
                return Err(PackError::ShortRadiusNotShorter {
                    target: itype,
                    atom: iatom,
                    short_radius: self.short_radii[iatom],
                    radius: self.radii[iatom],
                });
            }
        }
        Ok(())
    }

    /// Write atom `iatom` of the template onto system slot `icart`.
    fn stamp(&self, sys: &mut PackSystem, icart: usize, iatom: usize) {
        sys.radius[icart] = self.radii[iatom];
        sys.radius_ini[icart] = self.radii[iatom];
        sys.fscale[icart] = self.fscale[iatom];
        sys.short_radius[icart] = self.short_radii[iatom];
        sys.short_radius_scale[icart] = self.short_scale[iatom];
        sys.use_short_radius[icart] = self.use_short[iatom];
    }
}

fn reference_coords(target: &Target) -> &[[F; 3]] {
    match target.centering {
        CenteringMode::Center => &target.ref_coords,
        CenteringMode::Off => &target.input_coords,
        CenteringMode::Auto => {
            if target.fixed_at.is_some() {
                &target.input_coords
            } else {
                &target.ref_coords
            }
        }
    }
}

/// Stage ② of `pack_with_report`: lower targets into a fully built
/// [`PackSystem`] (counts, per-copy conformers, radii, restraints,
/// fixed placements, AoS sync, frame constants). Callable more than once
/// — the mixed-method composition builds one system for growth and one
/// for the rigid stage (spec Design §6).
pub(crate) fn build_system(
    knobs: &SystemKnobs,
    targets: &[Target],
) -> Result<BuiltSystem, PackError> {
    // Split into free and fixed targets
    let free_targets: Vec<&Target> = targets.iter().filter(|t| t.fixed_at.is_none()).collect();
    let fixed_targets: Vec<&Target> = targets.iter().filter(|t| t.fixed_at.is_some()).collect();

    let ntype = free_targets.len();
    let ntype_with_fixed = ntype + fixed_targets.len();

    // Count atoms
    let ntotmol_free: usize = free_targets.iter().map(|t| t.count).sum();
    let ntotat_free: usize = free_targets.iter().map(|t| t.count * t.natoms()).sum();
    let ntotat_fixed: usize = fixed_targets.iter().map(|t| t.natoms()).sum();
    let ntotat = ntotat_free + ntotat_fixed;

    // Build PackSystem
    let mut sys = PackSystem::new(ntotat, ntotmol_free, ntype);
    sys.ntype_with_fixed = ntype_with_fixed;
    sys.nfixedat = ntotat_fixed;
    // comptype is initialized with size ntype; resize to include fixed types
    sys.comptype = vec![true; ntype_with_fixed];

    // Fill nmols, natoms, idfirst for free types
    let mut cum_atoms = 0usize;
    let mut coor = Vec::new();
    let mut maxmove_per_type = vec![0usize; ntype];

    sys.nmols = vec![0; ntype_with_fixed];
    sys.natoms = vec![0; ntype_with_fixed];
    sys.idfirst = vec![0; ntype_with_fixed];
    sys.constrain_rot = vec![[false; 3]; ntype];
    sys.rot_bound = vec![[[0.0; 2]; 3]; ntype];

    for (itype, target) in free_targets.iter().enumerate() {
        sys.nmols[itype] = target.count;
        sys.natoms[itype] = target.natoms();
        sys.idfirst[itype] = cum_atoms;
        // One reference conformer **per copy**: in-loop optimizers relax
        // each copy independently, so copies must not share a block.
        // Layout matches `xcart` exactly (type-major, copy-major,
        // atom-minor), so a single index addresses both buffers.
        for _ in 0..target.count {
            coor.extend_from_slice(reference_coords(target));
        }
        cum_atoms += target.natoms() * target.count;

        maxmove_per_type[itype] = target.perturb_budget.unwrap_or(target.count);
        for k in 0..3 {
            if let Some((center, half_width)) = target.rotation_bound[k] {
                sys.constrain_rot[itype][k] = true;
                sys.rot_bound[itype][k][0] = center.radians();
                sys.rot_bound[itype][k][1] = half_width.radians();
            }
        }
    }

    for (fi, target) in fixed_targets.iter().enumerate() {
        let itype = ntype + fi;
        sys.nmols[itype] = 1;
        sys.natoms[itype] = target.natoms();
        sys.idfirst[itype] = cum_atoms;
        coor.extend_from_slice(reference_coords(target));
        cum_atoms += target.natoms();
    }
    sys.coor = coor;

    // Assign radii, element symbols, and per-atom (itype, imol) tags.
    //
    // Radii follow Packmol's layering (packmol.f90 lines 281-515): every
    // atom starts at `tolerance / 2` (line 283), a structure-level
    // `radius` covers the whole species, and an atom-specific `radius`
    // overrides selected atoms. `Target::resolved_radii` collapses those
    // layers into one per-type template, which is then broadcast to every
    // copy — Packmol does the same broadcast explicitly at lines 503-515.
    // VdW radii from the source file are deliberately not used.
    //
    // `ibtype` / `ibmol` are derivable from position in the sequential
    // atom layout, so we set them here once instead of having
    // `insert_atom_in_cell` rewrite the same constants on every eval.
    let default_radius = knobs.tolerance / 2.0;
    let (global_short_radius, global_short_scale) =
        knobs.short_tolerance.unwrap_or((default_radius / 2.0, 3.0));
    let short_on_globally = knobs.short_tolerance.is_some();
    let mut icart = 0usize;
    for (itype, target) in free_targets.iter().enumerate() {
        let props = AtomPropsTemplate::of(
            target,
            default_radius,
            global_short_radius,
            global_short_scale,
            short_on_globally,
        );
        props.validate(itype)?;
        for imol in 0..target.count {
            for iatom in 0..target.natoms() {
                props.stamp(&mut sys, icart, iatom);
                sys.ibtype[icart] = itype;
                sys.ibmol[icart] = imol;
                sys.elements[icart] = Element::by_symbol(&target.elements[iatom]);
                icart += 1;
            }
        }
    }
    for (fi, target) in fixed_targets.iter().enumerate() {
        let itype = ntype + fi;
        let props = AtomPropsTemplate::of(
            target,
            default_radius,
            global_short_radius,
            global_short_scale,
            short_on_globally,
        );
        props.validate(itype)?;
        for iatom in 0..target.natoms() {
            props.stamp(&mut sys, icart, iatom);
            sys.ibtype[icart] = itype;
            sys.ibmol[icart] = 0;
            sys.elements[icart] = Element::by_symbol(&target.elements[iatom]);
            icart += 1;
        }
    }

    // Assign restraints: per-atom
    let mut irest_pool = Vec::new();
    let mut iratom_lists: Vec<Vec<usize>> = vec![Vec::new(); ntotat];
    let mut icart = 0usize;
    for target in free_targets.iter() {
        for _imol in 0..target.count {
            for iatom in 0..target.natoms() {
                // molecule-level restraints applied to all atoms
                for r in &target.molecule_restraints {
                    let irest = irest_pool.len();
                    irest_pool.push(std::sync::Arc::clone(r));
                    iratom_lists[icart].push(irest);
                }
                // atom-subset restraints
                for (indices, restraint) in &target.atom_restraints {
                    if indices.contains(&iatom) {
                        let irest = irest_pool.len();
                        irest_pool.push(std::sync::Arc::clone(restraint));
                        iratom_lists[icart].push(irest);
                    }
                }
                icart += 1;
            }
        }
    }
    // Fixed atoms: no restraints needed (they are placed directly)
    sys.restraints = irest_pool;
    sys.iratom_offsets.clear();
    sys.iratom_offsets.reserve(ntotat + 1);
    sys.iratom_offsets.push(0);
    for atom_restraints in &iratom_lists {
        let next = sys.iratom_offsets.last().copied().unwrap_or(0) + atom_restraints.len();
        sys.iratom_offsets.push(next);
    }
    sys.iratom_indices.clear();
    sys.iratom_indices
        .reserve(sys.iratom_offsets.last().copied().unwrap_or(0));
    for atom_restraints in iratom_lists {
        sys.iratom_indices.extend(atom_restraints);
    }

    // Group-level (collective) restraints: one entry per (free type, restraint).
    // The free target's position is its 0-based type index `itype`.
    sys.collective.clear();
    for (itype, target) in free_targets.iter().enumerate() {
        for r in &target.collective_restraints {
            sys.collective.push((itype, std::sync::Arc::clone(r)));
        }
    }

    // Handle fixed molecules: place them using eulerfixed
    let free_atoms = ntotat_free;
    let mut fixed_icart = free_atoms;
    for target in fixed_targets.iter() {
        let fp = target.fixed_at.as_ref().unwrap();
        let (v1, v2, v3) = eulerfixed(
            fp.orientation[0].radians(),
            fp.orientation[1].radians(),
            fp.orientation[2].radians(),
        );
        let ref_coords = reference_coords(target);
        for ref_coord in ref_coords.iter().take(target.natoms()) {
            let pos = compcart(&fp.position, ref_coord, &v1, &v2, &v3);
            sys.xcart[fixed_icart] = pos;
            sys.fixedatom[fixed_icart] = true;
            fixed_icart += 1;
        }
    }
    // Populate the AoS `atom_props` mirror from the individual per-atom
    // Vecs now that every hot-loop field is finalized. `sync_atom_props`
    // also refreshes the `any_fixed_atoms` / `any_short_radius`
    // summary flags used by the hot-loop fast paths.
    sys.sync_atom_props();
    // Plumb the builder's opt-in parallel flag through to the
    // objective kernels.
    sys.parallel_pair_eval = knobs.parallel_eval;

    Ok(BuiltSystem {
        sys,
        maxmove_per_type,
        ntype,
        ntype_with_fixed,
        ntotmol_free,
        ntotat,
        ntotat_free,
    })
}

#[cfg(test)]
mod short_radius_tests {
    //! The short penalty is only meaningful as the tighter of the two radii,
    //! so system construction refuses a short radius that is not shorter —
    //! by target and atom index, never silently.

    use super::{SystemKnobs, build_system};
    use crate::PackError;
    use crate::Target;
    use molrs::op::F;

    fn knobs() -> SystemKnobs {
        SystemKnobs {
            tolerance: 4.0,
            short_tolerance: None,
            parallel_eval: false,
        }
    }

    fn two_atoms(count: usize) -> Target {
        Target::from_coords(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &[1.0, 1.0], count)
    }

    #[test]
    fn a_short_radius_below_the_radius_is_accepted() {
        let target = two_atoms(1).with_radius(2.0).with_short_radius(1.0);
        assert!(build_system(&knobs(), &[target]).is_ok());
    }

    #[test]
    fn a_short_radius_equal_to_the_radius_names_the_atom() {
        let target = two_atoms(1)
            .with_radius(2.0)
            .with_atom_short_radius(&[1], 2.0);

        let err = match build_system(&knobs(), &[target]) {
            Err(e) => e,
            Ok(_) => panic!("a short radius equal to the radius must be refused"),
        };

        match err {
            PackError::ShortRadiusNotShorter {
                target,
                atom,
                short_radius,
                radius,
            } => {
                assert_eq!((target, atom), (0, 1), "the error names target and atom");
                assert_eq!((short_radius, radius), (2.0 as F, 2.0 as F));
            }
            other => panic!("expected ShortRadiusNotShorter, got {other:?}"),
        }
    }

    #[test]
    fn an_atom_that_never_opted_in_is_not_checked() {
        // Default short radius (half the tolerance) exceeds this radius, but
        // the atom never asked for the short penalty, so nothing is refused.
        let target = two_atoms(1).with_radius(0.1);
        assert!(build_system(&knobs(), &[target]).is_ok());
    }
}
