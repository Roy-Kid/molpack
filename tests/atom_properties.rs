//! Per-atom packing properties: `radius`, `fscale`, `short_radius`,
//! `short_radius_scale`.
//!
//! All four follow the same Packmol scheme (`app/packmol.f90` lines 281-515):
//! a global default, a structure-level keyword covering every atom of every
//! copy, and an atom-specific keyword inside an `atoms ... end atoms` block
//! that overrides the selected atoms. Per-atom values are a **per-type
//! template** — every copy of a type gets the same ones.

use molpack::{F, GenCanPack, InsideBoxRestraint, PackEngine, Target};

fn coords(n: usize) -> Vec<[F; 3]> {
    (0..n).map(|i| [i as F * 3.0, 0.0, 0.0]).collect()
}

fn target(n: usize, count: usize) -> Target {
    Target::from_coords(&coords(n), &vec![1.5; n], count)
}

// ── layer 1: the global default ────────────────────────────────────────────

#[test]
fn radius_defaults_to_the_packer_default() {
    let t = target(3, 2);
    assert_eq!(t.resolved_radii(2.0), vec![2.0; 3]);
}

/// The frame's van der Waals radii are *not* the packing radii — Packmol packs
/// on `tolerance / 2` unless told otherwise.
#[test]
fn frame_vdw_radii_do_not_become_packing_radii() {
    let t = Target::from_coords(&coords(3), &[9.0, 9.0, 9.0], 1);
    assert_eq!(t.resolved_radii(2.0), vec![2.0; 3]);
}

// ── layer 2: structure-level ───────────────────────────────────────────────

#[test]
fn structure_radius_covers_every_atom() {
    let t = target(3, 2).with_radius(3.5);
    assert_eq!(t.resolved_radii(2.0), vec![3.5; 3]);
}

// ── layer 3: atom-specific ─────────────────────────────────────────────────

#[test]
fn atom_radius_applies_only_to_the_selected_atoms() {
    let t = target(3, 1).with_atom_radius(&[1], 5.0);
    assert_eq!(t.resolved_radii(2.0), vec![2.0, 5.0, 2.0]);
}

/// Packmol applies the structure-level pass first and lets the atom-specific
/// pass overwrite it.
#[test]
fn atom_radius_overrides_the_structure_radius() {
    let t = target(3, 1).with_radius(3.0).with_atom_radius(&[2], 6.0);
    assert_eq!(t.resolved_radii(2.0), vec![3.0, 3.0, 6.0]);
}

/// A structure-level radius set *after* an atom-specific one still wins for
/// every atom — the builder applies calls in order, so the caller controls
/// precedence explicitly rather than by a hidden rule.
#[test]
fn a_later_structure_radius_replaces_earlier_atom_radii() {
    let t = target(3, 1).with_atom_radius(&[0], 9.0).with_radius(3.0);
    assert_eq!(t.resolved_radii(2.0), vec![3.0; 3]);
}

#[test]
fn repeated_atom_radius_calls_accumulate() {
    let t = target(4, 1)
        .with_atom_radius(&[0, 1], 5.0)
        .with_atom_radius(&[3], 7.0);
    assert_eq!(t.resolved_radii(2.0), vec![5.0, 5.0, 2.0, 7.0]);
}

#[test]
#[should_panic(expected = "atom index")]
fn atom_radius_rejects_an_out_of_range_index() {
    let _ = target(3, 1).with_atom_radius(&[3], 5.0);
}

#[test]
#[should_panic(expected = "positive")]
fn radius_rejects_a_non_positive_value() {
    let _ = target(3, 1).with_radius(0.0);
}

// ── the packer honours them ────────────────────────────────────────────────

/// A bulky atom must hold its neighbours further away than a default one.
/// Without per-atom radii this is unexpressible, which is the whole point.
#[test]
fn a_bulky_atom_enforces_a_larger_separation() {
    let cell = InsideBoxRestraint::new([0.0; 3], [30.0, 30.0, 30.0], [false; 3]);
    let bulky = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], 1)
        .with_name("bulky")
        .with_radius(6.0)
        .with_restraint(cell);
    let small = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], 20)
        .with_name("small")
        .with_restraint(cell);

    let result = GenCanPack::new()
        .with_tolerance(2.0)
        .with_seed(5)
        .run(&[bulky, small], 80)
        .expect("pack");

    let pos = result.positions();
    let (b, rest) = pos.split_at(1);
    let closest = rest
        .iter()
        .map(|p| {
            ((p[0] - b[0][0]).powi(2) + (p[1] - b[0][1]).powi(2) + (p[2] - b[0][2]).powi(2)).sqrt()
        })
        .fold(F::MAX, F::min);

    // bulky radius 6.0 + small radius (tolerance/2 = 1.0) = 7.0
    assert!(
        closest > 6.5,
        "a radius-6 atom must keep radius-1 atoms ~7 A away, closest was {closest:.3}"
    );
}

/// Per-atom radii are a per-type template: every copy is packed with the same
/// ones, so a two-atom molecule with one bulky end keeps that end apart in
/// every copy.
#[test]
fn per_atom_radii_apply_to_every_copy() {
    let cell = InsideBoxRestraint::new([0.0; 3], [40.0, 40.0, 40.0], [false; 3]);
    let dumbbell = Target::from_coords(&[[0.0, 0.0, 0.0], [8.0, 0.0, 0.0]], &[1.5; 2], 6)
        .with_name("dumbbell")
        .with_atom_radius(&[0], 5.0)
        .with_restraint(cell);

    let result = GenCanPack::new()
        .with_tolerance(2.0)
        .with_seed(9)
        .run(&[dumbbell], 80)
        .expect("pack");

    let pos = result.positions();
    // Bulky atoms are index 0 of each copy; two of them need 10 A between
    // centres, while the default atoms need only 2 A.
    let mut closest_bulky = F::MAX;
    for i in 0..6 {
        for j in (i + 1)..6 {
            let (a, b) = (pos[i * 2], pos[j * 2]);
            let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
            closest_bulky = closest_bulky.min(d);
        }
    }
    assert!(
        closest_bulky > 9.0,
        "bulky-bulky separation must hold in every copy, closest was {closest_bulky:.3}"
    );
}

// ── Packmol `.inp` parity ──────────────────────────────────────────────────

mod script {
    use molpack::script::parse;

    fn plan(src: &str) -> molpack::script::ScriptPlan {
        parse(src)
            .expect("parse")
            .lower(std::path::Path::new("."))
            .expect("lower")
    }

    /// Structure-level `radius`, as Packmol reads it outside an `atoms` block.
    #[test]
    fn structure_radius_is_parsed() {
        let p = plan(
            "tolerance 2.0\noutput o.xyz\n\
             structure a.pdb\n  number 3\n  radius 3.5\n\
             inside box 0. 0. 0. 10. 10. 10.\nend structure\n",
        );
        assert_eq!(p.structures[0].radius, Some(3.5));
    }

    /// Atom-specific `radius`, inside an `atoms ... end atoms` block.
    #[test]
    fn atom_group_radius_is_parsed() {
        let p = plan(
            "tolerance 2.0\noutput o.xyz\n\
             structure a.pdb\n  number 1\n\
             inside box 0. 0. 0. 10. 10. 10.\n\
             atoms 1 3\n  radius 6.0\nend atoms\n\
             end structure\n",
        );
        let g = &p.structures[0].atom_groups[0];
        assert_eq!(g.atom_indices, vec![1, 3]);
        assert_eq!(g.radius, Some(6.0));
    }

    /// Lowering must reproduce the Rust API's layering, with the script's
    /// 1-based indices mapped to 0-based.
    #[test]
    fn lowering_applies_both_radius_layers() {
        let p = plan(
            "tolerance 4.0\noutput o.xyz\n\
             structure a.pdb\n  number 2\n  radius 3.0\n\
             inside box 0. 0. 0. 10. 10. 10.\n\
             atoms 3\n  radius 6.0\nend atoms\n\
             end structure\n",
        );
        let t = p.structures[0].apply(molpack::Target::from_coords(
            &[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            &[1.5; 3],
            2,
        ));
        assert_eq!(t.resolved_radii(2.0), vec![3.0, 3.0, 6.0]);
    }

    #[test]
    fn a_non_positive_radius_is_a_parse_error() {
        let src = "tolerance 2.0\noutput o.xyz\n\
                   structure a.pdb\n  number 1\n  radius -1.0\n\
                   inside box 0. 0. 0. 10. 10. 10.\nend structure\n";
        assert!(parse(src).is_err(), "negative radius must be rejected");
    }
}

// ── fscale ─────────────────────────────────────────────────────────────────
//
// Packmol's per-atom weight on the overlap penalty: the pair term is scaled by
// `fscale_i * fscale_j` (`fparc.f90`), so a species can be made "softer" than
// the rest without changing how far apart it is asked to sit.

#[test]
fn fscale_defaults_to_one() {
    assert_eq!(target(3, 1).resolved_fscale(), vec![1.0; 3]);
}

#[test]
fn structure_fscale_covers_every_atom() {
    let t = target(3, 1).with_fscale(0.25);
    assert_eq!(t.resolved_fscale(), vec![0.25; 3]);
}

#[test]
fn atom_fscale_overrides_the_structure_value() {
    let t = target(3, 1).with_fscale(0.5).with_atom_fscale(&[1], 2.0);
    assert_eq!(t.resolved_fscale(), vec![0.5, 2.0, 0.5]);
}

#[test]
#[should_panic(expected = "positive")]
fn fscale_rejects_a_non_positive_value() {
    let _ = target(3, 1).with_fscale(0.0);
}

/// A species with a small `fscale` is penalised less for the same overlap, so
/// when the box is too crowded to satisfy everyone the optimizer spends its
/// budget on the full-weight species and lets the soft one crowd together.
///
/// Note the two species: a *uniform* fscale on a single species only rescales
/// the whole objective, which leaves the minimiser untouched. The weight is
/// only observable as a relative preference.
#[test]
fn a_soft_species_absorbs_the_crowding() {
    let cell = InsideBoxRestraint::new([0.0; 3], [14.0, 14.0, 14.0], [false; 3]);
    let soft = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], 20)
        .with_name("soft")
        .with_fscale(0.02)
        .with_restraint(cell);
    let firm = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], 20)
        .with_name("firm")
        .with_restraint(cell);

    let result = GenCanPack::new()
        .with_tolerance(4.0)
        .with_seed(3)
        .run(&[soft, firm], 60)
        .expect("pack");

    let pos = result.positions();
    let closest_within = |slice: &[[F; 3]]| {
        let mut m = F::MAX;
        for i in 0..slice.len() {
            for j in (i + 1)..slice.len() {
                let (a, b) = (slice[i], slice[j]);
                let d =
                    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
                m = m.min(d);
            }
        }
        m
    };
    let soft_min = closest_within(&pos[..20]);
    let firm_min = closest_within(&pos[20..]);
    assert!(
        soft_min < firm_min,
        "the soft species should absorb the crowding: soft {soft_min:.3} vs firm {firm_min:.3}"
    );
}

// ── short radius ───────────────────────────────────────────────────────────
//
// Packmol's optional second, shorter-range penalty (`use_short_tol` /
// `short_tol_dist` / `short_tol_scale`, overridable per structure and per
// atom). It lets a pair approach past the main radius while still being
// stopped hard at a smaller one.

#[test]
fn short_radius_is_off_until_asked_for() {
    let t = target(3, 1);
    assert!(!t.uses_short_radius().iter().any(|&b| b));
}

#[test]
fn structure_short_radius_enables_it_for_every_atom() {
    let t = target(3, 1).with_short_radius(0.5);
    assert_eq!(t.resolved_short_radii(1.0), vec![0.5; 3]);
    assert_eq!(t.uses_short_radius(), vec![true; 3]);
}

#[test]
fn atom_short_radius_enables_only_the_selected_atoms() {
    let t = target(3, 1).with_atom_short_radius(&[2], 0.5);
    assert_eq!(t.uses_short_radius(), vec![false, false, true]);
}

/// Setting only the scale still opts the atom in, as Packmol does.
#[test]
fn short_radius_scale_alone_enables_the_penalty() {
    let t = target(2, 1).with_atom_short_radius_scale(&[0], 5.0);
    assert_eq!(t.uses_short_radius(), vec![true, false]);
    assert_eq!(t.resolved_short_radius_scale(3.0), vec![5.0, 3.0]);
}

/// Packmol refuses a short radius that is not shorter than the main one.
#[test]
fn a_short_radius_at_or_above_the_radius_is_rejected() {
    let cell = InsideBoxRestraint::new([0.0; 3], [20.0, 20.0, 20.0], [false; 3]);
    let bad = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], 2)
        .with_name("bad")
        .with_radius(2.0)
        .with_short_radius(2.0)
        .with_restraint(cell);
    let err = GenCanPack::new()
        .with_tolerance(4.0)
        .run(&[bad], 10)
        .expect_err("a short radius >= radius must be rejected");
    let msg = err.to_string();
    assert!(msg.contains("short radius"), "{msg}");
}

/// The global switch, mirroring `use_short_tol` + `short_tol_dist` +
/// `short_tol_scale`.
#[test]
fn global_short_tolerance_applies_to_every_atom() {
    let cell = InsideBoxRestraint::new([0.0; 3], [20.0, 20.0, 20.0], [false; 3]);
    let t = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.5], 8)
        .with_name("a")
        .with_restraint(cell);
    let result = GenCanPack::new()
        .with_tolerance(4.0)
        .with_short_tolerance(2.0, 3.0)
        .with_seed(1)
        .run(&[t], 30)
        .expect("pack");
    assert_eq!(result.natoms(), 8);
}

#[test]
#[should_panic(expected = "smaller than the tolerance")]
fn global_short_tolerance_must_be_below_the_tolerance() {
    let _ = GenCanPack::new()
        .with_tolerance(2.0)
        .with_short_tolerance(4.0, 3.0);
}

// ── the other three keywords, at both script levels ────────────────────────

mod script_extra {
    use molpack::Target;
    use molpack::script::parse;

    fn plan(body: &str) -> molpack::script::ScriptPlan {
        let src = format!(
            "tolerance 4.0\noutput o.xyz\nstructure a.pdb\n  number 1\n\
             inside box 0. 0. 0. 10. 10. 10.\n{body}end structure\n"
        );
        parse(&src)
            .expect("parse")
            .lower(std::path::Path::new("."))
            .expect("lower")
    }

    #[test]
    fn structure_level_keywords_are_parsed() {
        let p = plan("  fscale 0.5\n  short_radius 0.75\n  short_radius_scale 4.0\n");
        let s = &p.structures[0];
        assert_eq!(s.fscale, Some(0.5));
        assert_eq!(s.short_radius, Some(0.75));
        assert_eq!(s.short_radius_scale, Some(4.0));
    }

    #[test]
    fn atom_level_keywords_are_parsed() {
        let p = plan("atoms 2\n  fscale 0.5\n  short_radius 0.75\nend atoms\n");
        let g = &p.structures[0].atom_groups[0];
        assert_eq!(g.fscale, Some(0.5));
        assert_eq!(g.short_radius, Some(0.75));
    }

    #[test]
    fn lowering_applies_every_layer() {
        let p = plan(
            "  fscale 0.5\n  short_radius 0.75\n\
             atoms 3\n  fscale 2.0\n  short_radius_scale 9.0\nend atoms\n",
        );
        let t = p.structures[0].apply(Target::from_coords(
            &[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            &[1.5; 3],
            1,
        ));
        assert_eq!(t.resolved_fscale(), vec![0.5, 0.5, 2.0]);
        assert_eq!(t.resolved_short_radii(1.0), vec![0.75; 3]);
        assert_eq!(t.resolved_short_radius_scale(3.0), vec![3.0, 3.0, 9.0]);
        assert_eq!(t.uses_short_radius(), vec![true; 3]);
    }

    #[test]
    fn a_non_positive_fscale_is_a_parse_error() {
        let src = "tolerance 4.0\noutput o.xyz\nstructure a.pdb\n  number 1\n  fscale 0\n\
                   inside box 0. 0. 0. 10. 10. 10.\nend structure\n";
        assert!(parse(src).is_err());
    }
}
