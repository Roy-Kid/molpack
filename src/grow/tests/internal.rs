//! `InternalTree`: the bond-graph → internal-coordinate decomposition, and
//! the rebuild that has to invert it exactly.

use super::*;

/// Rebuilding a linear chain with the template's own torsion values must
/// reproduce the template coordinates exactly (spec Task 2 (i)).
#[test]
fn internal_roundtrip_linear() {
    let coords = zigzag_coords(12, 1.53);
    let tree = assert_roundtrip(&coords, &chain_bonds(12), "linear 12-chain");
    // Bonds 1-2 … 9-10 rotate; the two terminal bonds never do.
    assert_eq!(tree.n_vars(), 9, "linear 12-chain: 9 free torsions");
    assert_eq!(
        tree.n_steps(),
        9,
        "one step per free torsion, no rigid preamble"
    );
}

/// Same round-trip for a branched molecule: branch atoms are placed off the
/// backbone and must come back exactly, and the terminal branch bonds must
/// not add free variables.
#[test]
fn internal_roundtrip_branched() {
    let (coords, bonds) = branched_parts();
    let tree = assert_roundtrip(&coords, &bonds, "branched");
    // Backbone bonds 1-2 … 7-8 rotate; both chain-end bonds and both
    // one-atom branch bonds are terminal.
    assert_eq!(tree.n_vars(), 7, "branched: 7 free torsions");
}

/// A ring template cannot build an internal-coordinate tree: it is
/// `GrowError::RingTemplate`, not a silent open-ring decomposition.
#[test]
fn internal_roundtrip_ring() {
    let (coords, bonds) = ring_tail_parts();
    let err = InternalTree::from_frame(
        &frame_from_parts(&coords, &bonds),
        &BondDistanceWeights::from_exclusion_depth(3),
    )
    .expect_err("a ring template must not construct InternalTree");
    assert!(
        matches!(err, GrowError::RingTemplate),
        "expected GrowError::RingTemplate, got {err:?}"
    );
}

/// A molecule with no rotatable bonds (tetrahedral star: every bond is
/// terminal) still decomposes — zero variables, at least one step, exact
/// round-trip with an empty variable vector.
#[test]
fn internal_rigid_molecule_no_vars() {
    let s = 1.53 as F / (3.0 as F).sqrt();
    let coords = [[0.0, 0.0, 0.0], [s, s, s], [s, -s, -s], [-s, s, -s]];
    let bonds = [(0, 1), (0, 2), (0, 3)];
    let tree = assert_roundtrip(&coords, &bonds, "rigid star");
    assert_eq!(tree.n_vars(), 0, "terminal bonds are never rotatable");
    assert!(tree.n_steps() >= 1, "the rigid remainder is still a step");
}

/// Bonded-geometry detector on an acyclic template: under ARBITRARY
/// free-variable values, every template bond length and every bonded angle
/// must still match the template, and all coordinates must be finite. Free
/// torsions legitimately change and are not checked. Ring-closure detection
/// is the named `GrowError::RingTemplate` refusal (`internal_roundtrip_ring`,
/// `ring_template_is_refused`); this test only randomizes acyclic templates.
#[test]
fn internal_random_vars_preserve_bonded_geometry() {
    let pi = std::f64::consts::PI as F;
    let mut rng = SmallRng::seed_from_u64(20260828);
    for (what, (coords, bonds)) in [("branched", branched_parts())] {
        let tree = InternalTree::from_frame(
            &frame_from_parts(&coords, &bonds),
            &BondDistanceWeights::from_exclusion_depth(3),
        )
        .unwrap_or_else(|e| panic!("{what}: template must decompose, got {e}"));
        assert!(
            tree.n_vars() > 0,
            "{what}: needs free torsions to randomize"
        );
        let angles = bonded_angle_triples(coords.len(), &bonds);
        for trial in 0..5 {
            let vars: Vec<F> = (0..tree.n_vars())
                .map(|_| (rng.random::<f64>() as F * 2.0 - 1.0) * pi)
                .collect();
            let rebuilt = rebuild_coords(&tree, &coords, &vars);
            for (i, p) in rebuilt.iter().enumerate() {
                assert!(
                    p.iter().all(|v| v.is_finite()),
                    "{what} trial {trial}: atom {i} is not finite: {p:?}"
                );
            }
            for &(i, j) in &bonds {
                let (i, j) = (i as usize, j as usize);
                let want = vdist(coords[i], coords[j]);
                let got = vdist(rebuilt[i], rebuilt[j]);
                assert!(
                    (got - want).abs() < 1e-9,
                    "{what} trial {trial}: bond ({i},{j}) = {got} vs template {want} — \
                     a ring-closure or otherwise non-free bond classified as a free \
                     variable breaks exactly this"
                );
            }
            for &(i, j, k) in &angles {
                let want = vangle(coords[i], coords[j], coords[k]);
                let got = vangle(rebuilt[i], rebuilt[j], rebuilt[k]);
                assert!(
                    (got - want).abs() < 1e-9,
                    "{what} trial {trial}: angle ({i},{j},{k}) = {got} rad vs template {want} rad"
                );
            }
        }
    }
}

/// The 1-`EXCLUDE_BONDS` (depth 3, AA convention) exclusion table: sorted,
/// self included, 1-2/1-3/1-4 partners in, 1-5 partners out.
#[test]
fn internal_exclusions_depth() {
    let tree = InternalTree::from_frame(
        &chain_frame(12, 1.53, true),
        &BondDistanceWeights::from_exclusion_depth(3),
    )
    .expect("linear chain decomposes");
    assert_eq!(
        tree.exclusions(0),
        &[0u32, 1, 2, 3][..],
        "atom 0 of a linear chain excludes exactly self + its 1-2/1-3/1-4 partners"
    );
    assert!(
        !tree.exclusions(0).contains(&4),
        "1-5 partners must be scored, not excluded — the chain must not thread itself"
    );
}

/// Growth compiles a binary skip table. A fractional 1-4 weight is legal on
/// `Target.special_bonds` but `InternalTree::from_frame` refuses it by name
/// before building exclusions. Slot 2 is 1-4.
#[test]
fn internal_tree_rejects_non_binary_special_bond() {
    let weights = BondDistanceWeights::new(vec![0.0, 0.0, 0.5, 1.0])
        .expect("fractional 1-4 is a legal BondDistanceWeights table");
    let err = InternalTree::from_frame(&chain_frame(5, 1.53, true), &weights)
        .expect_err("a fractional special-bonds weight must not compile a growth tree");
    match &err {
        GrowError::NonBinarySpecialBond { index, weight } => {
            assert_eq!(*index, 2, "slot 2 is the 1-4 weight");
            assert_eq!(
                weight.to_bits(),
                (0.5 as F).to_bits(),
                "reported weight must be the stored 0.5"
            );
        }
        other => panic!(
            "expected GrowError::NonBinarySpecialBond {{ index: 2, weight: 0.5 }}, got {other:?}"
        ),
    }
    let msg = err.to_string();
    assert!(
        msg.contains("1-4"),
        "Display must name the 1-4 pair, got: {msg}"
    );
    assert!(
        msg.contains("0.5"),
        "Display must include the offending weight, got: {msg}"
    );
    assert!(
        msg.contains("with_atom_radius"),
        "Display must point at with_atom_radius, got: {msg}"
    );
}

/// Same 5-bead chain as `(0,1)…(3,4)` with bonds written out of file order
/// still decomposes; skip sets include the root.
#[test]
fn internal_tree_shuffled_bonds_still_decomposes() {
    let coords = zigzag_coords(5, 1.53);
    let bonds = [(3, 4), (1, 2), (0, 1), (2, 3)];
    let frame = frame_from_parts(&coords, &bonds);
    let tree = InternalTree::from_frame(&frame, &BondDistanceWeights::from_exclusion_depth(3))
        .expect("shuffled linear-chain bonds must still decompose");
    assert_eq!(tree.n_atoms(), 5, "5-bead chain");
    let root = tree.seed_atoms()[0];
    assert!(
        tree.exclusions(root).contains(&(root as u32)),
        "exclusions({root}) must include the root"
    );
}
