//! The incremental interaction field: cell bookkeeping, minimum image,
//! exclusions, and the block kinds an insertion has to distinguish.

use super::*;

// ── Section: Task 3 — OverlapField: cell list over the periodic box ────────
//
// Review tests for `src/grow/field.rs` (spec Design §4c): insert/remove
// idempotence (retraction is O(1) removal, not a rebuild), minimum-image
// distances across the periodic boundary, 27-cell stencil de-duplication
// when the grid has fewer than 3 cells per axis, exclusion-list handling,
// and the `hard_scale` softening knob (spec §4f).

/// Insert / remove / re-insert must leave the field exactly consistent:
/// counts, occupancy, and what probe/nearest can see.
#[test]
fn field_insert_remove_idempotent() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [10.0; 3],
        [true; 3],
        vec![1.0; 10],
        (0..10).collect(),
        vec![0; 10],
        2.5,
    );

    let placed: [[F; 3]; 5] = [
        [1.0, 1.0, 1.0],
        [4.0, 4.0, 4.0],
        [7.0, 7.0, 7.0],
        [1.0, 7.0, 1.0],
        [7.0, 1.0, 4.0],
    ];
    for (s, &p) in placed.iter().enumerate() {
        field.insert(s, p);
    }
    assert_eq!(field.n_placed(), 5);
    for (s, &p) in placed.iter().enumerate() {
        assert!(field.is_placed(s), "slot {s} must be placed");
        assert_eq!(field.position(s), p, "slot {s} position");
    }
    for s in 5..10 {
        assert!(!field.is_placed(s), "slot {s} was never inserted");
    }

    field.remove(1);
    field.remove(3);
    assert_eq!(field.n_placed(), 3);
    assert!(!field.is_placed(1) && !field.is_placed(3));
    field.remove(1); // documented no-op on an empty slot
    assert_eq!(field.n_placed(), 3, "double remove must be a no-op");

    // Re-insert the two slots elsewhere: probe/nearest must see the new
    // positions and have forgotten the old ones.
    field.insert(1, [2.2, 1.0, 1.0]);
    field.insert(3, [9.4, 1.0, 1.0]);
    assert_eq!(field.n_placed(), 5);
    assert!(field.is_placed(1) && field.is_placed(3));
    let d = field.nearest(9, [2.2, 1.0, 1.0], &[]);
    assert!(
        d < 1e-12,
        "nearest at slot 1's new position = {d}, must be 0"
    );
    assert_eq!(
        field.probe(9, [4.0, 4.0, 4.0], &[], 1.0, 0.5, 1e9),
        Probe::Room(0.0),
        "slot 1's old position must be empty after re-insertion"
    );
    assert_eq!(
        field.probe(9, [2.2, 1.0, 1.0], &[], 1.0, 0.5, 1e9),
        Probe::Blocked,
        "slot 1's new position must be hard-blocked"
    );

    for s in 0..5 {
        field.remove(s);
    }
    assert_eq!(field.n_placed(), 0);
    for p in [[0.0; 3], [5.0; 3], [9.9, 0.1, 5.0]] {
        assert_eq!(
            field.probe(9, p, &[], 1.0, 0.5, 1e9),
            Probe::Room(0.0),
            "an empty field must answer Room(0.0) everywhere"
        );
    }
}

/// Distances must use the minimum image: 0.5 and 9.9 in a 10 Å periodic box
/// are 0.6 Å apart, well inside the 2.0 Å contact.
#[test]
fn field_minimum_image_across_boundary() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [10.0; 3],
        [true; 3],
        vec![1.0, 1.0],
        vec![0, 1],
        vec![0, 0],
        2.5,
    );
    field.insert(0, [0.5, 5.0, 5.0]);
    assert_eq!(
        field.probe(1, [9.9, 5.0, 5.0], &[], 1.0, 0.4, 1e9),
        Probe::Blocked,
        "0.6 Å across the periodic boundary is inside the 2.0 Å contact"
    );
    let d = field.nearest(1, [9.9, 5.0, 5.0], &[]);
    assert!(
        (d - 0.6).abs() < 1e-12,
        "nearest across the boundary = {d}, expected 0.6"
    );
}

/// Box 4³ with cutoff 2.5 → a single cell per axis: all 27 stencil offsets
/// alias the same cell, and a naive stencil would count the one neighbour up
/// to 27 times. The penalty is hand-derived and binary-exact, so any double
/// counting shows up as an integer multiple.
#[test]
fn field_small_box_stencil_dedup() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [4.0; 3],
        [true; 3],
        vec![0.5, 0.5],
        vec![0, 1],
        vec![0, 0],
        2.5,
    );
    field.insert(0, [1.0, 2.0, 2.0]);
    // d = 1.5, contact = 1.0, shell = contact + soft_shell = 2.0. Exactly one
    // pair inside the soft shell: gap = (shell² − d²)/shell² = (4 − 2.25)/4
    // = 0.4375, penalty = gap² = 0.19140625.
    let probe = field.probe(1, [2.5, 2.0, 2.0], &[], 1.0, 1.0, 1e9);
    let Probe::Room(penalty) = probe else {
        panic!("d = 1.5 > contact = 1.0 must not be Blocked, got {probe:?}");
    };
    assert!(
        (penalty - 0.19140625).abs() < 1e-12,
        "penalty = {penalty}, expected exactly 0.19140625 — 2× means the \
         aliased cell was visited twice"
    );
    let d = field.nearest(1, [2.5, 2.0, 2.0], &[]);
    assert!((d - 1.5).abs() < 1e-12, "nearest = {d}, expected 1.5");
}

/// Same-molecule template atoms on the exclusion list are neither blocking
/// nor charged; a same-molecule atom OFF the list still blocks.
#[test]
fn field_exclusions_respected() {
    // Three slots of the same molecule (mol 7): template atoms 0, 1 and 4.
    let mut field = OverlapField::new(
        [0.0; 3],
        [20.0; 3],
        [true; 3],
        vec![1.0; 3],
        vec![7; 3],
        vec![0, 1, 4],
        2.5,
    );
    let excluded: &[u32] = &[0, 1, 2, 3]; // template atom 4 is NOT excluded
    field.insert(0, [5.0, 5.0, 5.0]);
    assert_eq!(
        field.probe(1, [5.5, 5.0, 5.0], excluded, 1.0, 0.5, 1e9),
        Probe::Room(0.0),
        "deep inside slot 0's hard core, but its template atom is excluded → Room"
    );
    field.insert(2, [6.5, 5.0, 5.0]);
    assert_eq!(
        field.probe(1, [6.0, 5.0, 5.0], excluded, 1.0, 0.5, 1e9),
        Probe::Blocked,
        "template atom 4 is not on the exclusion list and violates → Blocked"
    );
}

/// `hard_scale` shrinks the refusal core (the §4f soften fallback): a pair at
/// d = 1.8 with contact 2.0 is Blocked at scale 1.0 and merely charged at 0.8.
#[test]
fn field_hard_scale_softens() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [10.0; 3],
        [true; 3],
        vec![1.0, 1.0],
        vec![0, 1],
        vec![0, 0],
        3.0,
    );
    field.insert(0, [5.0, 5.0, 5.0]);
    let p = [6.8, 5.0, 5.0];
    assert_eq!(
        field.probe(1, p, &[], 1.0, 0.4, 1e9),
        Probe::Blocked,
        "hard_scale 1.0: d = 1.8 < contact = 2.0"
    );
    match field.probe(1, p, &[], 0.8, 0.4, 1e9) {
        Probe::Room(penalty) => assert!(
            penalty > 0.0,
            "d = 1.8 sits inside the degraded soft shell (1.6..2.0): allowed \
             but charged, got penalty {penalty}"
        ),
        Probe::Blocked => panic!("hard_scale 0.8 must admit d = 1.8 > 1.6"),
    }
}

/// `block_kind` classifies a hard-core hit as same-molecule self-block vs
/// inter-chain. Every radius is 1.0 Å so contact at `hard_scale = 1.0` is
/// 2.0 Å. Probe slot 1 (mol 7, template atom 1) sits 0.5 Å from each core
/// neighbour — well inside that contact. Template atom 4 is the 1-5 partner
/// (not excluded); template atoms 0–3 are on the exclusion list.
#[test]
fn field_block_kind_self_vs_inter() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [20.0; 3],
        [true; 3],
        vec![1.0; 4],
        vec![7, 7, 7, 8],
        vec![0, 1, 4, 0],
        2.5,
    );
    let excluded: &[u32] = &[0, 1, 2, 3];
    let p = [5.5, 5.0, 5.0];
    let hard_scale = 1.0;

    // Empty field / no hard-core hit → None.
    assert_eq!(
        field.block_kind(1, p, excluded, hard_scale),
        None,
        "empty field: no neighbour, so no hard-core hit"
    );

    // Same-molecule neighbour ON the exclusion list (template atom 0).
    field.insert(0, [5.0, 5.0, 5.0]);
    assert_eq!(
        field.block_kind(1, p, excluded, hard_scale),
        None,
        "d = 0.5 inside contact 2.0, but template atom 0 is excluded → None"
    );
    field.remove(0);

    // Far same-molecule non-excluded neighbour: d = 3.0 > contact 2.0.
    field.insert(2, [8.5, 5.0, 5.0]);
    assert_eq!(
        field.block_kind(1, p, excluded, hard_scale),
        None,
        "d = 3.0 > contact 2.0: no hard-core hit → None"
    );
    field.remove(2);

    // SelfBlocked: same mol, template atom 4 not excluded, inside the core.
    field.insert(2, [5.0, 5.0, 5.0]);
    assert_eq!(
        field.block_kind(1, p, excluded, hard_scale),
        Some(BlockKind::SelfBlocked),
        "d = 0.5, mol 7 == mol 7, template atom 4 not excluded → SelfBlocked"
    );
    field.remove(2);

    // InterChain: different molecule, inside the hard core.
    field.insert(3, [5.0, 5.0, 5.0]);
    assert_eq!(
        field.block_kind(1, p, excluded, hard_scale),
        Some(BlockKind::InterChain),
        "d = 0.5, mol 8 ≠ mol 7 → InterChain"
    );
    field.remove(3);

    // SelfBlocked wins: insert the same-mol neighbour first, then the
    // inter-chain neighbour so the cell-list head is InterChain (a first-hit
    // scan would otherwise return InterChain).
    field.insert(2, [5.0, 5.0, 5.0]);
    field.insert(3, [5.5, 5.5, 5.0]);
    assert_eq!(
        field.block_kind(1, p, excluded, hard_scale),
        Some(BlockKind::SelfBlocked),
        "same-mol non-excluded core hit must beat an inter-chain core hit \
         (d_self = 0.5, d_inter = 0.5, both < contact 2.0)"
    );
}

/// Cell-occupancy bookkeeping behind void-biased seeding: the empty-cell
/// count tracks insert/remove exactly, and `empty_cell_point` never offers
/// a point inside an occupied cell.
#[test]
fn field_empty_cell_bookkeeping() {
    // 4³ box, cutoff 2 → 2×2×2 grid, cell edge 2.
    let mut field = OverlapField::new(
        [0.0; 3],
        [4.0; 3],
        [true; 3],
        vec![0.5; 4],
        vec![0, 1, 2, 3],
        vec![0; 4],
        2.0,
    );
    assert_eq!(field.n_empty_cells(), 8);

    field.insert(0, [1.0, 1.0, 1.0]); // cell (0,0,0)
    assert_eq!(field.n_empty_cells(), 7);
    field.insert(1, [1.5, 1.5, 1.5]); // same cell
    assert_eq!(field.n_empty_cells(), 7);

    // Every pick must avoid the occupied cell [0,2)³.
    for i in 0..64 {
        let u0 = (i as F + 0.5) / 64.0;
        let p = field
            .empty_cell_point([u0, 0.5, 0.5, 0.5])
            .expect("empty cells exist");
        assert!(
            p.iter().any(|&x| x >= 2.0),
            "sampled point {p:?} landed in the occupied cell"
        );
    }

    field.remove(1);
    assert_eq!(field.n_empty_cells(), 7, "cell still holds slot 0");
    field.remove(0);
    assert_eq!(field.n_empty_cells(), 8);

    // A fully occupied grid offers no cavity.
    let mut full = OverlapField::new(
        [0.0; 3],
        [2.0; 3],
        [true; 3],
        vec![0.5],
        vec![0],
        vec![0],
        2.0,
    );
    full.insert(0, [1.0, 1.0, 1.0]);
    assert_eq!(full.n_empty_cells(), 0);
    assert!(full.empty_cell_point([0.3, 0.5, 0.5, 0.5]).is_none());
}
