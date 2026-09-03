//! Integration tests for the shared bond-graph leaf `src/topology.rs`
//! (`.claude/specs/stage-pipeline-01-topology.md`, Tasks 1 and 6).
//!
//! One file, one module: everything here exercises the public surface of
//! `molpack::topology` — `Topology::{from_frame, natoms, bonds, neighbors,
//! require_connected, exclusions}`, the free function `frame_positions`, and
//! `TopologyError`. Geometry is synthesized in-file, following the
//! `frame_from_parts` / `chain_frame` pattern of `tests/grow.rs`; no solver,
//! no I/O, no RNG.
//!
//! Two facts here are load-bearing for the extraction and are asserted
//! explicitly, because both are places where a "tidier" implementation would
//! silently change growth:
//!
//! 1. **Adjacency keeps bond-file insertion order** — never sorted. `pick_ref`
//!    and `bfs_order` in `src/grow/internal.rs` take the *first* neighbor, so
//!    sorting would reshape the growth tree.
//! 2. **`exclusions(depth)` includes the root atom** — the pre-refactor doc
//!    comment said "excludes self", the pre-refactor *code* seeds `touched`
//!    with the root. The spec keeps the code semantics (spec Design, ac-005).

use molpack::F;
use molpack::topology::{Topology, TopologyError, frame_positions};
use molrs::store::block::Block;
use molrs::store::frame::Frame;
use ndarray::Array1;

// ── synthesized geometry (mirrors the helpers in `tests/grow.rs`) ──────────

/// Planar zigzag bead-chain coordinates: tetrahedral (109.5°) bond angles in
/// the x–z plane, all torsions trans.
fn zigzag_coords(n: usize, bond_len: F) -> Vec<[F; 3]> {
    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());
    (0..n)
        .map(|i| [i as F * dx, 0.0, if i % 2 == 0 { 0.0 } else { dz }])
        .collect()
}

/// The `(i, i+1)` bond list of a linear chain, in file order.
fn chain_bonds(n: usize) -> Vec<(u32, u32)> {
    (0..n as u32 - 1).map(|i| (i, i + 1)).collect()
}

/// Coordinates + explicit bond list as a `molrs::Frame`: atoms block with
/// x/y/z columns and (unless `bonds` is empty) a bonds block with
/// atomi/atomj. Deliberately NO `bond_type` column — this is the shape a PDB
/// CONECT list or a hand-built coarse-grain frame has.
fn frame_from_parts(coords: &[[F; 3]], bonds: &[(u32, u32)]) -> Frame {
    let mut atoms = Block::new();
    for (name, k) in [("x", 0), ("y", 1), ("z", 2)] {
        let col: Vec<F> = coords.iter().map(|p| p[k]).collect();
        atoms
            .insert(name, Array1::from_vec(col).into_dyn())
            .expect("coordinate column");
    }

    let mut frame = Frame::new();
    frame.insert("atoms", atoms);

    if !bonds.is_empty() {
        frame.insert("bonds", bond_block(bonds));
    }

    frame
}

/// A bonds block with `atomi` / `atomj` columns in the given order.
fn bond_block(bonds: &[(u32, u32)]) -> Block {
    let mut block = Block::new();
    let ai: Vec<u32> = bonds.iter().map(|&(i, _)| i).collect();
    let aj: Vec<u32> = bonds.iter().map(|&(_, j)| j).collect();
    block
        .insert("atomi", Array1::from_vec(ai).into_dyn())
        .expect("atomi column");
    block
        .insert("atomj", Array1::from_vec(aj).into_dyn())
        .expect("atomj column");
    block
}

/// Zigzag bead chain as a `molrs::Frame` (see [`zigzag_coords`] /
/// [`frame_from_parts`]).
fn chain_frame(n: usize, bond_len: F, with_bonds: bool) -> Frame {
    let bonds = if with_bonds {
        chain_bonds(n)
    } else {
        Vec::new()
    };
    frame_from_parts(&zigzag_coords(n, bond_len), &bonds)
}

/// 10-atom zigzag backbone with two one-atom side branches off backbone atoms
/// 3 and 6 (atoms 10 and 11, bonded last — insertion order matters).
fn branched_parts() -> (Vec<[F; 3]>, Vec<(u32, u32)>) {
    let mut coords = zigzag_coords(10, 1.53);
    let mut bonds = chain_bonds(10);
    let c3 = coords[3];
    coords.push([c3[0], c3[1] + 1.44, c3[2] + 0.51]); // atom 10
    bonds.push((3, 10));
    let c6 = coords[6];
    coords.push([c6[0], c6[1] - 1.44, c6[2] - 0.51]); // atom 11
    bonds.push((6, 11));
    (coords, bonds)
}

/// Chair-ish (non-planar) 6-ring with ~1.46 Å bonds plus a 4-atom zigzag tail
/// off ring atom 0. The ring closure bond `(5, 0)` is the tenth entry and is
/// what makes the cyclic exclusion assertions below non-trivial.
fn ring_tail_parts() -> (Vec<[F; 3]>, Vec<(u32, u32)>) {
    let pi = std::f64::consts::PI as F;
    let (r, h) = (1.46 as F, 0.25 as F);
    let mut coords: Vec<[F; 3]> = (0..6)
        .map(|i| {
            let phi = i as F * pi / 3.0;
            [
                r * phi.cos(),
                r * phi.sin(),
                if i % 2 == 0 { h } else { -h },
            ]
        })
        .collect();
    let mut bonds: Vec<(u32, u32)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
    let mut tip = coords[0];
    for t in 0..4u32 {
        tip = [
            tip[0] + 1.25,
            tip[1],
            tip[2] + if t % 2 == 0 { 0.9 } else { -0.9 },
        ];
        coords.push(tip); // atoms 6, 7, 8, 9
        bonds.push(if t == 0 { (0, 6) } else { (5 + t, 6 + t) });
    }
    (coords, bonds)
}

/// A frame whose atoms block carries only x and y — the "unreadable atoms
/// block" input, distinct from a frame with no atoms block at all.
fn frame_missing_z_column(coords: &[[F; 3]]) -> Frame {
    let mut atoms = Block::new();
    for (name, k) in [("x", 0), ("y", 1)] {
        let col: Vec<F> = coords.iter().map(|p| p[k]).collect();
        atoms
            .insert(name, Array1::from_vec(col).into_dyn())
            .expect("coordinate column");
    }
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame
}

/// Every exclusion list must be strictly ascending and contain its own root.
fn assert_root_inclusive_sorted(lists: &[Vec<u32>], natoms: usize, what: &str) {
    assert_eq!(lists.len(), natoms, "{what}: one list per atom");
    for (root, list) in lists.iter().enumerate() {
        assert!(
            list.contains(&(root as u32)),
            "{what}: exclusions[{root}] = {list:?} must include the root atom itself"
        );
        assert!(
            list.windows(2).all(|w| w[0] < w[1]),
            "{what}: exclusions[{root}] = {list:?} must be strictly ascending"
        );
    }
}

// ── Basics — from_frame on linear / branched / ring templates ─────────────

/// A 12-bead linear chain reads back as 12 atoms, its 11 bonds in file order,
/// and CSR adjacency of degree 1 / 2 / 1.
#[test]
fn topology_from_frame_linear_chain() {
    let topo = Topology::from_frame(&chain_frame(12, 1.53, true)).expect("linear 12-chain");

    // The leaf is re-exported from the crate root (spec Task 2).
    let _: molpack::Topology = topo.clone();

    assert_eq!(topo.natoms(), 12, "12-bead chain has 12 atoms");
    let want: Vec<(u32, u32)> = chain_bonds(12);
    assert_eq!(topo.bonds(), want.as_slice(), "bonds in bond-file order");
    assert_eq!(topo.neighbors(0), [1], "chain head has one neighbor");
    assert_eq!(topo.neighbors(5), [4, 6], "interior bead, insertion order");
    assert_eq!(topo.neighbors(11), [10], "chain tail has one neighbor");
}

/// The branched template: two side atoms bonded *after* the backbone, so the
/// branch neighbor is last in the CSR slice of its backbone atom.
#[test]
fn topology_from_frame_branched_template() {
    let (coords, bonds) = branched_parts();
    let topo = Topology::from_frame(&frame_from_parts(&coords, &bonds)).expect("branched template");

    assert_eq!(topo.natoms(), 12, "10 backbone + 2 branch atoms");
    assert_eq!(topo.bonds(), bonds.as_slice(), "bonds in bond-file order");
    assert_eq!(
        topo.neighbors(3),
        [2, 4, 10],
        "backbone bonds first, branch bond last — file order, not sorted-by-luck"
    );
    assert_eq!(topo.neighbors(6), [5, 7, 11], "second branch point");
    assert_eq!(topo.neighbors(10), [3], "branch atom is terminal");
}

/// The ring + tail template: the closure bond `(5, 0)` is a real edge and the
/// ring atoms keep degree 2 (3 where the tail attaches).
#[test]
fn topology_from_frame_ring_template() {
    let (coords, bonds) = ring_tail_parts();
    let topo = Topology::from_frame(&frame_from_parts(&coords, &bonds)).expect("ring template");

    assert_eq!(topo.natoms(), 10, "6-ring + 4-atom tail");
    assert_eq!(topo.bonds(), bonds.as_slice(), "bonds in bond-file order");
    assert_eq!(
        topo.bonds()[5],
        (5, 0),
        "the ring closure bond is stored as written, not normalized"
    );
    assert_eq!(
        topo.neighbors(0),
        [1, 5, 6],
        "ring neighbor, closure neighbor, then the tail"
    );
    assert_eq!(topo.neighbors(3), [2, 4], "plain ring atom");
    assert_eq!(topo.neighbors(9), [8], "tail tip");
}

/// **Insertion order, never sorted.** The bonds block is deliberately
/// reversed / shuffled — `(3,2), (0,1), (1,2)` — so a sorting implementation
/// and an insertion-order implementation disagree on the *first* neighbor of
/// atom 2, which is exactly what `pick_ref` / `bfs_order` read.
#[test]
fn topology_neighbors_preserve_bond_file_order() {
    let shuffled = [(3, 2), (0, 1), (1, 2)];
    let topo = Topology::from_frame(&frame_from_parts(&zigzag_coords(4, 1.53), &shuffled))
        .expect("shuffled 4-chain");

    assert_eq!(
        topo.bonds(),
        shuffled.as_slice(),
        "bonds() replays the bonds block verbatim"
    );
    assert_eq!(
        topo.neighbors(2),
        [3, 1],
        "atom 2 met atom 3 first (bond 0) and atom 1 second (bond 2) — NOT sorted"
    );
    assert_eq!(
        topo.neighbors(2).first().copied(),
        Some(3),
        "the first neighbor is the growth tree's reference pick"
    );
    let mut sorted = topo.neighbors(2).to_vec();
    sorted.sort_unstable();
    assert_ne!(
        topo.neighbors(2),
        sorted.as_slice(),
        "this fixture only bites if the adjacency is left unsorted"
    );
    assert_eq!(topo.neighbors(1), [0, 2], "atom 1: bond 1 then bond 2");
    assert_eq!(topo.neighbors(3), [2], "atom 3 saw only the first bond");
}

// ── Domain — exclusions are root-inclusive and ascending ──────────────────

/// The load-bearing semantics (spec Design / ac-005): every list contains its
/// own root and is strictly ascending, at every depth. Growth uses these
/// lists directly as the `OverlapField` skip set, so dropping the root would
/// silently change which atom pairs get scored.
#[test]
fn topology_exclusions_are_root_inclusive_and_sorted() {
    let topo = Topology::from_frame(&chain_frame(12, 1.53, true)).expect("linear 12-chain");
    for depth in [1usize, 2, 3] {
        let lists = topo.exclusions(depth);
        assert_root_inclusive_sorted(&lists, 12, &format!("C12 depth {depth}"));
    }
}

/// Literal spot checks of the BFS shells on the 12-bead chain at depth 1/2/3.
#[test]
fn topology_exclusions_depth_literals_on_c12_chain() {
    let topo = Topology::from_frame(&chain_frame(12, 1.53, true)).expect("linear 12-chain");

    let d1 = topo.exclusions(1);
    assert_eq!(d1[5], vec![4, 5, 6], "depth 1: bonded partners + root");
    assert_eq!(d1[0], vec![0, 1], "depth 1 at the chain head");
    assert_eq!(d1[11], vec![10, 11], "depth 1 at the chain tail");

    let d2 = topo.exclusions(2);
    assert_eq!(d2[5], vec![3, 4, 5, 6, 7], "depth 2: 1-2 and 1-3 partners");
    assert_eq!(d2[0], vec![0, 1, 2], "depth 2 at the chain head");

    let d3 = topo.exclusions(3);
    assert_eq!(d3[0], vec![0, 1, 2, 3], "depth 3 at the chain head");
    assert_eq!(
        d3[5],
        vec![2, 3, 4, 5, 6, 7, 8],
        "depth 3: three bonds each way, root included"
    );
    assert_eq!(d3[11], vec![8, 9, 10, 11], "depth 3 at the chain tail");
}

/// A cycle reaches back on itself: at depth 3 ring atom 0 excludes the whole
/// ring plus three tail atoms, because atom 3 is three bonds away *both* ways.
#[test]
fn topology_exclusions_close_around_a_ring() {
    let (coords, bonds) = ring_tail_parts();
    let topo = Topology::from_frame(&frame_from_parts(&coords, &bonds)).expect("ring template");

    assert_root_inclusive_sorted(&topo.exclusions(2), 10, "ring depth 2");
    assert_eq!(
        topo.exclusions(1)[0],
        vec![0, 1, 5, 6],
        "depth 1 at the substituted ring atom: both ring bonds + the tail bond"
    );
    assert_eq!(
        topo.exclusions(2)[0],
        vec![0, 1, 2, 4, 5, 6, 7],
        "depth 2 skips ring atom 3, which is three bonds away either way"
    );
    assert_eq!(
        topo.exclusions(3)[0],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
        "depth 3 closes the ring: everything but the tail tip"
    );
}

/// A branch atom's shells run through its backbone attachment point only.
#[test]
fn topology_exclusions_on_branched_template() {
    let (coords, bonds) = branched_parts();
    let topo = Topology::from_frame(&frame_from_parts(&coords, &bonds)).expect("branched template");

    assert_root_inclusive_sorted(&topo.exclusions(3), 12, "branched depth 3");
    assert_eq!(
        topo.exclusions(1)[3],
        vec![2, 3, 4, 10],
        "depth 1 at the branch point: two backbone partners + the branch"
    );
    assert_eq!(
        topo.exclusions(3)[10],
        vec![1, 2, 3, 4, 5, 10],
        "depth 3 from the branch atom, funnelled through backbone atom 3"
    );
    assert_eq!(
        topo.exclusions(3)[3],
        vec![0, 1, 2, 3, 4, 5, 6, 10],
        "depth 3 at the branch point reaches both backbone directions"
    );
}

// ── Edge cases — the named rejections ─────────────────────────────────────

/// A frame with a bonds block but no atoms block: `NoAtomsBlock`.
#[test]
fn topology_from_frame_without_atoms_block_is_no_atoms_block() {
    let mut frame = Frame::new();
    frame.insert("bonds", bond_block(&chain_bonds(6)));

    let err = Topology::from_frame(&frame).expect_err("no atoms block must be refused");
    assert!(
        matches!(err, TopologyError::NoAtomsBlock),
        "expected NoAtomsBlock, got {err:?}"
    );
}

/// An atoms block missing the `z` column is just as unreadable as a missing
/// block — the atom count comes from the same x/y/z check.
#[test]
fn topology_from_frame_with_unreadable_atoms_block_is_no_atoms_block() {
    let frame = frame_missing_z_column(&zigzag_coords(6, 1.53));

    let err = Topology::from_frame(&frame).expect_err("missing z column must be refused");
    assert!(
        matches!(err, TopologyError::NoAtomsBlock),
        "expected NoAtomsBlock, got {err:?}"
    );
}

/// Atoms but no bonds block at all: `NoBonds`.
#[test]
fn topology_from_frame_without_bonds_block_is_no_bonds() {
    let err = Topology::from_frame(&chain_frame(6, 1.53, false))
        .expect_err("a template without connectivity must be refused");
    assert!(
        matches!(err, TopologyError::NoBonds),
        "expected NoBonds, got {err:?}"
    );
}

/// A present-but-empty bonds block is the same rejection as no block.
#[test]
fn topology_from_frame_with_empty_bonds_block_is_no_bonds() {
    let mut frame = chain_frame(6, 1.53, false);
    frame.insert("bonds", bond_block(&[]));

    let err = Topology::from_frame(&frame).expect_err("an empty bonds block must be refused");
    assert!(
        matches!(err, TopologyError::NoBonds),
        "expected NoBonds, got {err:?}"
    );
}

/// A bond naming an atom index `>= natoms` reports both endpoints and the
/// atom count, so the user can find the offending line.
#[test]
fn topology_from_frame_with_out_of_range_bond_reports_endpoints() {
    let mut bonds = chain_bonds(6);
    bonds.push((4, 6)); // atom 6 does not exist in a 6-atom template
    let frame = frame_from_parts(&zigzag_coords(6, 1.53), &bonds);

    let err = Topology::from_frame(&frame).expect_err("out-of-range bond must be refused");
    match err {
        TopologyError::BondOutOfRange { a, b, n } => {
            assert_eq!(a, 4, "first endpoint");
            assert_eq!(b, 6, "offending endpoint");
            assert_eq!(n, 6, "template atom count");
        }
        other => panic!("expected BondOutOfRange, got {other:?}"),
    }
}

/// A self-loop `(a, a)` is dropped silently — not an error, not a bond, and
/// not a neighbor of itself.
#[test]
fn topology_from_frame_drops_self_loop_bonds() {
    let mut bonds = chain_bonds(6);
    bonds.insert(2, (2, 2)); // between (1,2) and (2,3) in file order
    let topo = Topology::from_frame(&frame_from_parts(&zigzag_coords(6, 1.53), &bonds))
        .expect("self loop is not an error");

    assert_eq!(topo.natoms(), 6, "the self loop adds no atom");
    assert_eq!(
        topo.bonds(),
        chain_bonds(6).as_slice(),
        "the self loop is dropped, the remaining order is untouched"
    );
    assert_eq!(
        topo.neighbors(2),
        [1, 3],
        "an atom is never its own neighbor"
    );
}

/// `require_connected` accepts a chain in which every atom carries a bond.
#[test]
fn topology_require_connected_accepts_a_chain() {
    let topo = Topology::from_frame(&chain_frame(6, 1.53, true)).expect("linear 6-chain");

    assert!(
        topo.require_connected().is_ok(),
        "every atom of a linear chain has a bond"
    );
}

/// An atom with no bond at all is a `Disconnected` template — the zero-degree
/// screen that `InternalTree::from_frame_with_depth` runs before its BFS.
#[test]
fn topology_require_connected_rejects_an_isolated_atom() {
    // 6 atoms, but only 0-1-2-3-4 are bonded: atom 5 is isolated.
    let frame = frame_from_parts(&zigzag_coords(6, 1.53), &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    let topo = Topology::from_frame(&frame).expect("the frame itself is readable");

    assert_eq!(topo.natoms(), 6, "the isolated atom still counts");
    assert!(topo.neighbors(5).is_empty(), "atom 5 has no neighbor");
    let err = topo
        .require_connected()
        .expect_err("an isolated atom must be refused");
    assert!(
        matches!(err, TopologyError::Disconnected),
        "expected Disconnected, got {err:?}"
    );
}

// ── Error surface — trait bounds and verbatim message text ────────────────

/// `TopologyError` is a real error type: `std::error::Error + Display`, so
/// `GrowError::Topology(e)` can pass it through with `write!(f, "{e}")`.
#[test]
fn topology_error_is_a_std_error() {
    fn assert_std_error<E: std::error::Error + std::fmt::Display + std::fmt::Debug>(
        e: &E,
    ) -> String {
        e.to_string()
    }

    let msg = assert_std_error(&TopologyError::Disconnected);
    assert!(
        !msg.is_empty(),
        "Display must produce a user-facing message"
    );
}

/// The message text users see is unchanged by the extraction (spec ac-003):
/// each variant's `Display` is verbatim today's `GrowError` wording
/// (`src/grow/config.rs` Display impl).
#[test]
fn topology_error_display_text_is_verbatim_grow_error_text() {
    assert_eq!(
        TopologyError::NoBonds.to_string(),
        "the template frame carries no bonds; growth needs the bond graph — pack this \
         target with GenCanPack or supply connectivity",
        "GrowError::NoBonds wording must survive the move"
    );
    assert_eq!(
        TopologyError::NoAtomsBlock.to_string(),
        "the template frame has no readable atoms block",
        "GrowError::NoAtomsBlock wording must survive the move"
    );
    assert_eq!(
        TopologyError::BondOutOfRange { a: 4, b: 6, n: 6 }.to_string(),
        "bond (4, 6) references an atom outside the template (natoms = 6)",
        "GrowError::BondOutOfRange wording must survive the move"
    );
    assert_eq!(
        TopologyError::Disconnected.to_string(),
        "the template's bond graph does not connect all atoms",
        "GrowError::DisconnectedTemplate wording must survive the move"
    );
}

// ── frame_positions — the geometry reader that shares these errors ────────

/// `frame_positions` returns the template coordinates in atom order,
/// bit-for-bit (the inputs are synthesized, so exact equality is right).
#[test]
fn frame_positions_returns_template_coordinates() {
    let coords = zigzag_coords(5, 1.53);
    let frame = frame_from_parts(&coords, &chain_bonds(5));

    let got = frame_positions(&frame).expect("readable atoms block");
    assert_eq!(got, coords, "positions come back exactly as written");
}

/// `frame_positions` shares the atoms-block rejection with `from_frame`:
/// missing block, or a block without x/y/z, is `NoAtomsBlock`.
#[test]
fn frame_positions_without_coordinates_is_no_atoms_block() {
    let mut bare = Frame::new();
    bare.insert("bonds", bond_block(&chain_bonds(4)));
    let err = frame_positions(&bare).expect_err("no atoms block");
    assert!(
        matches!(err, TopologyError::NoAtomsBlock),
        "expected NoAtomsBlock, got {err:?}"
    );

    let partial = frame_missing_z_column(&zigzag_coords(4, 1.53));
    let err = frame_positions(&partial).expect_err("atoms block without z");
    assert!(
        matches!(err, TopologyError::NoAtomsBlock),
        "expected NoAtomsBlock, got {err:?}"
    );
}

// ── Regression scenario (spec Task 6 / ac-006) ────────────────────────────

/// Hard-coded topology goldens for the canonical 12-bead linear chain at
/// exclusion depth 3 (the all-atom default): the full root-inclusive
/// exclusion table and the CSR adjacency slices at both ends and in the
/// middle. Integer goldens — no tolerance, no third-party tool, no I/O.
///
/// Provenance: captured 2026-09-02 from
/// `InternalTree::from_frame_with_depth(frame, 3).exclusions(a)` on the
/// pre-refactor build (root-inclusive semantics).
#[test]
fn topology_regression_c12_chain_golden() {
    let topo = Topology::from_frame(&chain_frame(12, 1.53, true)).expect("linear 12-chain");

    let golden_exclusions: [&[u32]; 12] = [
        &[0, 1, 2, 3],
        &[0, 1, 2, 3, 4],
        &[0, 1, 2, 3, 4, 5],
        &[0, 1, 2, 3, 4, 5, 6],
        &[1, 2, 3, 4, 5, 6, 7],
        &[2, 3, 4, 5, 6, 7, 8],
        &[3, 4, 5, 6, 7, 8, 9],
        &[4, 5, 6, 7, 8, 9, 10],
        &[5, 6, 7, 8, 9, 10, 11],
        &[6, 7, 8, 9, 10, 11],
        &[7, 8, 9, 10, 11],
        &[8, 9, 10, 11],
    ];
    let got = topo.exclusions(3);
    assert_eq!(got.len(), golden_exclusions.len(), "one list per atom");
    for (atom, (list, golden)) in got.iter().zip(golden_exclusions.iter()).enumerate() {
        assert_eq!(
            list.as_slice(),
            *golden,
            "C12 depth-3 exclusions for atom {atom}"
        );
    }

    let golden_adjacency: [(usize, &[u32]); 4] =
        [(0, &[1]), (1, &[0, 2]), (5, &[4, 6]), (11, &[10])];
    for (atom, golden) in golden_adjacency {
        assert_eq!(
            topo.neighbors(atom),
            golden,
            "C12 CSR adjacency for atom {atom}, in bond-file insertion order"
        );
    }
}

// ── from_frame_with_positions — the full read its projections come from ───

/// `from_frame_with_positions` is the whole template read; `from_frame` is its
/// topology projection and `frame_positions` its geometry projection. On both
/// the linear and the branched fixture the pair is indistinguishable from the
/// two projections taken separately: coordinates bit-for-bit (synthesized
/// inputs, so exact equality is right), same atom count, same bond list in the
/// same order, same CSR slice for every atom.
#[test]
fn topology_from_frame_with_positions_matches_projections() {
    let (branched_coords, branched_bonds) = branched_parts();
    let cases = [
        ("12-bead linear chain", chain_frame(12, 1.53, true)),
        (
            "branched template",
            frame_from_parts(&branched_coords, &branched_bonds),
        ),
    ];

    for (what, frame) in &cases {
        let (topo, xyz) = Topology::from_frame_with_positions(frame)
            .unwrap_or_else(|e| panic!("{what} must read: {e}"));
        let projected =
            Topology::from_frame(frame).unwrap_or_else(|e| panic!("{what} via from_frame: {e}"));
        let positions =
            frame_positions(frame).unwrap_or_else(|e| panic!("{what} via frame_positions: {e}"));

        assert_eq!(
            xyz, positions,
            "{what}: the geometry half is exactly frame_positions"
        );
        assert_eq!(
            xyz.len(),
            topo.natoms(),
            "{what}: one coordinate triple per atom"
        );
        assert_eq!(
            topo.natoms(),
            projected.natoms(),
            "{what}: same atom count as from_frame"
        );
        assert_eq!(
            topo.bonds(),
            projected.bonds(),
            "{what}: same bonds, same bond-file order, as from_frame"
        );
        for atom in 0..topo.natoms() {
            assert_eq!(
                topo.neighbors(atom),
                projected.neighbors(atom),
                "{what}: CSR adjacency of atom {atom} must match from_frame"
            );
        }
    }
}

/// Both entry points refuse the same frames with the same error value — the
/// projection adds no check and drops none, and the payload survives it.
/// `TopologyError` is compared through its `Debug` rendering because the type
/// is deliberately not `PartialEq` (call sites match on the variant).
#[test]
fn topology_from_frame_with_positions_reports_same_errors() {
    let mut no_atoms = Frame::new();
    no_atoms.insert("bonds", bond_block(&chain_bonds(6)));

    let no_bonds = chain_frame(6, 1.53, false);

    let mut bad_bonds = chain_bonds(6);
    bad_bonds.push((4, 6)); // atom 6 does not exist in a 6-atom template
    let out_of_range = frame_from_parts(&zigzag_coords(6, 1.53), &bad_bonds);

    let err = Topology::from_frame_with_positions(&no_atoms).expect_err("no atoms block");
    assert!(
        matches!(err, TopologyError::NoAtomsBlock),
        "expected NoAtomsBlock, got {err:?}"
    );

    let err = Topology::from_frame_with_positions(&no_bonds).expect_err("no bonds block");
    assert!(
        matches!(err, TopologyError::NoBonds),
        "expected NoBonds, got {err:?}"
    );

    let err = Topology::from_frame_with_positions(&out_of_range).expect_err("out-of-range bond");
    match err {
        TopologyError::BondOutOfRange { a, b, n } => {
            assert_eq!(a, 4, "first endpoint");
            assert_eq!(b, 6, "offending endpoint");
            assert_eq!(n, 6, "template atom count");
        }
        other => panic!("expected BondOutOfRange, got {other:?}"),
    }

    for (what, frame) in [
        ("no atoms block", &no_atoms),
        ("no bonds block", &no_bonds),
        ("out-of-range bond", &out_of_range),
    ] {
        let paired = Topology::from_frame_with_positions(frame)
            .expect_err("from_frame_with_positions must refuse this frame");
        let projected = Topology::from_frame(frame).expect_err("from_frame must refuse this frame");
        assert_eq!(
            format!("{paired:?}"),
            format!("{projected:?}"),
            "{what}: both entry points must report the same error, payload included"
        );
    }
}

/// **Recorded precedence** (`.claude/notes/notes.md`, 2026-09-03): the leaf
/// never reports an atom-count error. A 2-atom frame with no bonds block is
/// `NoBonds` here, because `from_frame` reads the bond list in the same pass
/// that gives it the atom count. Growth's "at least 3 atoms" rule is applied
/// by `validate_template` / `InternalTree` *after* a successful read, so the
/// user-visible order is `NoAtomsBlock → NoBonds → TemplateTooSmall`; that
/// half belongs to `tests/grow.rs`, this assertion pins the leaf's half only.
#[test]
fn topology_error_precedence_no_bonds_before_too_small() {
    let err = Topology::from_frame(&chain_frame(2, 1.53, false))
        .expect_err("a 2-atom frame without a bonds block must be refused");

    assert!(
        matches!(err, TopologyError::NoBonds),
        "the leaf reports NoBonds, never an atom-count error: got {err:?}"
    );
}
