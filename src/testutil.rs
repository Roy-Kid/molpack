//! Test fixtures shared by the in-module unit tests (compiled only under
//! `cfg(test)`). Everything is deterministic: no seeds, files or clocks.

use std::sync::Arc;

use molrs::core::Atom;
use molrs::core::Atomistic;
use molrs::core::Block;
use molrs::core::Cuboid;
use molrs::core::Frame;
use molrs::op::F;
use ndarray::{Array1, Array2, array};

use crate::RegionRestraint;

/// "Stay inside the axis-aligned box `[min, max]`": a molrs [`Cuboid`] lifted
/// by [`RegionRestraint`].
pub(crate) fn inside_box(min: [F; 3], max: [F; 3]) -> RegionRestraint {
    RegionRestraint(Arc::new(Cuboid::new(
        array![min[0], min[1], min[2]],
        array![max[0] - min[0], max[1] - min[1], max[2] - min[2]],
    )))
}

/// Planar zigzag bead-chain coordinates: tetrahedral (109.5°) bond angles in
/// the x–z plane, all torsions trans. The zigzag is load-bearing: a collinear
/// chain makes every torsion a no-op.
pub(crate) fn zigzag_coords(n: usize, bond_len: F) -> Vec<[F; 3]> {
    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());
    (0..n)
        .map(|i| [i as F * dx, 0.0, if i % 2 == 0 { 0.0 } else { dz }])
        .collect()
}

/// The `(i, i+1)` bond list of a linear chain.
pub(crate) fn chain_bonds(n: usize) -> Vec<(u32, u32)> {
    (0..n.saturating_sub(1) as u32)
        .map(|i| (i, i + 1))
        .collect()
}

/// Coordinates + explicit bond list as a `molrs::core::Frame`: an atoms block with
/// `x`/`y`/`z` and, unless `bonds` is empty, a bonds block with `atomi`/`atomj`.
/// Deliberately no `bond_type` column — the shape a PDB CONECT list or a
/// hand-built coarse-grain frame has.
pub(crate) fn frame_from_parts(coords: &[[F; 3]], bonds: &[(u32, u32)]) -> Frame {
    let mut frame = Frame::new();
    frame
        .set_coords(Array2::from(coords.to_vec()).view())
        .expect("an N x 3 array always fits a fresh frame");
    if !bonds.is_empty() {
        let mut block = Block::new();
        let ai: Vec<u32> = bonds.iter().map(|&(i, _)| i).collect();
        let aj: Vec<u32> = bonds.iter().map(|&(_, j)| j).collect();
        block
            .insert("atomi", Array1::from_vec(ai).into_dyn())
            .expect("atomi column");
        block
            .insert("atomj", Array1::from_vec(aj).into_dyn())
            .expect("atomj column");
        frame.insert("bonds", block);
    }
    frame
}

/// A bonded zigzag bead chain as a `molrs::core::Frame`.
pub(crate) fn chain_frame(n: usize, bond_len: F) -> Frame {
    frame_from_parts(&zigzag_coords(n, bond_len), &chain_bonds(n))
}

/// A linear chain's bond graph: `n` bare atoms bonded `(i, i+1)`.
pub(crate) fn chain_graph(n: usize) -> Atomistic {
    let mut g = Atomistic::new();
    let ids: Vec<_> = (0..n).map(|_| g.add_atom(Atom::new())).collect();
    for pair in ids.windows(2) {
        g.add_bond(pair[0], pair[1]).expect("bond");
    }
    g
}
