//! Coarse-grain geometry for the translocation system, synthesized in process.

use molrs::op::types::F;
use molrs::store::Block;
use molrs::store::Frame;
use molrs::system::Atomistic;
use ndarray::Array1;

pub const MEMBRANE: &str = "Au";
pub const CIS_TAG: &str = "O";
pub const PORE_TAG: &str = "N";
pub const TRANS_TAG: &str = "S";
pub const BACKBONE: &str = "C";
pub const SOLVENT: &str = "Ar";

fn atoms_block(elements: Vec<String>, xs: Vec<F>, ys: Vec<F>, zs: Vec<F>) -> Block {
    let mut atoms = Block::new();
    for (k, v) in [("x", xs), ("y", ys), ("z", zs)] {
        atoms
            .insert(k, Array1::from_vec(v).into_dyn())
            .expect("coordinate column");
    }
    atoms
        .insert("element", Array1::from_vec(elements).into_dyn())
        .expect("element column");
    atoms
}

fn frame_of(atoms: Block) -> Frame {
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame
}

/// Which segment of the chain a bead belongs to. The three segments carry three
/// *different* restraint geometries on the same molecule — that is the point of
/// the example.
pub struct Segments {
    pub cis: Vec<usize>,
    pub pore: Vec<usize>,
    pub trans: Vec<usize>,
}

/// A linear coarse-grain chain, zigzag so its torsions have somewhere to go.
///
/// Beads are tagged by segment so the packed structure can be inspected in a
/// viewer, and so the three per-atom restraints have index sets to attach to.
pub fn chain(
    n_beads: usize,
    bond_len: F,
    cis_len: usize,
    pore_len: usize,
) -> (Frame, Atomistic, Segments) {
    assert!(
        cis_len + pore_len < n_beads,
        "cis and pore segments must leave room for a trans segment"
    );

    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());

    // cis | spacer | pore | spacer | trans, with the pore segment centred.
    let pore_start = (n_beads - pore_len) / 2;
    let segments = Segments {
        cis: (0..cis_len).collect(),
        pore: (pore_start..pore_start + pore_len).collect(),
        trans: (n_beads - cis_len..n_beads).collect(),
    };

    let elements: Vec<String> = (0..n_beads)
        .map(|i| {
            if segments.cis.contains(&i) {
                CIS_TAG.to_string()
            } else if segments.pore.contains(&i) {
                PORE_TAG.to_string()
            } else if segments.trans.contains(&i) {
                TRANS_TAG.to_string()
            } else {
                BACKBONE.to_string()
            }
        })
        .collect();

    let xs: Vec<F> = (0..n_beads).map(|i| i as F * dx).collect();
    let zs: Vec<F> = (0..n_beads)
        .map(|i| if i % 2 == 0 { 0.0 } else { dz })
        .collect();

    let mut bonds = Block::new();
    bonds
        .insert(
            "atomi",
            Array1::from_vec((0..n_beads as u32 - 1).collect::<Vec<u32>>()).into_dyn(),
        )
        .expect("atomi column");
    bonds
        .insert(
            "atomj",
            Array1::from_vec((1..n_beads as u32).collect::<Vec<u32>>()).into_dyn(),
        )
        .expect("atomj column");

    let mut frame = frame_of(atoms_block(elements, xs, vec![0.0; n_beads], zs));
    frame.insert("bonds", bonds);
    let graph = Atomistic::from_frame(&frame).expect("chain graph");
    (frame, graph, segments)
}

/// A rigid membrane: `layers` square lattice sheets with a circular hole
/// punched at each of `holes`. The rim is real excluded volume, so the only way
/// through a pore is through it — which is exactly what makes the threaded
/// state unreachable by dynamics that did not start there.
pub fn membrane(n: usize, spacing: F, layers: &[F], holes: &[[F; 2]], hole_radius: F) -> Frame {
    let origin = -0.5 * (n as F - 1.0) * spacing;
    let (mut xs, mut ys, mut zs) = (Vec::new(), Vec::new(), Vec::new());
    for &z in layers {
        for ix in 0..n {
            for iy in 0..n {
                let (x, y) = (origin + ix as F * spacing, origin + iy as F * spacing);
                let in_hole = holes.iter().any(|h| {
                    let (dx, dy) = (x - h[0], y - h[1]);
                    (dx * dx + dy * dy).sqrt() < hole_radius
                });
                if in_hole {
                    continue;
                }
                xs.push(x);
                ys.push(y);
                zs.push(z);
            }
        }
    }
    let count = xs.len();
    frame_of(atoms_block(vec![MEMBRANE.to_string(); count], xs, ys, zs))
}

pub fn solvent_bead() -> Frame {
    frame_of(atoms_block(
        vec![SOLVENT.to_string()],
        vec![0.0],
        vec![0.0],
        vec![0.0],
    ))
}
