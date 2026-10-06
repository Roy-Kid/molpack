//! Coarse-grain geometry, synthesized in process — no data files, no `io`.
//!
//! Every species is built as a `molrs::Frame` (an `atoms` block, plus a `bonds`
//! block for the chain) so targets come from [`Target::new`] and carry their
//! topology into the packed result.

use molpack::F;
use molrs::store::block::Block;
use molrs::store::frame::Frame;
use molrs::system::atomistic::Atomistic;
use ndarray::Array1;

/// Bead label for a chain backbone site.
pub const BACKBONE: &str = "C";
/// Bead label for a chain's surface-seeking (sticky) site.
pub const STICKY: &str = "N";
/// Bead label for a substrate site.
pub const SUBSTRATE: &str = "Au";
/// Bead label for a solvent bead.
pub const SOLVENT: &str = "Ar";

fn atoms_block(elements: Vec<String>, xs: Vec<F>, ys: Vec<F>, zs: Vec<F>) -> Block {
    let mut atoms = Block::new();
    atoms
        .insert("element", Array1::from_vec(elements).into_dyn())
        .expect("element column");
    atoms
        .insert("x", Array1::from_vec(xs).into_dyn())
        .expect("x column");
    atoms
        .insert("y", Array1::from_vec(ys).into_dyn())
        .expect("y column");
    atoms
        .insert("z", Array1::from_vec(zs).into_dyn())
        .expect("z column");
    atoms
}

fn frame_of(atoms: Block) -> Frame {
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame
}

/// A periodically functionalized coarse-grain chain.
///
/// Returns the frame, its graph (for rotatable-bond perception), and the
/// 0-based indices of the sticky beads.
///
/// # Why a zigzag
///
/// The beads are laid out with a tetrahedral valence angle rather than on a
/// straight line. This is load-bearing: in a collinear chain every rotatable
/// bond's axis passes through all of its downstream beads, so torsion rotation
/// is the identity and the in-loop optimizer becomes a silent no-op — it
/// reports accepted moves while the geometry never changes.
///
/// # Bond classes
///
/// The `bonds` block states connectivity only. molrs reads that faithfully, so
/// every bond arrives as `BondType::Unknown`; `TorsionMcOptimizer` supplies the
/// single-bond fallback that rotatable-bond perception needs. Nothing here has
/// to know about `bond_type`.
pub fn chain(n_beads: usize, bond_len: F, sticky_every: usize) -> (Frame, Atomistic, Vec<usize>) {
    assert!(n_beads >= 4, "a chain needs at least one rotatable bond");
    assert!(sticky_every >= 1, "sticky spacing must be positive");

    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());

    let sticky: Vec<usize> = (0..n_beads).step_by(sticky_every).collect();
    let elements: Vec<String> = (0..n_beads)
        .map(|i| {
            if sticky.contains(&i) {
                STICKY.to_string()
            } else {
                BACKBONE.to_string()
            }
        })
        .collect();

    let xs: Vec<F> = (0..n_beads).map(|i| i as F * dx).collect();
    let ys: Vec<F> = vec![0.0; n_beads];
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

    let mut frame = frame_of(atoms_block(elements, xs, ys, zs));
    frame.insert("bonds", bonds);

    let graph = Atomistic::from_frame(&frame).expect("chain graph");
    (frame, graph, sticky)
}

/// A flat square lattice of substrate beads centred on the xy origin.
///
/// The substrate is packed as a single fixed target, so it contributes real
/// excluded volume through the pair term — chains cannot sink through it. A
/// plane restraint alone would not do that: it only acts on the atoms it is
/// attached to.
pub fn substrate(nx: usize, ny: usize, spacing: F, z: F) -> Frame {
    let n = nx * ny;
    let (x0, y0) = (
        -0.5 * (nx as F - 1.0) * spacing,
        -0.5 * (ny as F - 1.0) * spacing,
    );
    let mut xs = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    for ix in 0..nx {
        for iy in 0..ny {
            xs.push(x0 + ix as F * spacing);
            ys.push(y0 + iy as F * spacing);
        }
    }
    frame_of(atoms_block(
        vec![SUBSTRATE.to_string(); n],
        xs,
        ys,
        vec![z; n],
    ))
}

/// A single solvent bead.
pub fn solvent_bead() -> Frame {
    frame_of(atoms_block(
        vec![SOLVENT.to_string()],
        vec![0.0],
        vec![0.0],
        vec![0.0],
    ))
}
