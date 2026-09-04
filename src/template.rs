//! Template-frame geometry: Cartesian coordinates from an `atoms` block.
//!
//! Crate-root leaf: `std` and molrs `Frame` only. Coordinates are in Å
//! (the crate's length unit, as carried by the frame — nothing is converted).
//! Bond graphs are `molrs::Topology`; this leaf does not read them.

use molrs::store::frame::Frame;
use molrs::types::F;

/// Why [`frame_positions`] cannot read a template's coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FramePositionsError {
    /// Missing `atoms`, or the block has no `x` / `y` / `z` float columns (Å).
    NoAtomsBlock,
}

/// Read a template frame's coordinates, in atom order, in Å.
///
/// The returned triples are coordinates in Å, one per atom, in the `atoms`
/// block's row order. Missing `atoms` or missing any of the `x` / `y` / `z`
/// float columns is [`FramePositionsError::NoAtomsBlock`]. No unit conversion.
pub(crate) fn frame_positions(frame: &Frame) -> Result<Vec<[F; 3]>, FramePositionsError> {
    let atoms = frame
        .get("atoms")
        .ok_or(FramePositionsError::NoAtomsBlock)?;
    let x = atoms
        .get_float("x")
        .ok_or(FramePositionsError::NoAtomsBlock)?;
    let y = atoms
        .get_float("y")
        .ok_or(FramePositionsError::NoAtomsBlock)?;
    let z = atoms
        .get_float("z")
        .ok_or(FramePositionsError::NoAtomsBlock)?;
    Ok(x.iter()
        .zip(y.iter())
        .zip(z.iter())
        .map(|((&a, &b), &c)| [a, b, c])
        .collect())
}

#[cfg(test)]
mod tests {
    use super::{FramePositionsError, frame_positions};
    use molrs::store::block::Block;
    use molrs::store::frame::Frame;
    use molrs::types::F;
    use ndarray::Array1;

    /// Planar zigzag bead-chain coordinates: tetrahedral (109.5°) bond angles
    /// in the x–z plane, all torsions trans. Lengths are in Å.
    fn zigzag_coords(n: usize, bond_len: F) -> Vec<[F; 3]> {
        let theta = 109.5 * std::f64::consts::PI as F / 180.0;
        let alpha = (std::f64::consts::PI as F - theta) / 2.0;
        let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());
        (0..n)
            .map(|i| [i as F * dx, 0.0, if i % 2 == 0 { 0.0 } else { dz }])
            .collect()
    }

    fn frame_from_coords(coords: &[[F; 3]]) -> Frame {
        let mut atoms = Block::new();
        for (name, k) in [("x", 0), ("y", 1), ("z", 2)] {
            let col: Vec<F> = coords.iter().map(|p| p[k]).collect();
            atoms
                .insert(name, Array1::from_vec(col).into_dyn())
                .expect("coordinate column");
        }
        let mut frame = Frame::new();
        frame.insert("atoms", atoms);
        frame
    }

    /// 4-bead zigzag, bond 1.53 Å — hard-coded x/y/z column goldens (Å).
    const ZIGZAG_4_BOND_153_A: [[F; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.249_461_579_397_369, 0.0, 0.883_032_140_756_967_3],
        [2.498_923_158_794_738, 0.0, 0.0],
        [3.748_384_738_192_107, 0.0, 0.883_032_140_756_967_3],
    ];

    fn assert_coords_match(got: &[[F; 3]], want: &[[F; 3]], what: &str) {
        assert_eq!(got.len(), want.len(), "{what}: atom count");
        for (i, (g, w)) in got.iter().zip(want).enumerate() {
            for k in 0..3 {
                let d = (g[k] - w[k]).abs();
                assert!(
                    d < 1e-12,
                    "{what}: atom {i}[{k}] got {} want {}, |Δ| = {d} Å (must be < 1e-12)",
                    g[k],
                    w[k]
                );
            }
        }
    }

    #[test]
    fn frame_positions_zigzag_round_trip() {
        let coords = zigzag_coords(4, 1.53);
        let frame = frame_from_coords(&coords);
        let got = frame_positions(&frame).expect("x/y/z present");
        assert_coords_match(&got, &coords, "written columns");
        assert_coords_match(&got, &ZIGZAG_4_BOND_153_A, "Å goldens");
    }

    #[test]
    fn frame_positions_missing_atoms_block() {
        let err = frame_positions(&Frame::new()).expect_err("no atoms block");
        assert!(
            matches!(err, FramePositionsError::NoAtomsBlock),
            "expected NoAtomsBlock, got {err:?}"
        );
    }

    #[test]
    fn frame_positions_missing_z_column() {
        let coords = zigzag_coords(4, 1.53);
        let mut atoms = Block::new();
        for (name, k) in [("x", 0), ("y", 1)] {
            let col: Vec<F> = coords.iter().map(|p| p[k]).collect();
            atoms
                .insert(name, Array1::from_vec(col).into_dyn())
                .expect("coordinate column");
        }
        let mut frame = Frame::new();
        frame.insert("atoms", atoms);
        let err = frame_positions(&frame).expect_err("atoms block without z");
        assert!(
            matches!(err, FramePositionsError::NoAtomsBlock),
            "expected NoAtomsBlock, got {err:?}"
        );
    }
}
