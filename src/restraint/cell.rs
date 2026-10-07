//! The primitive cell as a wall *and* as the lattice declaration — stated once.

use std::sync::Arc;

use molrs::core::Parallelepiped;
use molrs::core::SimBox;
use molrs::op::F;
use ndarray::array;

use super::{AtomRestraint, RegionRestraint};
use crate::PackError;

/// A cell from lengths (Å) and angles (degrees), origin at zero — the one
/// home of this conversion and its [`PackError::InvalidCell`] wording, shared
/// by [`CellRestraint::from_lengths_angles`] and the engines' cell builder.
pub(crate) fn simbox_from_lengths_angles(
    lengths: [F; 3],
    angles_deg: [F; 3],
    pbc: [bool; 3],
) -> Result<SimBox, PackError> {
    let h = SimBox::matrix_from_lengths_angles(lengths, angles_deg).map_err(|_| {
        PackError::InvalidCell {
            detail: format!("lengths {lengths:?} and angles {angles_deg:?} do not describe a cell"),
        }
    })?;
    SimBox::new(h, array![0.0, 0.0, 0.0], pbc).map_err(|_| singular_cell())
}

/// A cell from a lattice matrix whose columns are the lattice vectors (Å).
pub(crate) fn simbox_from_matrix(
    h: [[F; 3]; 3],
    origin: [F; 3],
    pbc: [bool; 3],
) -> Result<SimBox, PackError> {
    SimBox::from_matrix(h, origin, pbc).map_err(|_| singular_cell())
}

fn singular_cell() -> PackError {
    PackError::InvalidCell {
        detail: "lattice matrix is singular".to_string(),
    }
}

/// Confine every atom to a primitive cell and declare that cell as the
/// packing lattice.
///
/// The wall is a molrs [`Parallelepiped`] lifted through
/// [`RegionRestraint`], so the penalty is measured in Å perpendicular to the
/// bounding lattice planes and stays comparable with every other restraint
/// however tilted the cell is. The declaration is
/// [`AtomRestraint::declared_cell`]: the packer picks the lattice up from
/// the same object that confines the molecules, so a caller writes the cell
/// once and cannot declare a lattice it does not confine to.
///
/// Every axis carries a wall, periodic ones included. Under periodicity a
/// molecule at fractional `1.04` is the same configuration as one at `0.04`
/// and the pair kernel's minimum image cannot tell them apart, but a packer
/// has to emit coordinates someone can use: with nothing holding them,
/// molecules drift across hundreds of lattice images and every consumer then
/// has to wrap before the result means anything. The penalty is quadratic,
/// not a hard wall, so equilibrium leaves sub-tolerance excursions past a
/// face — the same behaviour as any region lift.
///
/// # Examples
///
/// ```
/// use molpack::{CellRestraint, Target};
/// # let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
///
/// // A hexagonal cell, periodic on every axis: one declaration.
/// let cell = CellRestraint::from_lengths_angles([26.0; 3], [90.0, 90.0, 120.0], [true; 3])?;
/// let target = Target::from_coords(pos, rad, 10).with_restraint(cell);
/// # Ok::<(), molpack::PackError>(())
/// ```
#[derive(Debug, Clone)]
pub struct CellRestraint {
    bx: SimBox,
    wall: RegionRestraint,
}

impl CellRestraint {
    /// Cell from lengths (Å) and angles (degrees), origin at zero.
    ///
    /// # Errors
    ///
    /// [`PackError::InvalidCell`] when the lengths and angles do not describe
    /// a cell.
    pub fn from_lengths_angles(
        lengths: [F; 3],
        angles_deg: [F; 3],
        pbc: [bool; 3],
    ) -> Result<Self, PackError> {
        Self::from_simbox(simbox_from_lengths_angles(lengths, angles_deg, pbc)?)
    }

    /// Cell from a lattice matrix whose **columns** are the lattice vectors (Å).
    ///
    /// # Errors
    ///
    /// [`PackError::InvalidCell`] when the matrix is singular.
    pub fn from_matrix(h: [[F; 3]; 3], origin: [F; 3], pbc: [bool; 3]) -> Result<Self, PackError> {
        Self::from_simbox(simbox_from_matrix(h, origin, pbc)?)
    }

    /// Cell from an existing [`SimBox`].
    ///
    /// # Errors
    ///
    /// [`PackError::InvalidCell`] when the box's lattice cannot bound a
    /// region (singular or non-finite).
    pub fn from_simbox(bx: SimBox) -> Result<Self, PackError> {
        let wall = Parallelepiped::new(bx.h_view().to_owned(), bx.origin_view().to_owned())
            .map_err(|detail| PackError::InvalidCell { detail })?;
        Ok(Self {
            bx,
            wall: RegionRestraint(Arc::new(wall)),
        })
    }

    /// The declared cell.
    pub fn cell(&self) -> &SimBox {
        &self.bx
    }
}

impl AtomRestraint for CellRestraint {
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F {
        self.wall.f(x, scale, scale2)
    }

    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        self.wall.fg(x, scale, scale2, g)
    }

    fn name(&self) -> &'static str {
        "CellRestraint"
    }

    fn declared_cell(&self) -> Option<SimBox> {
        Some(self.bx.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hexagonal() -> CellRestraint {
        CellRestraint::from_lengths_angles([26.0; 3], [90.0, 90.0, 120.0], [true; 3]).expect("cell")
    }

    /// Fractional 1.1 along `a` is 0.1 of the `a` plane spacing,
    /// `26 · sin 120°`, past the face — in Å, not in fractions.
    #[test]
    fn face_distance_is_perpendicular_angstrom() {
        let cell = hexagonal();
        let h = cell.cell().h_view();
        let frac = array![1.1, 0.5, 0.5];
        let x = h.dot(&frac);
        let d = 0.1 * 26.0 * 120.0_f64.to_radians().sin();
        let f = cell.f(&[x[0], x[1], x[2]], 1.0, 1.0);
        assert!((f - d * d).abs() < 1e-9, "f = {f}, expected {}", d * d);
    }

    #[test]
    fn gradient_is_the_unit_outward_normal_scaled_by_the_penalty() {
        let cell = hexagonal();
        let h = cell.cell().h_view();
        let frac = array![1.1, 0.5, 0.5];
        let x = [h.dot(&frac)[0], h.dot(&frac)[1], h.dot(&frac)[2]];
        let mut g = [0.0; 3];
        cell.fg(&x, 1.0, 1.0, &mut g);
        for k in 0..3 {
            let hh = 1e-6;
            let mut xp = x;
            xp[k] += hh;
            let mut xm = x;
            xm[k] -= hh;
            let fd = (cell.f(&xp, 1.0, 1.0) - cell.f(&xm, 1.0, 1.0)) / (2.0 * hh);
            assert!((g[k] - fd).abs() < 1e-6, "axis {k}: {} vs {fd}", g[k]);
        }
    }

    #[test]
    fn interior_is_free() {
        let cell = hexagonal();
        let h = cell.cell().h_view();
        let x = h.dot(&array![0.3, 0.6, 0.5]);
        assert_eq!(cell.f(&[x[0], x[1], x[2]], 1.0, 1.0), 0.0);
    }

    #[test]
    fn declared_cell_round_trips() {
        let cell = CellRestraint::from_matrix(
            [[10.0, 2.0, 0.0], [0.0, 9.0, 1.0], [0.0, 0.0, 8.0]],
            [1.0, 2.0, 3.0],
            [true, true, false],
        )
        .expect("cell");
        let bx = cell.declared_cell().expect("declares its cell");
        assert_eq!(
            bx.h_view(),
            array![[10.0, 2.0, 0.0], [0.0, 9.0, 1.0], [0.0, 0.0, 8.0]]
        );
        assert_eq!(bx.origin_view(), array![1.0, 2.0, 3.0]);
        assert_eq!(bx.pbc(), [true, true, false]);
        assert_eq!(cell.name(), "CellRestraint");
    }

    #[test]
    fn degenerate_cell_is_named() {
        assert!(matches!(
            CellRestraint::from_lengths_angles([10.0; 3], [170.0, 170.0, 170.0], [true; 3]),
            Err(PackError::InvalidCell { .. })
        ));
        assert!(matches!(
            CellRestraint::from_matrix(
                [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
                [true; 3]
            ),
            Err(PackError::InvalidCell { .. })
        ));
    }
}
