//! Per-copy centroid reduction: a species' atoms collapse to one site per
//! molecule, and a gradient on those sites scatters back onto the atoms.
//!
//! The packer hands a collective restraint one flat, **copy-major** coordinate
//! slice per species, so copy `c` owns `coords[c·m .. (c+1)·m]` for
//! `m = GroupCtx::natoms_per_copy`. That is the only fact this module needs; it
//! is what lets a term reason about *molecules* instead of atoms.
//!
//! A copy's atoms are always contiguous in space, so the centroid needs no
//! unwrapping: the packer builds `xcart` by rotating a reference conformer about
//! the molecule's own centre and translating it, and never folds the result back
//! into the cell. A molecule therefore cannot straddle a periodic boundary with
//! half its atoms an image away, which is the case a naive mean would get wrong.
//! (Distances *between* two centroids are a different matter — those are the
//! caller's to take under the minimum image.)
//!
//! The centroid is **geometric**, not mass-weighted. Two reasons, in order:
//! solvers on this seam consume geometry only, and a mass-weighted centre would
//! make the packed result depend on which force field labelled the input; and
//! masses are genuinely absent for coarse-grained beads and for targets built
//! from bare coordinates.

use molrs::op::superpose::centroid;
use molrs::op::types::F;

/// Geometric centroid of every copy: `R_c = (1/m) Σ_{i∈c} r_i`, molrs's
/// [`centroid`] at unit weights.
///
/// Returns one site per copy, in copy order. Empty when `m == 0`.
pub(super) fn centroids(coords: &[[F; 3]], m: usize) -> Vec<[F; 3]> {
    if m == 0 {
        return Vec::new();
    }
    let unit_weights = vec![1.0; m];
    coords
        .chunks_exact(m)
        .map(|copy| centroid(copy, &unit_weights).expect("m unit weights sum to m > 0"))
        .collect()
}

/// Scatter a per-site gradient back onto the atoms that formed it.
///
/// `∂R_c/∂r_i = 1/m` for every atom of copy `c`, so each atom takes an equal
/// share. Accumulates with `+=`, per the [`Restraint`](super::Restraint)
/// gradient convention.
pub(super) fn scatter(dsites: &[[F; 3]], m: usize, grads: &mut [[F; 3]]) {
    if m == 0 {
        return;
    }
    let inv_m = 1.0 / m as F;
    for (copy, ds) in grads.chunks_exact_mut(m).zip(dsites.iter()) {
        let share = [ds[0] * inv_m, ds[1] * inv_m, ds[2] * inv_m];
        for g in copy.iter_mut() {
            g[0] += share[0];
            g[1] += share[1];
            g[2] += share[2];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn centroid_of_one_atom_per_copy_is_the_atom() {
        let coords = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        assert_eq!(centroids(&coords, 1), coords.to_vec());
    }

    #[test]
    fn centroid_averages_within_a_copy_not_across_copies() {
        // Two copies of a 2-atom molecule; centroids must be (0.5,0,0) and
        // (10.5,0,0), never the global mean.
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [11.0, 0.0, 0.0],
        ];
        let sites = centroids(&coords, 2);
        assert_eq!(sites, vec![[0.5, 0.0, 0.0], [10.5, 0.0, 0.0]]);
    }

    #[test]
    fn scatter_splits_a_site_gradient_evenly_and_accumulates() {
        let mut grads = [[1.0, 0.0, 0.0]; 4]; // pre-loaded, must accumulate
        scatter(&[[4.0, 0.0, 0.0], [8.0, 0.0, 0.0]], 2, &mut grads);
        assert_eq!(
            grads,
            [
                [3.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [5.0, 0.0, 0.0]
            ]
        );
    }

    #[test]
    fn zero_atoms_per_copy_is_inert() {
        assert!(centroids(&[[0.0; 3]], 0).is_empty());
        let mut grads = [[7.0; 3]];
        scatter(&[[1.0; 3]], 0, &mut grads);
        assert_eq!(grads, [[7.0; 3]]);
    }
}
