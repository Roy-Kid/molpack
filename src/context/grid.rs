//! Cell-grid installation shared by every stage.
//!
//! The grid is stage infrastructure, not a GENCAN step: growth, lattice
//! growth, and a continuation that skips initial placement all install the
//! same box and the same coverage. Living here keeps those stages from
//! depending on the rigid-body driver.

use molrs::core::CellGrid;
use molrs::core::SimBox;
use molrs::op::F;

use crate::context::{NONE_IDX, PackSystem};

/// Install the resolved simulation box and cell grid on the system and bin
/// the fixed atoms.
///
/// A solver that skips initial placement still needs a populated grid for
/// the shared-objective evaluation.
pub(crate) fn install_simbox_and_grid(
    sys: &mut PackSystem,
    simbox: SimBox,
    radmax: F,
    discale: F,
    free_atoms: usize,
) {
    sys.simbox = simbox;
    let periodic = sys.simbox.pbc();

    let cell_side = if radmax > 0.0 {
        discale * 1.01 * radmax
    } else {
        1.0
    };
    log::debug!("setting up cell grid (cell_side={cell_side:.4})");
    // Cap the total cell count. With no spatial constraint the fallback box is
    // ±`sidemax` (default 1000 Å) wide, which drives an uncapped grid to ~10⁹
    // cells and OOMs `resize_cell_arrays`. There is no benefit to having far
    // more cells than atoms, so the budget scales with `ntotat` under a hard
    // ceiling. Coarser cells only slow the neighbor search — they never change
    // the packing result. `for_cutoff_capped` coarsens by a common factor and
    // pins a thin axis, so a floor of one cell on a short axis cannot blow
    // the budget.
    let max_total_cells = sys.ntotat.max(1).saturating_mul(64).clamp(1 << 16, 1 << 22);
    sys.grid = CellGrid::for_cutoff_capped(&sys.simbox, cell_side, max_total_cells);
    log::debug!("celldim={:?}  periodic={:?}", sys.grid.celldim(), periodic);

    sys.resize_cell_arrays();

    // Add fixed atoms to latomfix (Packmol lines 303-318)
    for icart in free_atoms..sys.ntotat {
        let pos = sys.xcart[icart];
        let icell = sys.grid.cell_of(&sys.simbox, pos);
        if sys.latomfix[icell] == NONE_IDX {
            sys.fixed_cells.push(icell);
        }
        sys.latomnext[icart] = sys.latomfix[icell];
        sys.latomfix[icell] = icart as u32;
    }
}

/// Derive `radmax` and the free-atom count from the system, then install
/// the resolved cell and its grid — the shared "box and its cell grid"
/// prelude for every stage that hands `run` an already-resolved [`SimBox`]
/// (growth, lattice growth, and a GENCAN stage that continues from existing
/// placements).
///
/// The coverage scale is [`coverage_radmax`], the same derivation initial
/// placement uses — there is one answer to "how wide must a cell be", and
/// both entries into the grid read it from the same place.
pub(crate) fn install_resolved_cell(sys: &mut PackSystem, cell: &SimBox, discale: F) {
    let radmax = coverage_radmax(sys);
    let free_atoms = sys.ntotat - sys.nfixedat;
    install_simbox_and_grid(sys, cell.clone(), radmax, discale, free_atoms);
}

/// The distance the cell grid has to cover: the largest **diameter** among
/// the unscaled packing radii (Packmol's `radmax`, `packmol.f90` 532-534).
///
/// Two properties, both load-bearing:
///
/// - **Diameter, not radius.** The pair kernel interacts out to
///   `radius_i + radius_j` on radii already multiplied by `discale`, so the
///   reach is `2 · discale · max(radius_ini)` and the cell side
///   (`discale · 1.01 · radmax`) covers it with 1% to spare. Sized from the
///   radius instead, the side is half the reach and a `±1` stencil never
///   enumerates the pairs in between — the objective then reports a clean
///   structure that overlaps.
/// - **`radius_ini`, not `radius`.** `radius` is GENCAN's transient working
///   copy, rescaled at every phase start, so reading it would size the grid
///   from whatever the previous stage happened to leave behind.
pub(crate) fn coverage_radmax(sys: &PackSystem) -> F {
    sys.radius_ini
        .iter()
        .copied()
        .map(|r| 2.0 * r)
        .fold(0.0 as F, F::max)
}

#[cfg(test)]
mod tests {
    //! The cell grid must cover the pair kernel's reach.
    //!
    //! A `±1` stencil finds every pair closer than one cell side, so the side
    //! has to be at least the largest interacting distance — `2 · discale ·
    //! max(radius_ini)`, since the kernel's cutoff is `radius_i + radius_j` on
    //! radii already scaled by `discale`. A grid sized from the *radius*
    //! instead of the *diameter* is half that, and the pairs in between are
    //! not merely found late: they are never enumerated, so the objective
    //! reports a clean structure that overlaps.

    use super::install_resolved_cell;
    use crate::PackSystem;
    use crate::objective::compute_f;
    use molrs::core::SimBox;
    use molrs::op::F;
    use molrs::op::F3;

    const DISCALE: F = 1.1;

    /// Two single-atom molecules `dx` apart on the x axis in a 20 Å free box,
    /// with the working radius already scaled by `discale` (what a phase
    /// start leaves behind).
    fn two_atoms(dx: F) -> (PackSystem, Vec<F>) {
        let mut sys = PackSystem::new(2, 2, 1);
        sys.ntype_with_fixed = 1;
        sys.nmols = vec![2];
        sys.natoms = vec![1];
        sys.idfirst = vec![0];
        sys.comptype = vec![true];
        sys.coor = vec![[0.0; 3]; 2];
        sys.radius_ini = vec![1.0; 2];
        sys.radius = vec![DISCALE; 2];
        sys.fscale = vec![1.0; 2];
        sys.ibmol = vec![0, 1];
        sys.iratom_offsets = vec![0, 0, 0];
        sys.sync_atom_props();

        let cell = SimBox::cube(20.0, F3::zeros(3), [false; 3]).expect("box");
        install_resolved_cell(&mut sys, &cell, DISCALE);

        // Placed off the cell boundary so the pair straddles two cell widths
        // under the under-sized grid.
        let x = vec![
            1.1,
            10.0,
            10.0,
            1.1 + dx,
            10.0,
            10.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        (sys, x)
    }

    #[test]
    fn an_overlapping_pair_two_cells_apart_is_still_seen() {
        // Contact is `radius_i + radius_j` = 2.2 Å; 2.15 Å is an overlap, and
        // it is further apart than the radius-sized cell (1.11 Å), so only a
        // diameter-sized grid enumerates it.
        let (mut sys, x) = two_atoms(2.15);
        assert!(
            compute_f(&x, &mut sys) > 0.0,
            "a 2.15 Å pair inside a 2.2 Å contact must be found: the grid has \
             to cover the kernel's reach, not half of it",
        );
    }

    #[test]
    fn a_pair_beyond_contact_stays_free() {
        // The complement: coverage is not an excuse to invent a penalty.
        let (mut sys, x) = two_atoms(2.25);
        assert_eq!(
            compute_f(&x, &mut sys),
            0.0,
            "2.25 Å is outside the 2.2 Å contact — no pair term is owed",
        );
    }
}
