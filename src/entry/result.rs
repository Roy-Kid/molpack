//! The engine-run outcome: [`PackResult`], the placement solution it
//! carries, and the coordinate reordering that assembles the frame in
//! target-declared order.

use molrs::spatial::simbox::SimBox;
use molrs::types::F;

use crate::target::Target;

/// The solver-native placement solution for the FREE copies, captured
/// verbatim at the end of a run (placement-seeding spec): the packed
/// (COM | Euler) vector, the per-copy centered reference conformers in
/// xcart order, a per-copy atom-count fingerprint for validation, and the
/// simbox the run installed. A later entry continues on this state with
/// zero conversion — reconstructing from the assembled frame would recompute
/// COMs and lose bitwise continuity ((p − com) + com ≠ p).
#[derive(Debug, Clone)]
pub(crate) struct Placements {
    /// `6 * n_free_mol`: COM block then Euler block (the solver `x`).
    pub(crate) x: Vec<F>,
    /// Per-copy centered reference coordinates for the free atoms.
    pub(crate) coor: Vec<[F; 3]>,
    /// Atoms per free copy, in xcart (declared) order.
    pub(crate) copy_atoms: Vec<usize>,
    /// The simbox this run installed (grid + minimum image).
    pub(crate) cell: SimBox,
}

/// The outcome of one engine run: the packed frame plus the shared
/// objective's verdict.
///
/// The `frame` contains an "atoms" block with x, y, z, element, mol_id
/// columns — moved from the packing context (zero-copy ownership transfer).
#[derive(Debug, Clone)]
pub struct PackResult {
    /// Atoms frame with x, y, z (f64), element (String), mol_id (i64).
    pub frame: molrs::Frame,
    /// The verbatim placement solution, for cross-entry seeding
    /// (`GenCanPack::seeded_from`).
    pub(crate) placements: Placements,
    /// Maximum inter-molecular distance violation at termination.
    pub fdist: F,
    /// Maximum constraint violation at termination.
    pub frest: F,
    /// Whether the packing converged (`fdist < precision && frest < precision`).
    pub converged: bool,
    /// How many times the growth solver had to relax its constructive
    /// hard-core guarantee. Always `0` on the GENCAN path; a grown structure
    /// is only `converged` when it is `0` there too.
    pub softened: usize,
}

impl PackResult {
    /// Extract atom positions as `Vec<[F; 3]>` (SoA→AoS conversion).
    pub fn positions(&self) -> Vec<[F; 3]> {
        let atoms = self.frame.get("atoms").expect("frame has no 'atoms' block");
        let x = atoms.get_float("x").expect("no 'x' column");
        let y = atoms.get_float("y").expect("no 'y' column");
        let z = atoms.get_float("z").expect("no 'z' column");
        x.iter()
            .zip(y.iter())
            .zip(z.iter())
            .map(|((&xi, &yi), &zi)| [xi, yi, zi])
            .collect()
    }

    /// Number of atoms in the result.
    #[inline]
    pub fn natoms(&self) -> usize {
        self.frame.get("atoms").and_then(|b| b.nrows()).unwrap_or(0)
    }
}

/// Reorder packed coordinates from the packer's internal `xcart` layout
/// (all free targets first, then all fixed targets) into the target-declared
/// order that [`crate::assemble::assemble_frame`] expects (target-by-target,
/// copy-by-copy, atom-by-atom).
///
/// Without this, a fixed target declared *before* a free target — e.g. a
/// `fixed` protein in a solvation box — has its topology replayed onto the
/// free atoms' coordinates, scrambling the rigid molecule. A no-op when no
/// target is fixed (the free blocks are already in declared order).
pub(crate) fn positions_in_target_order(
    targets: &[Target],
    xcart: &[[F; 3]],
    n_free_atoms: usize,
) -> Vec<[F; 3]> {
    let mut out = Vec::with_capacity(xcart.len());
    let mut free_cursor = 0usize;
    let mut fixed_cursor = n_free_atoms;
    for t in targets {
        if t.fixed_at.is_some() {
            let n = t.natoms();
            out.extend_from_slice(&xcart[fixed_cursor..fixed_cursor + n]);
            fixed_cursor += n;
        } else {
            let n = t.count * t.natoms();
            out.extend_from_slice(&xcart[free_cursor..free_cursor + n]);
            free_cursor += n;
        }
    }
    out
}
