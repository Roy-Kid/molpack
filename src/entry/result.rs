//! The engine-run outcome: [`PackResult`], the placement solution it
//! carries, the intra-molecular residual [`IntraResidual`], and the
//! coordinate reordering that assembles the frame in target-declared order.

use molrs::spatial::simbox::SimBox;
use molrs::types::F;

use crate::context::RigidView;
use crate::target::Target;

/// The solver-native placement solution for the FREE copies, captured
/// verbatim at the end of a run: the run's [`RigidView`] (the packed
/// COM + Euler placement vector), the per-copy centered reference conformers
/// in xcart order, a per-copy atom-count fingerprint for validation, and the
/// simbox the run installed. A later entry continues on this state with
/// zero conversion — reconstructing from the assembled frame would recompute
/// COMs and lose bitwise continuity ((p − com) + com ≠ p).
#[derive(Debug, Clone)]
pub(crate) struct Placements {
    /// The run's rigid degrees of freedom, as the solver left them.
    pub(crate) rigid: RigidView,
    /// Per-copy centered reference coordinates for the free atoms.
    pub(crate) coor: Vec<[F; 3]>,
    /// Atoms per free copy, in xcart (declared) order.
    pub(crate) copy_atoms: Vec<usize>,
    /// The simbox this run installed (grid + minimum image).
    pub(crate) cell: SimBox,
}

/// Intra-molecular residual of one configuration: the minimum same-copy
/// pair distance among **scored** pairs and among **exempted** pairs, in Å
/// (minimum image).
///
/// An empty class is `+∞` ([`F::INFINITY`]). There is no [`Default`] — a
/// residual is classified from coordinates and a skip table, not from a
/// placeholder. Distances are Euclidean lengths of
/// [`SimBox::shortest_vector_impl`].
///
/// Classification is driven by each target's
/// [`Target::special_bonds`]. This type does not pick an exclusion depth.
/// Analogous per-target table: [`crate::grow::internal::InternalTree::from_frame`].
#[derive(Debug, Clone, Copy)]
pub struct IntraResidual {
    /// Minimum same-copy pair distance among pairs the table scores (Å).
    pub scored: F,
    /// Minimum same-copy pair distance among pairs the table exempts (Å).
    pub exempted: F,
}

/// How one target's bond graph classifies same-copy pairs.
enum IntraSkip {
    /// No template, or a zero-edge graph: identity exemption only.
    Identity,
    /// `Topology::from_frame` `NotFound` / `Validation`: omit the target.
    Omit,
    /// Per-atom skip lists from [`molrs::Topology::exclusions`] (sorted,
    /// root-inclusive).
    Partners(Vec<Vec<usize>>),
}

impl IntraResidual {
    /// Classify every same-copy `i < j` pair of `targets` against each
    /// target's [`Target::special_bonds`] and return the two class minima
    /// (Å, minimum image).
    ///
    /// Spec 05 reads `Target.special_bonds` (no longer broadcasts depth 3).
    ///
    /// A missing template or a zero-edge bond graph is identity exemption
    /// only (`i == j` skipped; every `i != j` scored). [`molrs::Topology::from_frame`]
    /// errors `NotFound` and `Validation` omit that target — neither class is
    /// updated — so a broken 1-2 is not reported as scored.
    pub fn from_targets(targets: &[Target], positions: &[[F; 3]], simbox: &SimBox) -> Self {
        let mut scored = F::INFINITY;
        let mut exempted = F::INFINITY;
        let mut offset = 0usize;
        for target in targets {
            let table = &target.special_bonds;
            let n = target.natoms();
            let ncopy = if target.fixed_at.is_some() {
                1
            } else {
                target.count
            };
            let span = ncopy * n;
            let skip = match target.template.as_ref() {
                None => IntraSkip::Identity,
                Some(frame) => match molrs::Topology::from_frame(frame) {
                    Ok(topo) if topo.n_bonds() == 0 => IntraSkip::Identity,
                    Ok(topo) => IntraSkip::Partners(topo.exclusions(table)),
                    Err(_) => IntraSkip::Omit,
                },
            };
            if matches!(skip, IntraSkip::Omit) {
                offset += span;
                continue;
            }
            let tail_scored = table.weight(usize::MAX) != 0.0;
            for c in 0..ncopy {
                let base = offset + c * n;
                for i in 0..n {
                    for j in (i + 1)..n {
                        let dr =
                            simbox.shortest_vector_impl(positions[base + i], positions[base + j]);
                        let dist = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                        let pair_exempted = match &skip {
                            IntraSkip::Identity => false,
                            IntraSkip::Partners(lists) => {
                                lists[i].binary_search(&j).is_ok() || !tail_scored
                            }
                            IntraSkip::Omit => false,
                        };
                        if pair_exempted {
                            exempted = exempted.min(dist);
                        } else {
                            scored = scored.min(dist);
                        }
                    }
                }
            }
            offset += span;
        }
        Self { scored, exempted }
    }
}

/// The outcome of one engine run: the packed frame plus the shared
/// objective's verdict and the intra-molecular residual.
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
    /// Intra-molecular residual (Å, minimum image): scored vs exempted
    /// same-copy minima. `Pipeline::assemble` reads each target's
    /// [`crate::target::Target::special_bonds`].
    pub intra: IntraResidual,
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

#[cfg(test)]
mod tests {
    use super::IntraResidual;
    use crate::target::Target;
    use molrs::BondDistanceWeights;
    use molrs::spatial::simbox::SimBox;
    use molrs::store::block::Block;
    use molrs::store::frame::Frame;
    use molrs::types::{F, Idx};
    use ndarray::{Array1, array};

    /// Coordinates + explicit bond list as a `molrs::Frame`: atoms `x`/`y`/`z`
    /// and, unless `bonds` is empty, a bonds block with uint `atomi`/`atomj`.
    fn frame_from_parts(coords: &[[F; 3]], bonds: &[(Idx, Idx)]) -> Frame {
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
            let mut block = Block::new();
            let ai: Vec<Idx> = bonds.iter().map(|&(i, _)| i).collect();
            let aj: Vec<Idx> = bonds.iter().map(|&(_, j)| j).collect();
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

    fn chain_bonds(n: usize) -> Vec<(Idx, Idx)> {
        (0..n.saturating_sub(1) as Idx)
            .map(|i| (i, i + 1))
            .collect()
    }

    fn along_x(n: usize) -> Vec<[F; 3]> {
        (0..n).map(|i| [i as F, 0.0, 0.0]).collect()
    }

    fn bonded_target(coords: &[[F; 3]], count: usize) -> Target {
        Target::new(frame_from_parts(coords, &chain_bonds(coords.len())), count)
    }

    fn open_box() -> SimBox {
        SimBox::ortho(
            array![100.0, 100.0, 100.0],
            array![0.0, 0.0, 0.0],
            [false; 3],
        )
        .expect("orthorhombic box")
    }

    fn assert_angstrom(got: F, want: F, what: &str) {
        assert_eq!(
            got.to_bits(),
            want.to_bits(),
            "{what}: got {got} Å, want {want} Å"
        );
    }

    /// Depth-3 linear hexamer at (i, 0, 0) Å. Exempted 1-2/1-3/1-4 min is the
    /// 1.0 Å bond; scored min is the 1-5 contact at 4.0 Å.
    #[test]
    fn intra_residual_linear_hexamer_depth3() {
        let positions = along_x(6);
        let target = bonded_target(&positions, 1)
            .with_special_bonds(BondDistanceWeights::from_exclusion_depth(3));
        let intra = IntraResidual::from_targets(&[target], &positions, &open_box());
        assert_angstrom(intra.exempted, 1.0, "exempted");
        assert_angstrom(intra.scored, 4.0, "scored");
    }

    /// Folded hexamer, Å. Atoms 0-1-2 are a 3-4-5 right triangle; pair 0–5
    /// (1-6) is 3.0 Å. Depth-3 exempted pairs are all ≥ 3.0 Å.
    const FOLDED_345: [[F; 3]; 6] = [
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [3.0, 4.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 4.0, 3.0],
        [0.0, 0.0, 3.0],
    ];

    #[test]
    fn intra_residual_folded_345_scored_16() {
        let positions = FOLDED_345;
        let target = bonded_target(&positions, 1)
            .with_special_bonds(BondDistanceWeights::from_exclusion_depth(3));
        let intra = IntraResidual::from_targets(&[target], &positions, &open_box());
        assert_angstrom(intra.scored, 3.0, "scored 1-6");
        assert_angstrom(intra.exempted, 3.0, "exempted min");
    }

    /// Same 5-bead chain, two tables: depth 3 scores 1-5 (4.0 Å); depth 1
    /// scores 1-3 (2.0 Å). A constructor that nails depth 3 fails the depth-1
    /// half.
    #[test]
    fn intra_residual_five_bead_depth1_vs_3() {
        let positions = along_x(5);
        let cell = open_box();
        let d3 = IntraResidual::from_targets(
            &[bonded_target(&positions, 1)
                .with_special_bonds(BondDistanceWeights::from_exclusion_depth(3))],
            &positions,
            &cell,
        );
        let d1 = IntraResidual::from_targets(
            &[bonded_target(&positions, 1)
                .with_special_bonds(BondDistanceWeights::from_exclusion_depth(1))],
            &positions,
            &cell,
        );
        assert_angstrom(d3.scored, 4.0, "depth-3 scored");
        assert_angstrom(d3.exempted, 1.0, "depth-3 exempted");
        assert_angstrom(d1.scored, 2.0, "depth-1 scored");
        assert_angstrom(d1.exempted, 1.0, "depth-1 exempted");
    }

    /// Bonded dimer: the only pair is 1-2, so scored is empty (+∞) and
    /// exempted is the 1.0 Å bond.
    #[test]
    fn intra_residual_bonded_dimer_scored_infinite() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let target = bonded_target(&positions, 1)
            .with_special_bonds(BondDistanceWeights::from_exclusion_depth(3));
        let intra = IntraResidual::from_targets(&[target], &positions, &open_box());
        assert!(intra.scored.is_infinite(), "empty scored class is +∞");
        assert_angstrom(intra.exempted, 1.0, "exempted");
    }

    /// Coincident unbonded atoms are a real scored pair at 0.0 Å, not an empty
    /// class.
    #[test]
    fn intra_residual_coincident_unbonded_scored_zero() {
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let target = Target::from_coords(&positions, &[1.0, 1.0], 1);
        let intra = IntraResidual::from_targets(&[target], &positions, &open_box());
        assert_angstrom(intra.scored, 0.0, "scored coincident");
    }

    /// A second copy is a different molecule. Translating copy 1 by +0.25 Å
    /// must not pull intra below the single-copy 1.0 / 4.0 Å goldens.
    #[test]
    fn intra_residual_second_copy_ignored() {
        let copy0 = along_x(6);
        let target = bonded_target(&copy0, 2)
            .with_special_bonds(BondDistanceWeights::from_exclusion_depth(3));
        let mut positions = copy0.clone();
        positions.extend(copy0.iter().map(|&[x, y, z]| [x + 0.25, y, z]));
        let intra = IntraResidual::from_targets(&[target], &positions, &open_box());
        assert_angstrom(intra.exempted, 1.0, "exempted");
        assert_angstrom(intra.scored, 4.0, "scored");
    }

    /// Linear hexamer at x = 0,1,2,3,4,9.5 in a 10 Å PBC-x box: the 1-6
    /// contact wraps to 0.5 Å.
    #[test]
    fn intra_residual_pbc_scored_half() {
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [9.5, 0.0, 0.0],
        ];
        let target = bonded_target(&positions, 1)
            .with_special_bonds(BondDistanceWeights::from_exclusion_depth(3));
        let cell = SimBox::ortho(
            array![10.0, 20.0, 20.0],
            array![0.0, 0.0, 0.0],
            [true, false, false],
        )
        .expect("orthorhombic box");
        let intra = IntraResidual::from_targets(&[target], &positions, &cell);
        assert_angstrom(intra.scored, 0.5, "PBC scored");
    }

    /// `from_coords` has no template: identity exemption only (i==j skipped).
    /// The 1.0 Å pair is scored; exempted is empty (+∞).
    #[test]
    fn intra_residual_bondless_dimer_scores_the_pair() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let target = Target::from_coords(&positions, &[1.0, 1.0], 1);
        let intra = IntraResidual::from_targets(&[target], &positions, &open_box());
        assert_angstrom(intra.scored, 1.0, "scored");
        assert!(
            intra.exempted.is_infinite(),
            "identity-only has no i≠j exempted pair"
        );
    }

    /// Bond (0, 9) on a 3-atom frame is `Topology::from_frame` Validation.
    /// The target is omitted (both +∞), not scored as the 1-2 length 1.5 Å.
    #[test]
    fn intra_residual_validation_omit_not_identity() {
        let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let target = Target::new(frame_from_parts(&positions, &[(0, 9)]), 1);
        let intra = IntraResidual::from_targets(&[target], &positions, &open_box());
        assert!(
            intra.scored.is_infinite(),
            "Validation omit is +∞, not the 1-2 length"
        );
        assert!(intra.exempted.is_infinite(), "Validation omit is +∞");
        assert_ne!(
            intra.scored.to_bits(),
            (1.5 as F).to_bits(),
            "must not report the 1-2 length 1.5 Å as scored"
        );
    }
}
