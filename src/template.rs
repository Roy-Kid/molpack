//! What molpack reads off a template: Cartesian coordinate rows, and the
//! rotatable bonds of its bond graph.
//!
//! Crate-root leaf over molrs only. Coordinates are in Å (the crate's length
//! unit, as carried by the frame — nothing is converted).

use molrs::core::Atomistic;
use molrs::op::{F, Fnx3};
use molrs::perceive::{RotatableBond, UnknownBondPolicy, perceive_rotatable_bonds_with_downstream};

/// Coordinate rows of an `N × 3` array (as `molrs::core::Frame::coords` returns it), in
/// the array's row order — the `[x, y, z]` triples molpack's placement code
/// works in.
pub(crate) fn coord_rows(xyz: &Fnx3) -> Vec<[F; 3]> {
    xyz.rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect()
}

/// The rotatable bonds of a template graph, each with its downstream atoms.
///
/// Formats that carry connectivity without orders — PDB `CONECT`, GROMACS
/// `.top`, XYZ `Connct`, a hand-built coarse-grain frame — read back
/// `BondOrder::Unknown`, because molrs reports what the file said rather than
/// guessing. molpack's policy, held here once for every solver: an unclassed
/// bond is a rotatable single bond ([`UnknownBondPolicy::AsSingle`]). A class
/// the input did state is kept.
pub(crate) fn rotatable_bonds(graph: &Atomistic) -> Vec<RotatableBond> {
    perceive_rotatable_bonds_with_downstream(graph, UnknownBondPolicy::AsSingle)
}

#[cfg(test)]
mod rotatable_tests {
    //! The single-bond fallback: a bond whose class the input left unstated is
    //! treated as rotatable, and one the input did state is not overridden.
    //! Raw molrs perception supplies neither — this is molpack's own policy,
    //! so it is pinned here.

    use super::rotatable_bonds;
    use crate::testutil::{chain_frame, chain_graph};
    use molrs::core::Atomistic;
    use ndarray::Array1;

    /// A frame carrying connectivity but no `bond_type` column — a PDB `CONECT`
    /// list, a hand-built coarse-grain chain — reads back with every bond
    /// `BondOrder::Unknown`, and molrs keeps it that way on purpose. Raw perception
    /// therefore finds nothing; `rotatable_bonds` supplies the single-bond
    /// fallback so it finds the same bonds a graph-built chain has.
    #[test]
    fn unclassed_bonds_fall_back_to_single() {
        use molrs::perceive::{UnknownBondPolicy, perceive_rotatable_bonds_with_downstream};

        let from_frame = Atomistic::from_frame(&chain_frame(10, 1.54)).expect("graph");
        let from_graph = chain_graph(10);

        // The reader stays faithful: no stated class, so nothing is rotatable.
        assert_eq!(
            perceive_rotatable_bonds_with_downstream(&from_frame, UnknownBondPolicy::NotRotatable)
                .len(),
            0,
            "an unstated bond class must not be inferred by the reader"
        );

        // The consumer applies the policy and recovers the real chain topology.
        assert_eq!(
            rotatable_bonds(&from_frame).len(),
            rotatable_bonds(&from_graph).len(),
            "the fallback must recover the graph-built chain's bonds"
        );
        assert!(!rotatable_bonds(&from_frame).is_empty());
    }

    /// The fallback fills in only what was unstated. A bond the input explicitly
    /// classes as double must stay non-rotatable.
    #[test]
    fn fallback_does_not_override_a_stated_class() {
        let frame = chain_frame(6, 1.54);
        let mut with_double = frame.clone();
        {
            let bonds = with_double
                .get_mut("bonds")
                .expect("chain frame has a bonds block");
            let n = bonds.nrows().expect("bond rows");
            let mut classes = vec![1u32; n];
            classes[2] = 2; // one stated double bond
            bonds
                .insert("bond_type", Array1::from_vec(classes).into_dyn())
                .expect("bond_type column");
        }

        let all_unclassed = Atomistic::from_frame(&frame).expect("graph");
        let one_double = Atomistic::from_frame(&with_double).expect("graph");

        assert_eq!(
            rotatable_bonds(&one_double).len() + 1,
            rotatable_bonds(&all_unclassed).len(),
            "the stated double bond must stay out of the rotatable set"
        );
    }
}
