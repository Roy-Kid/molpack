//! In-loop optimizers: TorsionMcOptimizer + OptimizeSelect via with_optimizer.

#![cfg(feature = "ff")]

use molpack::{
    F, InsideBoxRestraint, Molpack, OptimizeSelect, Target, TorsionMcOptimizer,
    validate_from_targets,
};
use molrs::system::atomistic::Atomistic;
use molrs::system::molgraph::Atom;

fn chain(n: usize) -> (Atomistic, Vec<[F; 3]>, Vec<F>) {
    let mut g = Atomistic::new();
    let mut ids = Vec::new();
    for _ in 0..n {
        ids.push(g.add_atom(Atom::new()));
    }
    for i in 0..n - 1 {
        g.add_bond(ids[i], ids[i + 1]).expect("bond");
    }
    let bond_len: F = 1.54;
    let mut coords = Vec::with_capacity(n);
    for i in 0..n {
        coords.push([i as F * bond_len, 0.0, 0.0]);
    }
    let radii = vec![1.0; n];
    (g, coords, radii)
}

#[test]
fn torsion_mc_optimizer_packs() {
    let (graph, coords, radii) = chain(6);
    let target = Target::from_coords(&coords, &radii, 1)
        .with_name("chain")
        .with_restraint(InsideBoxRestraint::new(
            [0.0; 3],
            [40.0, 40.0, 40.0],
            [false; 3],
        ));

    let opt = TorsionMcOptimizer::new(&graph)
        .with_steps(5)
        .with_self_avoidance(1.0);

    let result = Molpack::new()
        .with_tolerance(2.0)
        .with_seed(3)
        .with_optimizer(OptimizeSelect::per_copy(["chain"]), opt)
        .pack_with_report(std::slice::from_ref(&target), 30)
        .expect("pack");

    let report = validate_from_targets(&[target], &result.positions(), 2.0, 1e-2);
    assert!(report.is_valid(), "{report:?}");
}

#[test]
fn soft_spec_per_copy_with_environment() {
    use molrs::ff::potential::soft::SoftSpec;
    use molrs::store::block::Block;
    use molrs::store::frame::Frame;
    use ndarray::Array1;

    // Two-atom "water"
    let mut atoms = Block::new();
    atoms
        .insert("x", Array1::from_vec(vec![0.0, 2.5]).into_dyn())
        .unwrap();
    atoms
        .insert("y", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
        .unwrap();
    atoms
        .insert("z", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
        .unwrap();
    let mut bonds = Block::new();
    bonds
        .insert("atomi", Array1::from_vec(vec![0u32]).into_dyn())
        .unwrap();
    bonds
        .insert("atomj", Array1::from_vec(vec![1u32]).into_dyn())
        .unwrap();
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame.insert("bonds", bonds);

    let target = Target::new(frame.clone(), 8)
        .with_name("water")
        .with_restraint(InsideBoxRestraint::new(
            [0.0; 3],
            [18.0, 18.0, 18.0],
            [true, true, true],
        ));

    let soft_opt = SoftSpec::from_frame(&frame)
        .with_repulsion(8.0)
        .into_optimizer(0.05, 40, 0.2, 8);

    let result = Molpack::new()
        .with_tolerance(2.0)
        .with_seed(5)
        .with_periodic_box([0.0; 3], [18.0; 3])
        .with_optimizer(
            OptimizeSelect::per_copy(["water"]).with_environment(5.0),
            soft_opt,
        )
        .pack_with_report(std::slice::from_ref(&target), 40)
        .expect("pack");

    assert_eq!(result.natoms(), 16);
}
