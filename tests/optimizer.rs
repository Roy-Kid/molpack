//! In-loop optimizers: TorsionMcOptimizer + OptimizeSelect via with_optimizer.

#![cfg(feature = "ff")]

use molpack::{
    F, InsideBoxRestraint, Molpack, OptimizeSelect, Target, TorsionMcOptimizer,
    validate_from_targets,
};
use molrs::store::block::Block;
use molrs::store::frame::Frame;
use molrs::system::atomistic::Atomistic;
use molrs::system::molgraph::Atom;
use ndarray::Array1;

fn chain(n: usize) -> (Atomistic, Vec<[F; 3]>, Vec<F>) {
    let mut g = Atomistic::new();
    let mut ids = Vec::new();
    for _ in 0..n {
        ids.push(g.add_atom(Atom::new()));
    }
    for i in 0..n - 1 {
        g.add_bond(ids[i], ids[i + 1]).expect("bond");
    }
    // Zigzag, not collinear — see `chain_frame` for why torsion optimizers
    // are a silent no-op on a perfectly straight chain.
    let bond_len: F = 1.54;
    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());
    let mut coords = Vec::with_capacity(n);
    for i in 0..n {
        coords.push([i as F * dx, 0.0, if i % 2 == 0 { 0.0 } else { dz }]);
    }
    let radii = vec![1.0; n];
    (g, coords, radii)
}

/// Zigzag bead chain as a `molrs::Frame` (atoms + bonds) plus its graph, so
/// targets are built the normal way (`Target::new`) and keep their topology.
///
/// The zigzag is load-bearing, not cosmetic: in a *collinear* chain every
/// rotatable bond's axis runs through all of its downstream atoms, so torsion
/// rotation is the identity and any torsion optimizer is silently a no-op.
/// Tetrahedral bond angles give the torsions something to move.
fn chain_frame(n: usize, bond_len: F) -> (Frame, Atomistic) {
    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());
    let mut xs = Vec::with_capacity(n);
    let mut zs = Vec::with_capacity(n);
    for i in 0..n {
        xs.push(i as F * dx);
        zs.push(if i % 2 == 0 { 0.0 } else { dz });
    }

    let mut atoms = Block::new();
    atoms
        .insert("x", Array1::from_vec(xs).into_dyn())
        .expect("x column");
    atoms
        .insert("y", Array1::from_vec(vec![0.0 as F; n]).into_dyn())
        .expect("y column");
    atoms
        .insert("z", Array1::from_vec(zs).into_dyn())
        .expect("z column");

    let mut bonds = Block::new();
    let ai: Vec<u32> = (0..n as u32 - 1).collect();
    let aj: Vec<u32> = (1..n as u32).collect();
    bonds
        .insert("atomi", Array1::from_vec(ai).into_dyn())
        .expect("atomi column");
    bonds
        .insert("atomj", Array1::from_vec(aj).into_dyn())
        .expect("atomj column");
    // Deliberately no `bond_type` column — this is the shape a PDB CONECT list
    // or a hand-built coarse-grain frame has. The bonds read back `Unknown`;
    // supplying the single-bond fallback is `TorsionMcOptimizer`'s job.

    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame.insert("bonds", bonds);

    let graph = Atomistic::from_frame(&frame).expect("atomistic from frame");
    (frame, graph)
}

/// Single-bead filler species, used to keep the per-type phases from
/// converging so the all-type phase (where optimizers run) is reached.
fn filler_frame() -> Frame {
    let mut atoms = Block::new();
    for k in ["x", "y", "z"] {
        atoms
            .insert(k, Array1::from_vec(vec![0.0 as F]).into_dyn())
            .expect("column");
    }
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame
}

/// End-to-end distance of copy `c` in a packed chain target laid out
/// copy-by-copy, atom-by-atom from `positions[0]`.
fn end_to_end(positions: &[[F; 3]], copy: usize, natoms: usize) -> F {
    let a = positions[copy * natoms];
    let b = positions[copy * natoms + natoms - 1];
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[test]
fn torsion_mc_optimizer_packs() {
    let (frame, graph) = chain_frame(6, 1.54);
    // Optimizers run only in the all-type phase, and a phase-0 convergence
    // short-circuits before it. A second species plus a box tight enough to
    // stay unsolved through the per-type phases is what makes this test
    // actually reach the optimizer rather than silently skip it.
    let cell = InsideBoxRestraint::new([0.0; 3], [16.0, 16.0, 16.0], [false; 3]);
    let chains = Target::new(frame, 4)
        .with_name("chain")
        .with_restraint(cell);
    let filler = Target::new(filler_frame(), 40)
        .with_name("filler")
        .with_restraint(cell);

    let opt = TorsionMcOptimizer::new(&graph)
        .with_steps(5)
        .with_self_avoidance(1.0);

    let targets = [chains, filler];
    let result = Molpack::new()
        .with_tolerance(2.0)
        .with_seed(3)
        .with_optimizer(OptimizeSelect::per_copy(["chain"]), opt)
        .pack_with_report(&targets, 30)
        .expect("pack");

    let report = validate_from_targets(&targets, &result.positions(), 2.0, 1e-2);
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

    let cell = InsideBoxRestraint::new([0.0; 3], [18.0, 18.0, 18.0], [true, true, true]);
    let target = Target::new(frame.clone(), 8)
        .with_name("water")
        .with_restraint(cell);
    // Second species: without it the per-type phases converge and the
    // all-type phase — the only one that runs optimizers — is skipped.
    let filler = Target::new(filler_frame(), 90)
        .with_name("filler")
        .with_restraint(cell);

    let soft_opt = SoftSpec::from_frame(&frame)
        .with_repulsion(8.0)
        .into_optimizer(0.05, 40, 0.2, 8);

    let result = Molpack::new()
        .with_tolerance(2.0)
        .with_seed(5)
        .with_periodic_box([0.0; 3], [18.0; 3], [true; 3])
        .with_optimizer(
            OptimizeSelect::per_copy(["water"]).with_environment(5.0),
            soft_opt,
        )
        .pack_with_report(&[target, filler], 40)
        .expect("pack");

    assert_eq!(result.natoms(), 16 + 90);
}

// ── per-copy conformations (regression: coor index space) ──────────────────

/// `PerCopy` over a target with `count > 1` must address each copy's own
/// reference conformer. Regression: `run_optimizer_bindings` indexed the
/// per-type `coor` buffer with a per-copy stride and panicked.
#[test]
fn per_copy_handles_multiple_copies() {
    let n = 8;
    let copies = 6;
    let (frame, graph) = chain_frame(n, 1.54);
    let cell = InsideBoxRestraint::new([0.0; 3], [16.0, 16.0, 16.0], [false; 3]);

    let chains = Target::new(frame, copies)
        .with_name("chain")
        .with_restraint(cell);
    // Second species keeps the per-type phases from converging, so the
    // all-type phase — the only one that runs optimizers — is reached.
    let filler = Target::new(filler_frame(), 40)
        .with_name("filler")
        .with_restraint(cell);

    let opt = TorsionMcOptimizer::new(&graph)
        .with_steps(5)
        .with_self_avoidance(1.0);

    let result = Molpack::new()
        .with_tolerance(2.0)
        .with_seed(7)
        .with_optimizer(OptimizeSelect::per_copy(["chain"]), opt)
        .pack_with_report(&[chains, filler], 20)
        .expect("pack");

    assert_eq!(result.natoms(), copies * n + 40);
}

/// Each copy carries its own conformer, so per-copy torsion MC must be able
/// to drive them apart. Regression: all copies shared one `coor` block, so
/// every copy was forced into an identical conformation.
#[test]
fn per_copy_yields_distinct_conformations() {
    let n = 10;
    let copies = 5;
    let (frame, graph) = chain_frame(n, 1.54);
    let cell = InsideBoxRestraint::new([0.0; 3], [18.0, 18.0, 18.0], [false; 3]);

    let chains = Target::new(frame, copies)
        .with_name("chain")
        .with_restraint(cell);
    let filler = Target::new(filler_frame(), 60)
        .with_name("filler")
        .with_restraint(cell);

    let opt = TorsionMcOptimizer::new(&graph)
        .with_steps(20)
        .with_self_avoidance(1.0)
        .with_temperature(1.0);

    let result = Molpack::new()
        .with_tolerance(2.0)
        .with_seed(11)
        .with_optimizer(OptimizeSelect::per_copy(["chain"]), opt)
        .pack_with_report(&[chains, filler], 30)
        .expect("pack");

    let pos = result.positions();
    let e2e: Vec<F> = (0..copies).map(|c| end_to_end(&pos, c, n)).collect();

    let spread =
        e2e.iter().cloned().fold(F::MIN, F::max) - e2e.iter().cloned().fold(F::MAX, F::min);
    assert!(
        spread > 1e-6,
        "per-copy optimization must let conformations diverge, got identical \
         end-to-end distances {e2e:?}"
    );
}

// ── consumer-side fallback for unclassed bonds ─────────────────────────────

/// A frame carrying connectivity but no `bond_type` column — a PDB `CONECT`
/// list, a hand-built coarse-grain chain — reads back with every bond
/// `BondType::Unknown`, and molrs keeps it that way on purpose. Raw perception
/// therefore finds nothing; `TorsionMcOptimizer` supplies the single-bond
/// fallback so it finds the same bonds a graph-built chain has.
#[test]
fn optimizer_falls_back_to_single_for_unclassed_bonds() {
    use molrs::perceive::rotatable::detect_rotatable_bonds_with_downstream;

    let (_frame, from_frame) = chain_frame(10, 1.54);
    let (from_graph, _, _) = chain(10);

    // The reader stays faithful: no stated class, so nothing is rotatable.
    assert_eq!(
        detect_rotatable_bonds_with_downstream(&from_frame).len(),
        0,
        "an unstated bond class must not be inferred by the reader"
    );

    // The consumer applies the policy and recovers the real chain topology.
    assert_eq!(
        TorsionMcOptimizer::new(&from_frame).rotatable_bond_count(),
        TorsionMcOptimizer::new(&from_graph).rotatable_bond_count(),
        "the optimizer's fallback must recover the graph-built chain's bonds"
    );
    assert!(TorsionMcOptimizer::new(&from_frame).rotatable_bond_count() > 0);
}

/// The fallback fills in only what was unstated. A bond the input explicitly
/// classes as double must stay non-rotatable.
#[test]
fn optimizer_fallback_does_not_override_a_stated_class() {
    let (frame, _) = chain_frame(6, 1.54);
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
        TorsionMcOptimizer::new(&one_double).rotatable_bond_count() + 1,
        TorsionMcOptimizer::new(&all_unclassed).rotatable_bond_count(),
        "the stated double bond must stay out of the rotatable set"
    );
}

// ── non-harm gate (regression: stale geometry cache) ───────────────────────

/// An optimizer that always makes things worse: it inflates the molecule 5×.
/// Every proposal it makes must be rejected.
#[derive(Debug)]
struct Exploder;

impl molrs::optimize::Optimizer for Exploder {
    fn run(
        &mut self,
        frame: &mut molrs::store::frame::Frame,
    ) -> Result<molpack::OptReport, String> {
        let flat = molrs::ff::potential::extract_coords(frame)?;
        let blown: Vec<F> = flat.iter().map(|c| c * 5.0).collect();
        molrs::ff::potential::write_coords(frame, &blown)?;
        Ok(molpack::OptReport {
            converged: true,
            n_steps: 1,
            final_energy: 0.0,
            final_fmax: 0.0,
        })
    }
}

/// The packer accepts an optimizer's conformer only if the packing objective
/// did not get worse. Regression: both sides of that comparison were evaluated
/// at the same `x`, and the geometry cache is keyed on `x` and not on `coor` —
/// so the second evaluation returned the first one's value, the gate never
/// fired, and an arbitrarily harmful conformer was written through.
#[test]
fn non_harm_gate_rejects_a_harmful_conformer() {
    let n = 8;
    let bond_len: F = 1.54;
    let (frame, _graph) = chain_frame(n, bond_len);
    let cell = InsideBoxRestraint::new([0.0; 3], [16.0, 16.0, 16.0], [false; 3]);

    let chains = Target::new(frame, 4)
        .with_name("chain")
        .with_restraint(cell);
    let filler = Target::new(filler_frame(), 40)
        .with_name("filler")
        .with_restraint(cell);

    let result = Molpack::new()
        .with_tolerance(2.0)
        .with_seed(17)
        .with_optimizer(OptimizeSelect::per_copy(["chain"]), Exploder)
        .pack_with_report(&[chains, filler], 20)
        .expect("pack");

    // Bond lengths survive: a 5× inflation would show up immediately.
    let pos = result.positions();
    let mut longest = 0.0 as F;
    for c in 0..4 {
        for i in 0..n - 1 {
            let (a, b) = (pos[c * n + i], pos[c * n + i + 1]);
            let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
            longest = longest.max(d);
        }
    }

    assert!(
        longest < bond_len * 1.5,
        "harmful conformers must be reverted; longest bond {longest:.3} \
         (input {bond_len}) means the gate let a 5x inflation through"
    );
}
