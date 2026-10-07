//! Polymer translocation initial state — a configuration molecular dynamics
//! cannot reach, only build.
//!
//! # Why this one is different
//!
//! A packing problem MD could also solve, given enough steps, only makes setup
//! *faster*. This one is a **construction** problem, not a sampling one: the
//! target state has chains threaded through pores in a rigid membrane — part of
//! each chain on the cis side, part inside a channel, part on the trans side —
//! and dynamics started from unthreaded chains does not find it.
//!
//! That claim is measured, not assumed. The same system was built as a LAMMPS
//! CG model (WCA + harmonic bonds, frozen membrane, Langevin at T*=1) and run
//! for 2e6 steps from each starting point:
//!
//! | pore | from this packed state | from unthreaded chains |
//! |---|---|---|
//! | wide (rim radius 9 A) | leaks out over ~1e6 steps | threads spontaneously |
//! | narrow (rim radius 6.5 A) | leaks out over ~7e5 steps | **0 threading events** |
//! | narrow + attractive channel | persists, 3 -> 1-2 chains | **0-1 transient events** |
//!
//! So the honest statement is the second column, not a blanket "MD cannot do
//! this": with a channel narrow enough to admit one bead at a time, dynamics
//! essentially never *finds* the threaded state, while packing produces it in
//! seconds. Whether the state then *persists* is set by the physics applied
//! afterwards — a bare repulsive channel drains, a functionalized one holds —
//! which is a property of the model, not of how the state was built.
//!
//! This is what packing is for: reaching a configuration that is legitimate but
//! kinetically inaccessible, so a study of what happens *next* can start from
//! it. Translocation work builds such states (or pulls them through with a
//! steered run) for exactly this reason.
//!
//! # Heterogeneity
//!
//! | axis | in this system |
//! |---|---|
//! | species | rigid membrane / threaded chains / free cis chains / solvent |
//! | restraint scope | whole-molecule cell; **three different geometries on three different bead subsets of the same molecule** |
//! | restraint geometry | half-space (cis), cylinder (pore), half-space (trans) |
//! | degrees of freedom | fixed placement / rigid body / rigid body + conformation |
//!
//! The free cis chains are the same species as the threaded ones and differ
//! only in what they are asked to satisfy — the clearest statement that a
//! "species" here is a restraint set, not a molecule type.
//!
//! Run with:
//! ```sh
//! cargo run --release --example pack_translocation --features io
//! ```
//! `MOLPACK_TRANSLOCATION_XYZ=path` dumps the structure;
//! `MOLPACK_TRANSLOCATION_LOOPS=n` sets the outer iteration count.

mod geometry;

use molpack::{
    CenteringMode, GencanPack, OptimizeSelect, PackEngine, RegionRestraint, Target,
    TorsionMcOptimizer,
};
use molrs::op::F;
use std::sync::Arc;

use molrs::core::{Cuboid, Cylinder, HalfSpace, NotRegion};
use ndarray::array;

// ── molrs regions lifted to "stay inside" (the one geometric restraint) ─────

fn inside_box(min: [F; 3], max: [F; 3]) -> RegionRestraint {
    RegionRestraint(Arc::new(Cuboid::new(
        array![min[0], min[1], min[2]],
        array![max[0] - min[0], max[1] - min[1], max[2] - min[2]],
    )))
}

/// `n · x >= d`: the complement of the half-space behind the plane.
fn above_plane(normal: [F; 3], distance: F) -> RegionRestraint {
    let n = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    let point = [
        distance * normal[0] / n,
        distance * normal[1] / n,
        distance * normal[2] / n,
    ];
    RegionRestraint(Arc::new(NotRegion::new(Arc::new(
        HalfSpace::new(normal, point).expect("plane"),
    ))))
}

/// `n · x <= d`: the half-space behind the plane.
fn below_plane(normal: [F; 3], distance: F) -> RegionRestraint {
    let n = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    let point = [
        distance * normal[0] / n,
        distance * normal[1] / n,
        distance * normal[2] / n,
    ];
    RegionRestraint(Arc::new(HalfSpace::new(normal, point).expect("plane")))
}

fn inside_cylinder(base: [F; 3], axis: [F; 3], radius: F, length: F) -> RegionRestraint {
    RegionRestraint(Arc::new(
        Cylinder::new(base, axis, radius, length).expect("cylinder"),
    ))
}

// ── system ─────────────────────────────────────────────────────────────────

const TOLERANCE: F = 4.0;
const BEAD_R: F = TOLERANCE / 2.0;

/// Membrane: two lattice sheets with a hole punched through both.
const MEM_N: usize = 16;
const MEM_SPACING: F = 4.0;
const MEM_Z: [F; 2] = [0.0, 4.0];
/// Lattice sites inside this radius are removed. A chain bead centre can reach
/// `HOLE_R - TOLERANCE` from the axis before it touches the rim.
const HOLE_R: F = 6.5;

const N_BEADS: usize = 24;
const BOND_LEN: F = 4.0;
/// Beads held on the cis side / inside the channel / on the trans side.
const CIS_LEN: usize = 5;
const PORE_LEN: usize = 3;

/// One pore per threaded chain. A single pore cannot hold several chains: four
/// bead segments need more excluded volume than a channel this narrow has, so
/// over-subscribing one hole makes the restraint unsatisfiable rather than
/// merely hard. Each chain therefore gets its own target with its own cylinder —
/// same molecule template, different restraint set.
const PORES: [[F; 2]; 3] = [[-16.0, -16.0], [16.0, -16.0], [0.0, 16.0]];
const N_FREE: usize = 6;
const N_SOLVENT: usize = 300;

/// Cis chamber is below the membrane, trans above it.
const CIS_Z: F = -6.0;
const TRANS_Z: F = 10.0;
const Z_LO: F = -34.0;
const Z_HI: F = 38.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::try_init();

    let span = MEM_N as F * MEM_SPACING;
    let half = span / 2.0;
    let cell = inside_box([-half, -half, Z_LO], [half, half, Z_HI]);

    let (chain_frame, chain_graph, seg) = geometry::chain(N_BEADS, BOND_LEN, CIS_LEN, PORE_LEN);

    let membrane = Target::new(
        geometry::membrane(MEM_N, MEM_SPACING, &MEM_Z, &PORES, HOLE_R),
        1,
    )
    .with_name("membrane")
    .with_centering(CenteringMode::Off)
    .fixed_at([0.0, 0.0, 0.0]);

    // The threaded chains. Three restraints, three bead subsets, three
    // geometries — all on one molecule, and a different pore per chain.
    let pore_len_axial = (PORE_LEN as F - 1.0) * BOND_LEN;
    let pore_z0 = 0.5 * (MEM_Z[0] + MEM_Z[1]) - 0.5 * pore_len_axial;
    let threaded: Vec<Target> = PORES
        .iter()
        .enumerate()
        .map(|(i, hole)| {
            let pore = inside_cylinder(
                [hole[0], hole[1], pore_z0],
                [0.0, 0.0, 1.0],
                HOLE_R - TOLERANCE,
                pore_len_axial,
            );
            Target::new(chain_frame.clone(), 1)
                .with_name(format!("threaded{i}"))
                .with_restraint(cell.clone())
                .with_atom_restraint(&seg.cis, below_plane([0.0, 0.0, 1.0], CIS_Z))
                .with_atom_restraint(&seg.pore, pore)
                .with_atom_restraint(&seg.trans, above_plane([0.0, 0.0, 1.0], TRANS_Z))
        })
        .collect();

    // Same molecule, different ask: stay in the cis chamber, unthreaded.
    let free = Target::new(chain_frame, N_FREE)
        .with_name("free")
        .with_restraint(cell.clone())
        .with_restraint(below_plane([0.0, 0.0, 1.0], CIS_Z));

    let solvent = Target::new(geometry::solvent_bead(), N_SOLVENT)
        .with_name("solvent")
        .with_restraint(cell);

    let torsion = TorsionMcOptimizer::new(&chain_graph)
        .with_steps(30)
        .with_temperature(0.6)
        .with_self_avoidance(BEAD_R);
    assert!(
        torsion.rotatable_bond_count() > 0,
        "no rotatable bonds perceived — the optimizer would be a silent no-op"
    );

    println!(
        "chain {N_BEADS} beads | cis {:?} pore {:?} trans {:?} | {} rotatable bonds",
        seg.cis,
        seg.pore,
        seg.trans,
        torsion.rotatable_bond_count()
    );

    let max_loops = std::env::var("MOLPACK_TRANSLOCATION_LOOPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(400);

    let mut targets = vec![membrane];
    targets.extend(threaded);
    targets.push(free);
    targets.push(solvent);
    let mut names: Vec<String> = (0..PORES.len()).map(|i| format!("threaded{i}")).collect();
    names.push("free".to_string());
    let result = GencanPack::new()
        .with_tolerance(TOLERANCE)
        .with_seed(20_260_807)
        .with_periodic_box(
            [-half, -half, Z_LO],
            [half, half, Z_HI],
            [true, true, false],
        )
        .with_optimizer(
            OptimizeSelect::per_copy(names).with_environment(8.0),
            torsion,
        )
        .run(&targets, max_loops)?;

    report(&result, &seg)?;
    Ok(())
}

// ── reporting ──────────────────────────────────────────────────────────────

fn report(
    result: &molpack::State,
    seg: &geometry::Segments,
) -> Result<(), Box<dyn std::error::Error>> {
    let pore_len_axial = (PORE_LEN as F - 1.0) * BOND_LEN;
    let pore_z0 = 0.5 * (MEM_Z[0] + MEM_Z[1]) - 0.5 * pore_len_axial;
    let pos = result.positions();
    let n_mem = result.natoms() - (PORES.len() + N_FREE) * N_BEADS - N_SOLVENT;

    println!(
        "\npack: {} atoms ({n_mem} membrane), fdist={:.3e}, frest={:.3e}",
        result.natoms(),
        result.fdist,
        result.frest
    );

    // Report the worst residual rather than a pass/fail at the exact boundary:
    // the restraints are soft penalties, so "inside" is a matter of degree and
    // a boolean would call a 0.5 A overhang a failure.
    println!("\n── threading check (worst residual per chain, A) ───────");
    println!(
        "  {:>5}  {:>9}  {:>9}  {:>9}  {:>9}",
        "chain", "cis", "pore r", "pore z", "trans"
    );
    let mut worst_overall = 0.0 as F;
    for c in 0..PORES.len() {
        let bead = |i: usize| pos[n_mem + c * N_BEADS + i];
        let cis_v = seg
            .cis
            .iter()
            .map(|&i| (bead(i)[2] - CIS_Z).max(0.0))
            .fold(0.0 as F, F::max);
        let trans_v = seg
            .trans
            .iter()
            .map(|&i| (TRANS_Z - bead(i)[2]).max(0.0))
            .fold(0.0 as F, F::max);
        let (mut pore_r, mut pore_z) = (0.0 as F, 0.0 as F);
        for &i in &seg.pore {
            let p = bead(i);
            let r = ((p[0] - PORES[c][0]).powi(2) + (p[1] - PORES[c][1]).powi(2)).sqrt();
            pore_r = pore_r.max((r - (HOLE_R - TOLERANCE)).max(0.0));
            pore_z = pore_z
                .max((pore_z0 - p[2]).max(0.0))
                .max((p[2] - (pore_z0 + pore_len_axial)).max(0.0));
        }
        let worst = cis_v.max(trans_v).max(pore_r).max(pore_z);
        worst_overall = worst_overall.max(worst);
        println!("  {c:>5}  {cis_v:>9.3}  {pore_r:>9.3}  {pore_z:>9.3}  {trans_v:>9.3}");
    }
    println!(
        "\n  worst residual over all {} threaded chains: {worst_overall:.3} A ({:.0}% of a bead radius)",
        PORES.len(),
        100.0 * worst_overall / BEAD_R
    );

    if let Some(path) = std::env::var_os("MOLPACK_TRANSLOCATION_XYZ") {
        molrs::io::write_xyz(&path, &result.frame)?;
        println!("wrote {}", std::path::Path::new(&path).display());
    }
    Ok(())
}
