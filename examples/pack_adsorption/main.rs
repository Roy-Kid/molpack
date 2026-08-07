//! Polymer adsorption onto a flat substrate — a heterogeneous system packed in
//! one call, with the chains relaxing their conformation *while* they pack.
//!
//! # The problem
//!
//! A coarse-grain chain carries a sticky site every `STICKY_EVERY` beads, and
//! every sticky site must end up in a thin slab just above a solid substrate.
//! An extended chain cannot satisfy that: its sticky beads are spread along the
//! chain axis, so no rigid-body placement puts them all on one plane. The chain
//! has to change shape — and it has to do so against the same objective that
//! places it, or the two fight each other.
//!
//! That is what an in-loop optimizer is for. `TorsionMcOptimizer` runs inside
//! the packing loop; the packer accepts its result only if the packing
//! objective did not get worse, so conformational search and placement move
//! together instead of taking turns.
//!
//! # What this demonstrates
//!
//! Heterogeneity along three independent axes at once:
//!
//! | axis | in this system |
//! |---|---|
//! | species | fixed substrate lattice / flexible chains / single-bead solvent |
//! | restraint scope | molecule-wide cell, per-atom-subset slab on sticky beads only |
//! | degrees of freedom | fixed placement / rigid body / rigid body + conformation |
//!
//! The measurable outcome is the classical train–loop–tail decomposition of an
//! adsorbed polymer layer, printed at the end. Geometry is synthesized in
//! process, so no data files and no `io` feature.
//!
//! # Reading the result
//!
//! The run is not expected to reach `converged=true`. Driving every sticky bead
//! into a one-bead-thick slab is a stiff demand on a connected chain, and a
//! real adsorbed layer does not put every site on the surface either — the
//! residual is the entropic cost of the loops and tails, which is the physics
//! this example is about. Judge it by the adsorbed fraction and the
//! train/loop/tail morphology, not by the residual alone.
//!
//! Run with:
//! ```sh
//! cargo run --release --example pack_adsorption --features ff
//! ```
//!
//! | env var | effect |
//! |---|---|
//! | `MOLPACK_ADSORPTION_RIGID=1` | control run: same system, chains kept rigid |
//! | `MOLPACK_ADSORPTION_LOOPS=n` | outer iterations (default 300) |
//! | `MOLPACK_ADSORPTION_XYZ=path` | dump the packed structure |
//! | `MOLPACK_ADSORPTION_PROGRESS=1` | per-iteration progress |
//!
//! The rigid control is the point of comparison: with the chains frozen, the
//! restraints can only be met by placement, and the layer comes out with long
//! tails and a lower adsorbed fraction.

mod analysis;
mod geometry;

use molpack::{
    AbovePlaneRestraint, BelowPlaneRestraint, CenteringMode, F, InsideBoxRestraint, Molpack,
    OptimizeSelect, ProgressHandler, Target, TorsionMcOptimizer,
};

// ── system definition ──────────────────────────────────────────────────────

/// Bead contact distance; atom radii are `TOLERANCE / 2`.
const TOLERANCE: F = 4.0;
/// Substrate lattice: `SUB_N × SUB_N` sites at `SUB_SPACING`, at `z = 0`.
const SUB_N: usize = 16;
const SUB_SPACING: F = 4.0;

const N_CHAINS: usize = 10;
const N_BEADS: usize = 20;
const BOND_LEN: F = 4.0;
/// One sticky site every this many beads.
const STICKY_EVERY: usize = 4;

const N_SOLVENT: usize = 400;

/// Adsorption slab for the sticky beads — the first bead layer above the
/// substrate. The lower edge is not a free choice: substrate beads sit at
/// `z = 0` with radius `TOLERANCE / 2`, so any chain bead centre must clear
/// `TOLERANCE` to avoid overlapping them. A slab below that is infeasible by
/// construction and the restraint can only ever be violated.
const SLAB_LO: F = TOLERANCE;
const SLAB_HI: F = TOLERANCE + BOND_LEN;
/// A bead below this z counts as adsorbed in the train/loop/tail analysis.
const ADSORBED_BELOW: F = SLAB_HI + 1.0;
/// Soft wall keeping everything out of the substrate plane — same floor as the
/// slab, for the same excluded-volume reason.
const WALL_Z: F = TOLERANCE;
const CELL_Z: F = 40.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::try_init();

    let span = SUB_N as F * SUB_SPACING;
    let half = span / 2.0;
    // Periodic in x/y, open in z: a slab geometry with one solid face.
    let cell = InsideBoxRestraint::new(
        [-half, -half, 0.0],
        [half, half, CELL_Z],
        [true, true, false],
    );
    let wall = AbovePlaneRestraint::new([0.0, 0.0, 1.0], WALL_Z);

    let (chain_frame, chain_graph, sticky) = geometry::chain(N_BEADS, BOND_LEN, STICKY_EVERY);

    // 1. Substrate — one fixed copy, kept exactly where it was built.
    let substrate = Target::new(geometry::substrate(SUB_N, SUB_N, SUB_SPACING, 0.0), 1)
        .with_name("substrate")
        .with_centering(CenteringMode::Off)
        .fixed_at([0.0, 0.0, 0.0]);

    // 2. Chains — one target, `N_CHAINS` copies. Each copy carries its own
    //    conformer, so the optimizer can drive them apart; the slab restraint
    //    applies to the sticky beads *only*, leaving the rest free to loop.
    let chains = Target::new(chain_frame, N_CHAINS)
        .with_name("chain")
        .with_restraint(cell)
        .with_restraint(wall)
        .with_atom_restraint(&sticky, AbovePlaneRestraint::new([0.0, 0.0, 1.0], SLAB_LO))
        .with_atom_restraint(&sticky, BelowPlaneRestraint::new([0.0, 0.0, 1.0], SLAB_HI));

    // 3. Solvent — rigid single beads filling the rest of the cell.
    let solvent = Target::new(geometry::solvent_bead(), N_SOLVENT)
        .with_name("solvent")
        .with_restraint(cell)
        .with_restraint(wall);

    // Torsion MC on each chain copy independently, with the substrate and
    // neighbours frozen but *visible* inside the cutoff, so a chain folds
    // against its real local environment rather than an abstract plane.
    let torsion = TorsionMcOptimizer::new(&chain_graph)
        .with_steps(25)
        .with_temperature(0.6)
        .with_self_avoidance(TOLERANCE / 2.0);
    assert!(
        torsion.rotatable_bond_count() > 0,
        "no rotatable bonds perceived — the optimizer would be a silent no-op"
    );
    println!(
        "chain: {N_BEADS} beads, {} sticky sites, {} rotatable bonds",
        sticky.len(),
        torsion.rotatable_bond_count()
    );

    let mut packer = Molpack::new()
        .with_tolerance(TOLERANCE)
        .with_seed(20_260_807);
    // Control switch: packing the same system with rigid chains is what shows
    // the in-loop optimizer is doing the work, rather than the restraints being
    // satisfiable by placement alone.
    if std::env::var_os("MOLPACK_ADSORPTION_RIGID").is_none() {
        packer = packer.with_optimizer(
            OptimizeSelect::per_copy(["chain"]).with_environment(8.0),
            torsion,
        );
    }
    if std::env::var_os("MOLPACK_ADSORPTION_PROGRESS").is_some() {
        packer = packer.with_handler(ProgressHandler::new());
    }

    let targets = [substrate, chains, solvent];
    let max_loops = std::env::var("MOLPACK_ADSORPTION_LOOPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(300);
    let result = packer.pack_with_report(&targets, max_loops)?;

    report(&result, &sticky)?;
    Ok(())
}

// ── reporting ──────────────────────────────────────────────────────────────

fn report(
    result: &molpack::PackResult,
    sticky: &[usize],
) -> Result<(), Box<dyn std::error::Error>> {
    let positions = result.positions();
    let n_sub = SUB_N * SUB_N;
    let n_chain_atoms = N_CHAINS * N_BEADS;
    // Targets come back in declared order: substrate, chains, solvent.
    let chain_positions = &positions[n_sub..n_sub + n_chain_atoms];

    let stats = analysis::layer_stats(chain_positions, N_CHAINS, N_BEADS, sticky, ADSORBED_BELOW);

    println!(
        "\npack: {} atoms, converged={}, fdist={:.3e}, frest={:.3e}",
        result.natoms(),
        result.converged,
        result.fdist,
        result.frest
    );

    println!("\n── adsorbed layer ─────────────────────────────");
    println!(
        "  sticky beads in the slab : {}/{} ({:.0}%)",
        stats.sticky_adsorbed,
        stats.sticky_total,
        100.0 * stats.sticky_fraction()
    );
    println!(
        "  trains : {:>3}   mean length {:.1}",
        stats.trains.len(),
        analysis::LayerStats::mean(&stats.trains)
    );
    println!(
        "  loops  : {:>3}   mean length {:.1}",
        stats.loops.len(),
        analysis::LayerStats::mean(&stats.loops)
    );
    println!(
        "  tails  : {:>3}   mean length {:.1}",
        stats.tails.len(),
        analysis::LayerStats::mean(&stats.tails)
    );

    println!("\n── chain-bead z density ───────────────────────");
    let profile = analysis::z_profile(chain_positions, 0.0, CELL_Z, 16);
    print!("{}", analysis::render_profile(&profile, 40));

    if std::env::var_os("MOLPACK_ADSORPTION_RIGID").is_some() {
        println!("\n(rigid control run — chains were not relaxed in-loop)");
    } else {
        println!("\n(re-run with MOLPACK_ADSORPTION_RIGID=1 for the rigid-chain control)");
    }

    if let Some(path) = std::env::var_os("MOLPACK_ADSORPTION_XYZ") {
        write_xyz(std::path::Path::new(&path), result)?;
        println!("\nwrote {}", std::path::Path::new(&path).display());
    }
    Ok(())
}

/// Minimal XYZ writer — keeps the example free of the `io` feature.
fn write_xyz(
    path: &std::path::Path,
    result: &molpack::PackResult,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::io::Write;

    let positions = result.positions();
    let atoms = result
        .frame
        .get("atoms")
        .ok_or("result has no atoms block")?;
    let elements = atoms.get_string("element");

    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(out, "{}", positions.len())?;
    writeln!(out, "molpack pack_adsorption")?;
    for (i, p) in positions.iter().enumerate() {
        let sym = elements.map(|c| c[[i]].as_str()).unwrap_or("X");
        writeln!(out, "{sym} {:.4} {:.4} {:.4}", p[0], p[1], p[2])?;
    }
    Ok(())
}
