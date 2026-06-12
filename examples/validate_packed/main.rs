//! Score an externally-produced packed structure against the targets implied by
//! a Packmol `.inp` script, using molpack's own validator. Prints one JSON line
//! of packing-quality metrics. Lets the molpack-jcc benchmark score *both*
//! Packmol's and molpack's outputs through the identical validator at a common
//! tolerance and precision.
//!
//! Usage:
//!   validate_packed <script.inp> <packed_file> [tolerance=2.0] [precision=1e-2]

use std::path::Path;

use molpack::frame_to_coords;
use molpack::script::{self, read_frame};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 2 {
        eprintln!("usage: validate_packed <script.inp> <packed_file> [tolerance] [precision]");
        std::process::exit(2);
    }
    let inp = Path::new(&args[0]);
    let packed = Path::new(&args[1]);
    let tolerance: f64 = args
        .get(2)
        .map(|s| s.parse().expect("tolerance"))
        .unwrap_or(2.0);
    let precision: f64 = args
        .get(3)
        .map(|s| s.parse().expect("precision"))
        .unwrap_or(1e-2);

    // Parse the official .inp and lower it to targets — the same path the CLI
    // uses, so the molecule counts and declared order match what both packers
    // wrote.
    let src = std::fs::read_to_string(inp).expect("read .inp");
    let base = inp
        .canonicalize()
        .ok()
        .and_then(|p| p.parent().map(|q| q.to_path_buf()))
        .unwrap_or_else(|| Path::new(".").to_path_buf());
    let ast = script::parse(&src).expect("parse .inp");
    let built = ast.build(&base).expect("build targets");

    // Load the packed coordinates (PDB/XYZ inferred from extension).
    let frame = read_frame(packed, None).expect("read packed file");
    let (coords, _radii) = frame_to_coords(&frame);

    let r = molpack::validate_from_targets(&built.targets, &coords, tolerance, precision);
    let m = &r.metrics;
    println!(
        "{{\"file\":\"{}\",\"expected_atoms\":{},\"actual_atoms\":{},\"atom_count_ok\":{},\
\"expected_molecules\":{},\"molecule_count_ok\":{},\
\"max_distance_violation\":{:.6e},\"max_constraint_penalty\":{:.6e},\
\"violating_pairs\":{},\"violating_atoms\":{},\"valid\":{}}}",
        packed.display(),
        r.expected_atoms,
        r.actual_atoms,
        r.atom_count_ok,
        r.expected_molecules,
        r.molecule_count_ok,
        m.max_distance_violation,
        m.max_constraint_penalty,
        m.violating_pairs,
        m.violating_atoms,
        r.is_valid(),
    );
}
