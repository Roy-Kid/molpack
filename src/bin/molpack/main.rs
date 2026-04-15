//! `molpack` — Packmol-compatible molecular packing CLI.
//!
//! Usage:
//!   molpack [INPUT]            # read from file
//!   molpack < input.inp        # read from stdin (Packmol-compatible)
//!   cat input.inp | molpack    # pipe

mod io;
mod parser;
mod runner;

use std::io::Read;
use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    name = "molpack",
    version,
    about = "Packmol-compatible molecular packing (molpack)",
    long_about = "\
Pack molecules into a simulation box using a Packmol-compatible .inp script.\n\
Supports all Packmol restraint types plus additional input formats via molrs-io\n\
(SDF/MOL, LAMMPS dump, LAMMPS data) beyond Packmol's PDB/XYZ.\n\
\n\
Examples:\n\
  molpack mixture.inp\n\
  molpack < mixture.inp\n\
  cat mixture.inp | molpack"
)]
struct Args {
    /// Path to the .inp input file. Reads from stdin when omitted.
    input: Option<PathBuf>,
}

fn main() {
    let args = Args::parse();

    // Determine base directory for resolving relative paths in the .inp file.
    // When a file argument is given, paths are relative to its parent directory.
    // When reading from stdin, paths are relative to the current working directory.
    let (src, base_dir) = match args.input {
        Some(ref path) => {
            let src = std::fs::read_to_string(path).unwrap_or_else(|e| {
                eprintln!("Error: cannot read `{}`: {e}", path.display());
                std::process::exit(1);
            });
            let base = path
                .canonicalize()
                .unwrap_or_else(|_| path.to_path_buf())
                .parent()
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| PathBuf::from("."));
            (src, base)
        }
        None => {
            let mut buf = String::new();
            std::io::stdin()
                .read_to_string(&mut buf)
                .unwrap_or_else(|e| {
                    eprintln!("Error: cannot read stdin: {e}");
                    std::process::exit(1);
                });
            let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
            (buf, cwd)
        }
    };

    let parsed = parser::parse_input(&src).unwrap_or_else(|e| {
        eprintln!("Parse error: {e}");
        std::process::exit(1);
    });

    runner::run(parsed, &base_dir).unwrap_or_else(|e| {
        eprintln!("Error: {e}");
        std::process::exit(1);
    });
}
