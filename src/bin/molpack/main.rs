//! `molpack` — molecular packing CLI.
//!
//! Usage:
//!   molpack [INPUT]            # read from file
//!   molpack < input.inp        # read from stdin
//!   cat input.inp | molpack    # pipe

use std::io::Read;
use std::path::PathBuf;

use clap::Parser;

use molpack::ProgressHandler;
use molpack::script::{self, BuildResult, ScriptError};

#[derive(Parser, Debug)]
#[command(
    name = "molpack",
    version,
    about = "Molecular packing CLI (molpack)",
    long_about = "\
Pack molecules into a simulation box from a molpack `.inp` script.\n\
Accepts the full `.inp` keyword set plus additional input formats via\n\
molrs-io (SDF/MOL, LAMMPS dump, LAMMPS data) alongside the standard\n\
PDB and XYZ.\n\
\n\
Examples:\n\
  molpack mixture.inp\n\
  molpack < mixture.inp\n\
  cat mixture.inp | molpack"
)]
struct Args {
    /// Path to the .inp script. Reads from stdin when omitted.
    input: Option<PathBuf>,
}

fn main() {
    let args = Args::parse();

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

    if let Err(e) = run(&src, &base_dir) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

fn run(src: &str, base_dir: &std::path::Path) -> Result<(), ScriptError> {
    let script_ast = script::parse(src)?;
    let BuildResult {
        mut packer,
        targets,
        output,
        nloop,
    } = script_ast.build(base_dir)?;

    // Attach the progress handler only for the CLI; the library stays headless.
    packer = packer.with_handler(ProgressHandler::new());

    let result = packer.pack(&targets, nloop)?;

    println!(
        "Packing complete — converged: {}, fdist: {:.6}, frest: {:.6}, natoms: {}",
        result.converged,
        result.fdist,
        result.frest,
        result.natoms()
    );

    if !result.converged {
        eprintln!(
            "Warning: packing did not fully converge (fdist={:.6}, frest={:.6}). \
             Consider increasing `nloop` or adjusting restraints.",
            result.fdist, result.frest
        );
    }

    script::write_frame(&output, &result.frame)?;
    println!("Output written to: {}", output.display());

    Ok(())
}
