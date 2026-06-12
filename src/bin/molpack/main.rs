//! `molpack` — molecular packing CLI.
//!
//! Usage:
//!   molpack [INPUT]            # read from file
//!   molpack < input.inp        # read from stdin
//!   cat input.inp | molpack    # pipe

use std::io::Read;
use std::path::PathBuf;

use clap::Parser;

use molpack::MolpackLogLevel;
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

    /// Run the pair-gradient evaluation on rayon worker threads. Requires a
    /// binary built with the `rayon` feature; errors otherwise.
    #[arg(short, long)]
    parallel: bool,

    /// Number of rayon worker threads (implies --parallel). Defaults to the
    /// rayon global pool size (CPU count or `RAYON_NUM_THREADS`).
    #[arg(short, long)]
    threads: Option<usize>,
}

fn main() {
    let args = Args::parse();

    let parallel = args.parallel || args.threads.is_some();
    if parallel {
        if let Err(e) = configure_parallel(args.threads) {
            eprintln!("Error: {e}");
            std::process::exit(1);
        }
    }

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

    if let Err(e) = run(&src, &base_dir, parallel) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

/// Configure the rayon global thread pool and verify the binary can actually
/// run in parallel. Fail-fast mirrors the Python binding: requesting
/// parallelism from a serial build is an error, not a silent no-op.
fn configure_parallel(threads: Option<usize>) -> Result<(), String> {
    #[cfg(not(feature = "rayon"))]
    let _ = threads;
    if !cfg!(feature = "rayon") {
        return Err(
            "--parallel/--threads requested but this binary was built without the \
             `rayon` feature; rebuild with `cargo build --features cli,rayon`"
                .to_string(),
        );
    }
    #[cfg(feature = "rayon")]
    if let Some(n) = threads {
        if n == 0 {
            return Err("--threads must be >= 1".to_string());
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| format!("failed to configure {n} rayon threads: {e}"))?;
    }
    Ok(())
}

fn run(src: &str, base_dir: &std::path::Path, parallel: bool) -> Result<(), ScriptError> {
    let script_ast = script::parse(src)?;
    let BuildResult {
        mut packer,
        targets,
        output,
        nloop,
    } = script_ast.build(base_dir)?;

    // CLI defaults to screen output; library callers stay headless unless
    // configured on the builder.
    packer = packer.with_log_level(MolpackLogLevel::Progress);
    if parallel {
        packer = packer.with_parallel_eval(true);
    }

    let frame = packer.pack(&targets, nloop)?;
    script::write_frame(&output, &frame)?;
    println!("Output written to: {}", output.display());

    Ok(())
}
