//! Converts [`ParsedInput`] into molpack API calls and writes the result.

use std::path::{Path, PathBuf};

use molpack::{
    AbovePlaneRestraint, Angle, BelowPlaneRestraint, CenteringMode, InsideBoxRestraint,
    InsideSphereRestraint, Molpack, OutsideSphereRestraint, ProgressHandler, Target,
};

use crate::io;
use crate::parser::{ParsedAtomGroup, ParsedInput, ParsedRestraint, ParsedStructure};

/// Resolve `path` against `base` if `path` is relative.
fn resolve(base: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
}

/// Run packing from a parsed input.
///
/// `base_dir` is used to resolve relative paths inside the .inp file.
/// Pass the directory containing the .inp file, or the current working directory
/// when reading from stdin.
pub fn run(input: ParsedInput, base_dir: &Path) -> Result<(), String> {
    if input.structures.is_empty() {
        return Err("no `structure` blocks found in input".into());
    }

    // Build Molpack.
    let mut packer = Molpack::new()
        .with_tolerance(input.tolerance)
        .with_handler(ProgressHandler::new());
    if let Some(seed) = input.seed {
        packer = packer.with_seed(seed);
    }

    // Build Target list.
    let targets: Vec<Target> = input
        .structures
        .iter()
        .map(|s| build_target(s, input.filetype.as_deref(), base_dir))
        .collect::<Result<_, _>>()?;

    // Run packing.
    let result = packer
        .pack(&targets, input.nloop)
        .map_err(|e| format!("packing failed: {e}"))?;

    // Report convergence.
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

    // Write output (resolve relative to base_dir).
    let output_path = resolve(base_dir, &input.output);
    io::write_frame(&output_path, &result.frame)?;
    println!("Output written to: {}", output_path.display());

    Ok(())
}

fn build_target(
    s: &ParsedStructure,
    global_filetype: Option<&str>,
    base_dir: &Path,
) -> Result<Target, String> {
    let path = resolve(base_dir, &s.filepath);
    let frame = io::read_frame(&path, global_filetype)?;
    let mut target = Target::new(frame, s.number);

    for r in &s.mol_restraints {
        target = apply_mol_restraint(target, r);
    }

    for group in &s.atom_groups {
        target = apply_atom_group(target, group)?;
    }

    if s.center {
        target = target.with_centering(CenteringMode::Center);
    }

    if let Some((pos, euler)) = s.fixed {
        target = target.fixed_at(pos).with_orientation([
            Angle::from_radians(euler[0]),
            Angle::from_radians(euler[1]),
            Angle::from_radians(euler[2]),
        ]);
    }

    Ok(target)
}

fn apply_mol_restraint(target: Target, r: &ParsedRestraint) -> Target {
    match *r {
        ParsedRestraint::InsideBox { min, max } => {
            target.with_restraint(InsideBoxRestraint::new(min, max, [false; 3]))
        }
        ParsedRestraint::InsideSphere { center, radius } => {
            target.with_restraint(InsideSphereRestraint::new(center, radius))
        }
        ParsedRestraint::OutsideSphere { center, radius } => {
            target.with_restraint(OutsideSphereRestraint::new(center, radius))
        }
        ParsedRestraint::AbovePlane { normal, distance } => {
            target.with_restraint(AbovePlaneRestraint::new(normal, distance))
        }
        ParsedRestraint::BelowPlane { normal, distance } => {
            target.with_restraint(BelowPlaneRestraint::new(normal, distance))
        }
    }
}

fn apply_atom_group(mut target: Target, group: &ParsedAtomGroup) -> Result<Target, String> {
    // Parsed atom indices follow Packmol's 1-based convention in the
    // `.inp` file; Target::with_atom_restraint expects 0-based.
    let zero_indexed: Vec<usize> = group
        .atom_indices
        .iter()
        .map(|&i| i.saturating_sub(1))
        .collect();
    for r in &group.restraints {
        let indices = zero_indexed.as_slice();
        target = match *r {
            ParsedRestraint::InsideBox { min, max } => {
                target.with_atom_restraint(indices, InsideBoxRestraint::new(min, max, [false; 3]))
            }
            ParsedRestraint::InsideSphere { center, radius } => {
                target.with_atom_restraint(indices, InsideSphereRestraint::new(center, radius))
            }
            ParsedRestraint::OutsideSphere { center, radius } => {
                target.with_atom_restraint(indices, OutsideSphereRestraint::new(center, radius))
            }
            ParsedRestraint::AbovePlane { normal, distance } => {
                target.with_atom_restraint(indices, AbovePlaneRestraint::new(normal, distance))
            }
            ParsedRestraint::BelowPlane { normal, distance } => {
                target.with_atom_restraint(indices, BelowPlaneRestraint::new(normal, distance))
            }
        };
    }
    Ok(target)
}
