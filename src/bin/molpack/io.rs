//! Format-aware molecular file reader and writer.
//!
//! Extends the formats Packmol supports (pdb, xyz) with the additional readers
//! available in molrs-io: SDF/MOL, LAMMPS dump, LAMMPS data.
//!
//! Format resolution order:
//! 1. Explicit `filetype` keyword from the .inp file
//! 2. File extension of the path

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use molrs::frame::Frame;
use molrs_io::reader::FrameReader;
use molrs_io::sdf::SDFReader;

/// Derive a format string from a file extension. Returns `None` if unrecognised.
fn ext_format(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
        .filter(|e| {
            matches!(
                e.as_str(),
                "pdb" | "xyz" | "sdf" | "mol" | "lammpstrj" | "data"
            )
        })
}

/// Read the first frame from `path`.
///
/// `filetype_hint` comes from the `.inp` `filetype` keyword; when `None`, the
/// format is inferred from the file extension.
pub fn read_frame(path: &Path, filetype_hint: Option<&str>) -> Result<Frame, String> {
    let fmt = filetype_hint
        .map(str::to_ascii_lowercase)
        .or_else(|| ext_format(path))
        .ok_or_else(|| {
            format!(
                "cannot determine format for `{}` — set `filetype` in the .inp file or use a \
                 recognised extension (.pdb, .xyz, .sdf, .mol, .lammpstrj, .data)",
                path.display()
            )
        })?;

    match fmt.as_str() {
        "pdb" => molrs_io::pdb::read_pdb_frame(path)
            .map_err(|e| format!("reading PDB `{}`: {e}", path.display())),

        "xyz" => molrs_io::xyz::read_xyz_frame(path)
            .map_err(|e| format!("reading XYZ `{}`: {e}", path.display())),

        "sdf" | "mol" => {
            let file =
                File::open(path).map_err(|e| format!("opening SDF `{}`: {e}", path.display()))?;
            let mut reader = SDFReader::new(std::io::BufReader::new(file));
            reader
                .read_frame()
                .map_err(|e| format!("reading SDF `{}`: {e}", path.display()))?
                .ok_or_else(|| format!("SDF file `{}` contains no records", path.display()))
        }

        "lammps_dump" | "lammpstrj" => {
            let mut frames = molrs_io::lammps_dump::read_lammps_dump(path)
                .map_err(|e| format!("reading LAMMPS dump `{}`: {e}", path.display()))?;
            if frames.is_empty() {
                return Err(format!(
                    "LAMMPS dump `{}` contains no frames",
                    path.display()
                ));
            }
            Ok(frames.swap_remove(0))
        }

        "lammps_data" | "data" => molrs_io::lammps_data::read_lammps_data(path)
            .map_err(|e| format!("reading LAMMPS data `{}`: {e}", path.display())),

        other => Err(format!("unsupported input format `{other}`")),
    }
}

/// Write `frame` to `path`.
///
/// Output format is inferred from the file extension of `path`.
/// Supported output formats: `.pdb`, `.xyz`, `.lammpstrj`.
pub fn write_frame(path: &Path, frame: &Frame) -> Result<(), String> {
    let fmt = ext_format(path).ok_or_else(|| {
        format!(
            "cannot determine output format for `{}` — use a recognised extension \
             (.pdb, .xyz, .lammpstrj)",
            path.display()
        )
    })?;

    match fmt.as_str() {
        "pdb" => {
            let file = File::create(path)
                .map_err(|e| format!("creating PDB `{}`: {e}", path.display()))?;
            let mut writer = BufWriter::new(file);
            molrs_io::pdb::write_pdb_frame(&mut writer, frame)
                .map_err(|e| format!("writing PDB `{}`: {e}", path.display()))
        }

        "xyz" => {
            let file = File::create(path)
                .map_err(|e| format!("creating XYZ `{}`: {e}", path.display()))?;
            let mut writer = BufWriter::new(file);
            molrs_io::xyz::write_xyz_frame(&mut writer, frame)
                .map_err(|e| format!("writing XYZ `{}`: {e}", path.display()))
        }

        "lammpstrj" => molrs_io::lammps_dump::write_lammps_dump(path, std::slice::from_ref(frame))
            .map_err(|e| format!("writing LAMMPS dump `{}`: {e}", path.display())),

        other => Err(format!("unsupported output format `{other}`")),
    }
}
