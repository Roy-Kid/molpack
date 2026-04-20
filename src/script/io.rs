//! Format-aware molecular file readers and writers for the script loader.
//!
//! Covers the formats the underlying molrs I/O crate supports today:
//! `.pdb`, `.xyz`, `.sdf`/`.mol`, `.lammpstrj`, `.data`. Input format can
//! be chosen explicitly via the script's `filetype` keyword; otherwise
//! it is inferred from the file extension.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use molrs::frame::Frame;
use molrs_io::reader::FrameReader;
use molrs_io::sdf::SDFReader;

use super::error::ScriptError;

fn io_err(path: &Path, message: impl Into<String>) -> ScriptError {
    ScriptError::Io {
        path: path.to_path_buf(),
        message: message.into(),
    }
}

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
/// `filetype_hint` comes from the script's `filetype` keyword; when
/// `None`, the format is inferred from the file extension.
pub fn read_frame(path: &Path, filetype_hint: Option<&str>) -> Result<Frame, ScriptError> {
    let fmt = filetype_hint
        .map(str::to_ascii_lowercase)
        .or_else(|| ext_format(path))
        .ok_or_else(|| {
            io_err(
                path,
                "cannot determine format — set `filetype` in the script or use a recognised \
                 extension (.pdb, .xyz, .sdf, .mol, .lammpstrj, .data)",
            )
        })?;

    match fmt.as_str() {
        "pdb" => molrs_io::pdb::read_pdb_frame(path)
            .map_err(|e| io_err(path, format!("reading PDB: {e}"))),

        "xyz" => molrs_io::xyz::read_xyz_frame(path)
            .map_err(|e| io_err(path, format!("reading XYZ: {e}"))),

        "sdf" | "mol" => {
            let file = File::open(path).map_err(|e| io_err(path, format!("opening SDF: {e}")))?;
            let mut reader = SDFReader::new(std::io::BufReader::new(file));
            reader
                .read_frame()
                .map_err(|e| io_err(path, format!("reading SDF: {e}")))?
                .ok_or_else(|| io_err(path, "SDF file contains no records"))
        }

        "lammps_dump" | "lammpstrj" => {
            let mut frames = molrs_io::lammps_dump::read_lammps_dump(path)
                .map_err(|e| io_err(path, format!("reading LAMMPS dump: {e}")))?;
            if frames.is_empty() {
                return Err(io_err(path, "LAMMPS dump contains no frames"));
            }
            Ok(frames.swap_remove(0))
        }

        "lammps_data" | "data" => molrs_io::lammps_data::read_lammps_data(path)
            .map_err(|e| io_err(path, format!("reading LAMMPS data: {e}"))),

        other => Err(io_err(path, format!("unsupported input format `{other}`"))),
    }
}

/// Write `frame` to `path`.
///
/// Output format is inferred from the file extension. Supported formats:
/// `.pdb`, `.xyz`, `.lammpstrj`.
pub fn write_frame(path: &Path, frame: &Frame) -> Result<(), ScriptError> {
    let fmt = ext_format(path).ok_or_else(|| {
        io_err(
            path,
            "cannot determine output format — use a recognised extension (.pdb, .xyz, .lammpstrj)",
        )
    })?;

    match fmt.as_str() {
        "pdb" => {
            let file =
                File::create(path).map_err(|e| io_err(path, format!("creating PDB: {e}")))?;
            let mut writer = BufWriter::new(file);
            molrs_io::pdb::write_pdb_frame(&mut writer, frame)
                .map_err(|e| io_err(path, format!("writing PDB: {e}")))
        }

        "xyz" => {
            let file =
                File::create(path).map_err(|e| io_err(path, format!("creating XYZ: {e}")))?;
            let mut writer = BufWriter::new(file);
            molrs_io::xyz::write_xyz_frame(&mut writer, frame)
                .map_err(|e| io_err(path, format!("writing XYZ: {e}")))
        }

        "lammpstrj" => molrs_io::lammps_dump::write_lammps_dump(path, std::slice::from_ref(frame))
            .map_err(|e| io_err(path, format!("writing LAMMPS dump: {e}"))),

        other => Err(io_err(path, format!("unsupported output format `{other}`"))),
    }
}
