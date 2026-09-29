//! Format-aware molecular file readers and writers for the script loader.
//!
//! Covers the formats the underlying molrs I/O crate supports today:
//! `.pdb`, `.xyz`, `.sdf`/`.mol`, `.lammpstrj`, `.data`. Input format can
//! be chosen explicitly via the script's `filetype` keyword; otherwise
//! it is inferred from the file extension.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use molrs::io::data::sdf::SDFReader;
use molrs::io::reader::FrameReader;
use molrs::store::frame::Frame;

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
        "pdb" => molrs::io::data::pdb::read_pdb_frame(path)
            .map_err(|e| io_err(path, format!("reading PDB: {e}"))),

        "xyz" => molrs::io::data::xyz::read_xyz_frame(path)
            .map_err(|e| io_err(path, format!("reading XYZ: {e}"))),

        "sdf" | "mol" => {
            let file = File::open(path).map_err(|e| io_err(path, format!("opening SDF: {e}")))?;
            let mut reader = SDFReader::new(std::io::BufReader::new(file));
            reader
                .read()
                .map_err(|e| io_err(path, format!("reading SDF: {e}")))?
                .ok_or_else(|| io_err(path, "SDF file contains no records"))
        }

        "lammps_dump" | "lammpstrj" => {
            // The template is the first snapshot; the rest are never read.
            molrs::io::trajectory::lammps_dump::open_lammps_dump(path)
                .map_err(|e| io_err(path, format!("opening LAMMPS dump: {e}")))?
                .read()
                .map_err(|e| io_err(path, format!("reading LAMMPS dump: {e}")))?
                .ok_or_else(|| io_err(path, "LAMMPS dump contains no frames"))
        }

        "lammps_data" | "data" => molrs::io::data::lammps_data::read_lammps_data(path)
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
            molrs::io::data::pdb::write_pdb_frame(&mut writer, frame)
                .map_err(|e| io_err(path, format!("writing PDB: {e}")))
        }

        "xyz" => {
            let file =
                File::create(path).map_err(|e| io_err(path, format!("creating XYZ: {e}")))?;
            let mut writer = BufWriter::new(file);
            molrs::io::data::xyz::write_xyz_frame(&mut writer, frame)
                .map_err(|e| io_err(path, format!("writing XYZ: {e}")))
        }

        "lammpstrj" => molrs::io::trajectory::lammps_dump::write_lammps_dump(
            path,
            std::slice::from_ref(frame),
            None,
        )
        .map_err(|e| io_err(path, format!("writing LAMMPS dump: {e}"))),

        other => Err(io_err(path, format!("unsupported output format `{other}`"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molrs::spatial::simbox::SimBox;
    use molrs::store::block::Block;
    use molrs::types::F;
    use ndarray::{Array1, array};

    fn one_atom_at(x: F) -> Frame {
        let mut atoms = Block::new();
        atoms
            .insert("id", Array1::from_vec(vec![1u32]).into_dyn())
            .unwrap();
        atoms
            .insert("x", Array1::from_vec(vec![x]).into_dyn())
            .unwrap();
        atoms
            .insert("y", Array1::from_vec(vec![0.0 as F]).into_dyn())
            .unwrap();
        atoms
            .insert("z", Array1::from_vec(vec![0.0 as F]).into_dyn())
            .unwrap();
        let mut frame = Frame::new();
        frame.insert("atoms", atoms);
        frame.simbox = Some(SimBox::cube(10.0, array![0.0, 0.0, 0.0], [true; 3]).unwrap());
        frame
    }

    /// A multi-frame dump is a template: its first snapshot, nothing else.
    #[test]
    fn a_lammps_dump_reads_its_first_frame() {
        let path =
            std::env::temp_dir().join(format!("molpack-io-{}.lammpstrj", std::process::id()));
        molrs::io::trajectory::lammps_dump::write_lammps_dump(
            &path,
            &[one_atom_at(1.0), one_atom_at(5.0)],
            None,
        )
        .unwrap();
        let frame = read_frame(&path, None);
        let _ = std::fs::remove_file(&path);
        let x = frame
            .unwrap()
            .get("atoms")
            .unwrap()
            .get_float("x")
            .unwrap()
            .iter()
            .copied()
            .collect::<Vec<F>>();
        assert_eq!(x, [1.0]);
    }
}
