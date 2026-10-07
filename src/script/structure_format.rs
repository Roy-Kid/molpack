//! The `.inp` `filetype` vocabulary: which structure-file format a template
//! or the output file is in.
//!
//! molrs names every file door after its format (`read_pdb`, `write_xyz`, …)
//! and picks no format for the caller, so the choice the `.inp` grammar makes
//! — the script's `filetype`, else the file name — is molpack's, and lives
//! here once. Reading and writing a format is molrs's (feature `io`); the
//! names, extensions and resolution are always compiled, so an embedding host
//! without `io` (the PyO3 wheel) resolves a template's format by the same
//! rule.

use std::fmt;
use std::path::Path;

use crate::script::ScriptError;

/// A structure-file format the `.inp` `filetype` keyword (or a file
/// extension) can name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StructureFormat {
    Pdb,
    Xyz,
    Sdf,
    Mol2,
    Gro,
    Cif,
    VaspPoscar,
    Xsf,
    Cube,
    AmberInpcrd,
    LammpsData,
    LammpsDump,
}

impl StructureFormat {
    /// Every format, in the order error messages list them.
    pub const ALL: [StructureFormat; 12] = [
        StructureFormat::Pdb,
        StructureFormat::Xyz,
        StructureFormat::Sdf,
        StructureFormat::Mol2,
        StructureFormat::Gro,
        StructureFormat::Cif,
        StructureFormat::VaspPoscar,
        StructureFormat::Xsf,
        StructureFormat::Cube,
        StructureFormat::AmberInpcrd,
        StructureFormat::LammpsData,
        StructureFormat::LammpsDump,
    ];

    /// The `filetype` value naming this format.
    pub fn name(self) -> &'static str {
        match self {
            StructureFormat::Pdb => "pdb",
            StructureFormat::Xyz => "xyz",
            StructureFormat::Sdf => "sdf",
            StructureFormat::Mol2 => "mol2",
            StructureFormat::Gro => "gro",
            StructureFormat::Cif => "cif",
            StructureFormat::VaspPoscar => "poscar",
            StructureFormat::Xsf => "xsf",
            StructureFormat::Cube => "cube",
            StructureFormat::AmberInpcrd => "inpcrd",
            StructureFormat::LammpsData => "lammps_data",
            StructureFormat::LammpsDump => "lammps_dump",
        }
    }

    /// The file extensions (lower case, no dot) that name this format.
    pub fn extensions(self) -> &'static [&'static str] {
        match self {
            StructureFormat::Pdb => &["pdb", "ent"],
            StructureFormat::Xyz => &["xyz", "extxyz"],
            StructureFormat::Sdf => &["sdf", "mol"],
            StructureFormat::Mol2 => &["mol2"],
            StructureFormat::Gro => &["gro"],
            StructureFormat::Cif => &["cif"],
            StructureFormat::VaspPoscar => &["poscar", "vasp"],
            StructureFormat::Xsf => &["xsf"],
            StructureFormat::Cube => &["cube", "cub"],
            StructureFormat::AmberInpcrd => &["inpcrd", "rst7", "restrt", "crd"],
            StructureFormat::LammpsData => &["data", "lmp"],
            StructureFormat::LammpsDump => &["lammpstrj", "dump"],
        }
    }

    /// Whether molrs writes this format (SDF and inpcrd are read-only).
    pub fn is_writable(self) -> bool {
        !matches!(self, StructureFormat::Sdf | StructureFormat::AmberInpcrd)
    }

    /// The format a `filetype` value names: its name or one of its
    /// extensions, case-insensitive, a leading dot ignored.
    pub fn from_name(name: &str) -> Option<StructureFormat> {
        let name = name.trim().trim_start_matches('.').to_ascii_lowercase();
        Self::ALL
            .into_iter()
            .find(|f| f.name() == name || f.extensions().contains(&name.as_str()))
    }

    /// The format a file name names: its extension, else a VASP
    /// `POSCAR*` / `CONTCAR*` stem.
    pub fn from_path(path: &Path) -> Option<StructureFormat> {
        let by_extension = path
            .extension()
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
            .and_then(|e| {
                Self::ALL
                    .into_iter()
                    .find(|f| f.extensions().contains(&e.as_str()))
            });
        if by_extension.is_some() {
            return by_extension;
        }
        let stem = path.file_name()?.to_str()?.to_ascii_uppercase();
        (stem.starts_with("POSCAR") || stem.starts_with("CONTCAR"))
            .then_some(StructureFormat::VaspPoscar)
    }

    /// The format of the file at `path`: `filetype` when the script states
    /// one, else the file name. Neither naming a known format is a
    /// [`ScriptError::Io`] listing the known names.
    pub fn resolve(path: &Path, filetype: Option<&str>) -> Result<StructureFormat, ScriptError> {
        let found = match filetype {
            Some(name) => Self::from_name(name),
            None => Self::from_path(path),
        };
        found.ok_or_else(|| {
            let known = Self::ALL.map(StructureFormat::name).join(", ");
            let message = match filetype {
                Some(name) => format!("unknown filetype {name:?}; known: {known}"),
                None => format!(
                    "cannot tell the format from the file name; state a filetype (one of {known})"
                ),
            };
            ScriptError::Io {
                path: path.to_path_buf(),
                message,
            }
        })
    }

    /// Read the first structure of the file at `path` through this format's
    /// molrs reader.
    #[cfg(feature = "io")]
    pub fn read(self, path: &Path) -> Result<molrs::core::Frame, ScriptError> {
        use molrs::io;
        let io_err = |message: String| ScriptError::Io {
            path: path.to_path_buf(),
            message,
        };
        let read = match self {
            StructureFormat::Pdb => io::read_pdb(path).map_err(|e| e.to_string()),
            StructureFormat::Xyz => io::read_xyz(path).map_err(|e| e.to_string()),
            StructureFormat::Sdf => io::read_sdf(path).map_err(|e| e.to_string()),
            StructureFormat::Mol2 => io::read_mol2(path).map_err(|e| e.to_string()),
            StructureFormat::Gro => io::read_gro(path).map_err(|e| e.to_string()),
            StructureFormat::Cif => io::read_cif(path).map_err(|e| e.to_string()),
            StructureFormat::VaspPoscar => io::read_vasp_poscar(path).map_err(|e| e.to_string()),
            StructureFormat::Xsf => io::read_xsf(path).map_err(|e| e.to_string()),
            StructureFormat::Cube => io::read_cube(path).map_err(|e| e.to_string()),
            StructureFormat::AmberInpcrd => io::read_amber_inpcrd(path).map_err(|e| e.to_string()),
            StructureFormat::LammpsData => io::read_lammps_data(path).map_err(|e| e.to_string()),
            StructureFormat::LammpsDump => io::read_lammps_dump_trajectory(path)
                .map_err(|e| e.to_string())
                .and_then(|frames| {
                    frames
                        .into_iter()
                        .next()
                        .ok_or_else(|| "the lammps_dump file holds no frame".to_string())
                }),
        };
        read.map_err(|e| io_err(format!("reading {self}: {e}")))
    }

    /// Write `frame` to `path` through this format's molrs writer.
    #[cfg(feature = "io")]
    pub fn write(self, path: &Path, frame: &molrs::core::Frame) -> Result<(), ScriptError> {
        use molrs::io;
        let written = match self {
            StructureFormat::Pdb => io::write_pdb(path, frame).map_err(|e| e.to_string()),
            StructureFormat::Xyz => io::write_xyz(path, frame).map_err(|e| e.to_string()),
            StructureFormat::Mol2 => io::write_mol2(path, frame).map_err(|e| e.to_string()),
            StructureFormat::Gro => io::write_gro(path, frame).map_err(|e| e.to_string()),
            StructureFormat::Cif => io::write_cif(path, frame).map_err(|e| e.to_string()),
            StructureFormat::VaspPoscar => {
                io::write_vasp_poscar(path, frame).map_err(|e| e.to_string())
            }
            StructureFormat::Xsf => io::write_xsf(path, frame).map_err(|e| e.to_string()),
            StructureFormat::Cube => io::write_cube(path, frame).map_err(|e| e.to_string()),
            StructureFormat::LammpsData => {
                io::write_lammps_data(path, frame).map_err(|e| e.to_string())
            }
            StructureFormat::LammpsDump => {
                io::write_lammps_dump_trajectory(path, std::slice::from_ref(frame), None)
                    .map_err(|e| e.to_string())
            }
            StructureFormat::Sdf | StructureFormat::AmberInpcrd => {
                Err("molrs has no writer for this format".to_string())
            }
        };
        written.map_err(|e| ScriptError::Io {
            path: path.to_path_buf(),
            message: format!("writing {self}: {e}"),
        })
    }
}

impl fmt::Display for StructureFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

#[cfg(test)]
mod tests {
    use super::StructureFormat;
    use std::path::Path;

    #[test]
    fn filetype_wins_over_the_file_name() {
        let format = StructureFormat::resolve(Path::new("water.xyz"), Some("pdb")).unwrap();
        assert_eq!(format, StructureFormat::Pdb);
    }

    #[test]
    fn the_file_name_names_the_format_without_a_filetype() {
        let resolve = |name: &str| StructureFormat::resolve(Path::new(name), None).unwrap();
        assert_eq!(resolve("a.PDB"), StructureFormat::Pdb);
        assert_eq!(resolve("conf.gro"), StructureFormat::Gro);
        assert_eq!(resolve("system.data"), StructureFormat::LammpsData);
        assert_eq!(resolve("CONTCAR"), StructureFormat::VaspPoscar);
    }

    #[test]
    fn every_name_round_trips() {
        for format in StructureFormat::ALL {
            assert_eq!(StructureFormat::from_name(format.name()), Some(format));
        }
    }

    #[test]
    fn an_unknown_format_is_named_in_the_error() {
        let err = StructureFormat::resolve(Path::new("a.unknown"), None).unwrap_err();
        assert!(err.to_string().contains("state a filetype"), "{err}");
        let err = StructureFormat::resolve(Path::new("a.pdb"), Some("tinker")).unwrap_err();
        assert!(err.to_string().contains("\"tinker\""), "{err}");
    }
}
