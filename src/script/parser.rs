//! Parser for molpack's line-oriented `.inp` script format.
//!
//! The format mirrors the keyword-driven input Packmol has used for two
//! decades, but the parser and types here belong to molpack and are
//! exposed as first-class library API. Unrecognised keywords are skipped
//! so scripts written against newer features remain readable.

use std::path::PathBuf;

use super::error::ScriptError;

// ──────────────────────────────────────────────────────────────────────────────
// Public AST
// ──────────────────────────────────────────────────────────────────────────────

/// A fully parsed packing script.
#[derive(Debug, Clone)]
pub struct Script {
    /// Minimum atom-atom distance tolerance (`tolerance`). Default 2.0 Å.
    pub tolerance: f64,
    /// Random seed (`seed` keyword).
    pub seed: Option<u64>,
    /// Global format hint (`filetype` keyword). When set, it overrides
    /// per-file extension inference; otherwise each input file's format
    /// is inferred from its extension.
    pub filetype: Option<String>,
    /// Output file path (`output` keyword). Required.
    pub output: PathBuf,
    /// Maximum outer-loop iterations (`nloop` keyword). Default 400.
    pub nloop: usize,
    /// Whether to enforce minimum-distance overlap avoidance
    /// (`avoid_overlap`). Parsed for compatibility; not yet wired to the
    /// Molpack API.
    #[allow(dead_code)]
    pub avoid_overlap: bool,
    pub structures: Vec<Structure>,
}

/// One `structure … end structure` block.
#[derive(Debug, Clone)]
pub struct Structure {
    /// Template molecule file path, as written in the script (may be relative).
    pub filepath: PathBuf,
    /// Number of copies to pack.
    pub number: usize,
    /// Molecule-wide restraints (applied to every atom of every copy).
    pub mol_restraints: Vec<RestraintSpec>,
    /// Atom-subset restraints (`atoms … end atoms` blocks).
    pub atom_groups: Vec<AtomGroup>,
    /// Whether the `center` keyword was present.
    pub center: bool,
    /// Fixed placement: `(position [x,y,z], euler [ex,ey,ez])`.
    pub fixed: Option<([f64; 3], [f64; 3])>,
}

/// One `atoms … end atoms` sub-block.
#[derive(Debug, Clone)]
pub struct AtomGroup {
    /// Atom indices as written in the script (1-based).
    pub atom_indices: Vec<usize>,
    pub restraints: Vec<RestraintSpec>,
}

/// Restraint as it appears in the script, before being mapped to a
/// concrete `Restraint` implementation.
#[derive(Debug, Clone)]
pub enum RestraintSpec {
    InsideBox { min: [f64; 3], max: [f64; 3] },
    InsideSphere { center: [f64; 3], radius: f64 },
    OutsideSphere { center: [f64; 3], radius: f64 },
    /// `over plane` — atom must lie above the plane.
    AbovePlane { normal: [f64; 3], distance: f64 },
    /// `below plane` — atom must lie below the plane.
    BelowPlane { normal: [f64; 3], distance: f64 },
}

// ──────────────────────────────────────────────────────────────────────────────
// Parser
// ──────────────────────────────────────────────────────────────────────────────

/// Parse a script from its textual source.
pub fn parse(src: &str) -> Result<Script, ScriptError> {
    let mut tolerance = 2.0_f64;
    let mut seed: Option<u64> = None;
    let mut filetype: Option<String> = None;
    let mut output: Option<PathBuf> = None;
    let mut nloop: usize = 400;
    let mut avoid_overlap = true;
    let mut structures: Vec<Structure> = Vec::new();

    enum State {
        TopLevel,
        InStructure(Structure),
        InAtoms { structure: Structure, group: AtomGroup },
    }

    let mut state = State::TopLevel;

    for (lineno, raw) in src.lines().enumerate() {
        let lineno = lineno + 1;

        // Strip inline comments and surrounding whitespace.
        let line = match raw.find('#') {
            Some(pos) => &raw[..pos],
            None => raw,
        }
        .trim();

        if line.is_empty() {
            continue;
        }

        let tokens: Vec<&str> = line.split_whitespace().collect();
        let keyword = tokens[0].to_ascii_lowercase();

        state = match state {
            // ── Top-level keywords ──────────────────────────────────────────
            State::TopLevel => match keyword.as_str() {
                "tolerance" => {
                    tolerance = parse_f64(&tokens, 1, "tolerance", lineno)?;
                    State::TopLevel
                }
                "seed" => {
                    seed = Some(parse_u64(&tokens, 1, "seed", lineno)?);
                    State::TopLevel
                }
                "filetype" => {
                    filetype = Some(
                        tokens
                            .get(1)
                            .ok_or_else(|| parse_err(lineno, "missing filetype value"))?
                            .to_ascii_lowercase(),
                    );
                    State::TopLevel
                }
                "output" => {
                    output = Some(PathBuf::from(
                        tokens
                            .get(1)
                            .ok_or_else(|| parse_err(lineno, "missing output path"))?,
                    ));
                    State::TopLevel
                }
                "nloop" => {
                    nloop = parse_usize(&tokens, 1, "nloop", lineno)?;
                    State::TopLevel
                }
                "avoid_overlap" => {
                    let val = tokens
                        .get(1)
                        .map(|s| s.to_ascii_lowercase())
                        .unwrap_or_default();
                    avoid_overlap = val != "no" && val != "false" && val != "0";
                    State::TopLevel
                }
                "structure" => {
                    let path = tokens
                        .get(1)
                        .ok_or_else(|| parse_err(lineno, "`structure` requires a file path"))?;
                    State::InStructure(Structure {
                        filepath: PathBuf::from(path),
                        number: 0,
                        mol_restraints: Vec::new(),
                        atom_groups: Vec::new(),
                        center: false,
                        fixed: None,
                    })
                }
                _ => State::TopLevel, // unknown top-level keyword — skip
            },

            // ── Inside a structure block ────────────────────────────────────
            State::InStructure(mut s) => match keyword.as_str() {
                "end"
                    if tokens.get(1).map(|t| t.to_ascii_lowercase()).as_deref()
                        == Some("structure") =>
                {
                    structures.push(s);
                    State::TopLevel
                }
                "number" => {
                    s.number = parse_usize(&tokens, 1, "number", lineno)?;
                    State::InStructure(s)
                }
                "center" => {
                    s.center = true;
                    State::InStructure(s)
                }
                "fixed" => {
                    let pos = parse_vec3(&tokens, 1, "fixed position", lineno)?;
                    let euler = parse_vec3(&tokens, 4, "fixed euler", lineno)?;
                    s.fixed = Some((pos, euler));
                    State::InStructure(s)
                }
                "inside" => {
                    let r = parse_inside(&tokens, lineno)?;
                    s.mol_restraints.push(r);
                    State::InStructure(s)
                }
                "outside" => {
                    let r = parse_outside(&tokens, lineno)?;
                    s.mol_restraints.push(r);
                    State::InStructure(s)
                }
                "over" | "above" => {
                    let r = parse_plane_above(&tokens, lineno)?;
                    s.mol_restraints.push(r);
                    State::InStructure(s)
                }
                "below" => {
                    let r = parse_plane_below(&tokens, lineno)?;
                    s.mol_restraints.push(r);
                    State::InStructure(s)
                }
                "atoms" => {
                    let indices = tokens[1..]
                        .iter()
                        .map(|t| {
                            t.parse::<usize>()
                                .map_err(|_| parse_err(lineno, format!("invalid atom index `{t}`")))
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    if indices.is_empty() {
                        return Err(parse_err(lineno, "`atoms` requires at least one index"));
                    }
                    State::InAtoms {
                        structure: s,
                        group: AtomGroup {
                            atom_indices: indices,
                            restraints: Vec::new(),
                        },
                    }
                }
                _ => State::InStructure(s),
            },

            // ── Inside an atoms sub-block ───────────────────────────────────
            State::InAtoms {
                structure: s,
                mut group,
            } => match keyword.as_str() {
                "end"
                    if tokens.get(1).map(|t| t.to_ascii_lowercase()).as_deref()
                        == Some("atoms") =>
                {
                    let mut s = s;
                    s.atom_groups.push(group);
                    State::InStructure(s)
                }
                "inside" => {
                    let r = parse_inside(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms { structure: s, group }
                }
                "outside" => {
                    let r = parse_outside(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms { structure: s, group }
                }
                "over" | "above" => {
                    let r = parse_plane_above(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms { structure: s, group }
                }
                "below" => {
                    let r = parse_plane_below(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms { structure: s, group }
                }
                _ => State::InAtoms { structure: s, group },
            },
        };
    }

    match state {
        State::InStructure(_) => return Err(ScriptError::UnclosedBlock("structure")),
        State::InAtoms { .. } => return Err(ScriptError::UnclosedBlock("atoms")),
        State::TopLevel => {}
    }

    Ok(Script {
        tolerance,
        seed,
        filetype,
        output: output.ok_or(ScriptError::MissingOutput)?,
        nloop,
        avoid_overlap,
        structures,
    })
}

// ──────────────────────────────────────────────────────────────────────────────
// Restraint helpers
// ──────────────────────────────────────────────────────────────────────────────

fn parse_inside(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let shape = tokens
        .get(1)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, "`inside` requires a shape keyword"))?;
    match shape.as_str() {
        "box" => {
            let min = parse_vec3(tokens, 2, "inside box min", lineno)?;
            let max = parse_vec3(tokens, 5, "inside box max", lineno)?;
            Ok(RestraintSpec::InsideBox { min, max })
        }
        "sphere" => {
            let center = parse_vec3(tokens, 2, "inside sphere center", lineno)?;
            let radius = parse_f64(tokens, 5, "inside sphere radius", lineno)?;
            Ok(RestraintSpec::InsideSphere { center, radius })
        }
        _ => Err(parse_err(
            lineno,
            format!("unsupported `inside` shape `{shape}`"),
        )),
    }
}

fn parse_outside(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let shape = tokens
        .get(1)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, "`outside` requires a shape keyword"))?;
    match shape.as_str() {
        "sphere" => {
            let center = parse_vec3(tokens, 2, "outside sphere center", lineno)?;
            let radius = parse_f64(tokens, 5, "outside sphere radius", lineno)?;
            Ok(RestraintSpec::OutsideSphere { center, radius })
        }
        _ => Err(parse_err(
            lineno,
            format!("unsupported `outside` shape `{shape}`"),
        )),
    }
}

fn parse_plane_above(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let normal = parse_vec3(tokens, 2, "over plane normal", lineno)?;
    let distance = parse_f64(tokens, 5, "over plane distance", lineno)?;
    Ok(RestraintSpec::AbovePlane { normal, distance })
}

fn parse_plane_below(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let normal = parse_vec3(tokens, 2, "below plane normal", lineno)?;
    let distance = parse_f64(tokens, 5, "below plane distance", lineno)?;
    Ok(RestraintSpec::BelowPlane { normal, distance })
}

// ──────────────────────────────────────────────────────────────────────────────
// Numeric helpers
// ──────────────────────────────────────────────────────────────────────────────

fn parse_f64(tokens: &[&str], idx: usize, ctx: &str, lineno: usize) -> Result<f64, ScriptError> {
    let tok = tokens
        .get(idx)
        .ok_or_else(|| parse_err(lineno, format!("`{ctx}` — missing value at position {idx}")))?;
    tok.parse::<f64>()
        .map_err(|_| parse_err(lineno, format!("`{ctx}` — `{tok}` is not a valid number")))
}

fn parse_u64(tokens: &[&str], idx: usize, ctx: &str, lineno: usize) -> Result<u64, ScriptError> {
    let tok = tokens
        .get(idx)
        .ok_or_else(|| parse_err(lineno, format!("`{ctx}` — missing value at position {idx}")))?;
    tok.parse::<u64>()
        .map_err(|_| parse_err(lineno, format!("`{ctx}` — `{tok}` is not a valid integer")))
}

fn parse_usize(tokens: &[&str], idx: usize, ctx: &str, lineno: usize) -> Result<usize, ScriptError> {
    let tok = tokens
        .get(idx)
        .ok_or_else(|| parse_err(lineno, format!("`{ctx}` — missing value at position {idx}")))?;
    tok.parse::<usize>()
        .map_err(|_| parse_err(lineno, format!("`{ctx}` — `{tok}` is not a valid integer")))
}

fn parse_vec3(
    tokens: &[&str],
    start: usize,
    ctx: &str,
    lineno: usize,
) -> Result<[f64; 3], ScriptError> {
    Ok([
        parse_f64(tokens, start, ctx, lineno)?,
        parse_f64(tokens, start + 1, ctx, lineno)?,
        parse_f64(tokens, start + 2, ctx, lineno)?,
    ])
}

fn parse_err(lineno: usize, message: impl Into<String>) -> ScriptError {
    ScriptError::Parse {
        line: lineno,
        message: message.into(),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_mixture() {
        let src = "\
tolerance 2.0
seed 1234567
filetype pdb
output mixture.xyz

structure water.pdb
  number 1000
  inside box 0. 0. 0. 40. 40. 40.
end structure

structure urea.pdb
  number 400
  inside box 0. 0. 0. 40. 40. 40.
end structure
";
        let inp = parse(src).expect("parse failed");
        assert_eq!(inp.tolerance, 2.0);
        assert_eq!(inp.seed, Some(1_234_567));
        assert_eq!(inp.filetype.as_deref(), Some("pdb"));
        assert_eq!(inp.output, PathBuf::from("mixture.xyz"));
        assert_eq!(inp.structures.len(), 2);
        assert_eq!(inp.structures[0].number, 1000);
        assert_eq!(inp.structures[1].number, 400);
        assert_eq!(inp.structures[0].mol_restraints.len(), 1);
    }

    #[test]
    fn parse_bilayer_atoms_blocks() {
        let src = "\
tolerance 2.0
seed 1234567
filetype pdb
output bilayer.xyz

structure palmitoil.pdb
  number 10
  inside box 0. 0. 0. 40. 40. 14.
  atoms 31 32
    below plane 0. 0. 1. 2.
  end atoms
  atoms 1 2
    over plane 0. 0. 1. 12.
  end atoms
end structure
";
        let inp = parse(src).expect("parse failed");
        let s = &inp.structures[0];
        assert_eq!(s.number, 10);
        assert_eq!(s.mol_restraints.len(), 1);
        assert_eq!(s.atom_groups.len(), 2);
        assert_eq!(s.atom_groups[0].atom_indices, vec![31, 32]);
        assert_eq!(s.atom_groups[1].atom_indices, vec![1, 2]);
    }

    #[test]
    fn parse_interface_fixed_center() {
        let src = "\
tolerance 2.0
seed 1234567
filetype pdb
output interface.xyz

structure t3.pdb
  number 1
  center
  fixed 0. 20. 20. 1.57 1.57 1.57
end structure
";
        let inp = parse(src).expect("parse failed");
        let s = &inp.structures[0];
        assert!(s.center);
        let (pos, euler) = s.fixed.unwrap();
        assert_eq!(pos, [0.0, 20.0, 20.0]);
        assert!((euler[0] - 1.57).abs() < 1e-9);
    }

    #[test]
    fn parse_comments_and_blank_lines() {
        let src = "\
# global settings
tolerance 2.0
# seed comment
seed 42
output out.pdb

# molecule block
structure mol.pdb
  number 10 # trailing comment
  inside sphere 0. 0. 0. 20.
end structure
";
        let inp = parse(src).expect("parse failed");
        assert_eq!(inp.seed, Some(42));
        let s = &inp.structures[0];
        assert_eq!(s.number, 10);
        assert!(matches!(s.mol_restraints[0], RestraintSpec::InsideSphere { .. }));
    }

    #[test]
    fn parse_avoid_overlap_no() {
        let src = "\
tolerance 2.0
output out.pdb
avoid_overlap no

structure mol.pdb
  number 1
  inside sphere 0. 0. 0. 50.
end structure
";
        let inp = parse(src).expect("parse failed");
        assert!(!inp.avoid_overlap);
    }

    #[test]
    fn parse_missing_output_errors() {
        let src = "\
tolerance 2.0

structure mol.pdb
  number 1
  inside sphere 0. 0. 0. 50.
end structure
";
        let err = parse(src).expect_err("should fail");
        assert!(matches!(err, ScriptError::MissingOutput));
    }

    #[test]
    fn parse_unclosed_structure_errors() {
        let src = "\
output out.pdb

structure mol.pdb
  number 1
  inside sphere 0. 0. 0. 50.
";
        let err = parse(src).expect_err("should fail");
        assert!(matches!(err, ScriptError::UnclosedBlock("structure")));
    }
}
