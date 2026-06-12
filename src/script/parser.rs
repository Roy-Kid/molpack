//! Parser for molpack's line-oriented `.inp` script format.
//!
//! The format mirrors the keyword-driven input Packmol has used for two
//! decades, but the parser and types here belong to molpack and are
//! exposed as first-class library API. Unrecognised keywords are a
//! hard error: silently dropping them risks wrong semantics (e.g., an
//! ignored `pbc` directive that leaves the packer guessing a box size
//! from the initial random placement, blowing the cell grid to
//! gigabytes of RAM).

use std::path::PathBuf;

use super::error::ScriptError;
use super::parser_profile::{
    ProfileDistribution, ProfileGeometry, ProfileInputKind, parse_profile,
};

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
    /// Maximum outer-loop iterations (`nloop` keyword). When unset, Packmol's
    /// default of `200 * ntype` is used (getinp.f90:537-539).
    pub nloop: usize,
    /// Whether to reject initial random placements that overlap a fixed
    /// molecule (`avoid_overlap`, default on). Wired through `Script::build`
    /// to `Molpack::with_avoid_overlap`.
    pub avoid_overlap: bool,
    /// Periodic-boundary box (`pbc` keyword). When set, it seeds the
    /// packer's cell grid so the initial ±`sidemax` random placement
    /// never drives `ncells` to anything astronomical.
    pub pbc: Option<PbcSpec>,
    pub structures: Vec<Structure>,
}

/// Periodic-boundary box declared at the top level via `pbc`.
///
/// Packmol accepts two forms (`getinp.f90`):
/// - `pbc X Y Z`                           → `min = [0,0,0]`, `max = [X,Y,Z]`
/// - `pbc X0 Y0 Z0  X1 Y1 Z1`              → `min = [X0,Y0,Z0]`, `max = [X1,Y1,Z1]`
///
/// Both forms declare periodicity on every axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PbcSpec {
    pub min: [f64; 3],
    pub max: [f64; 3],
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
    InsideBox {
        min: [f64; 3],
        max: [f64; 3],
    },
    InsideSphere {
        center: [f64; 3],
        radius: f64,
    },
    OutsideSphere {
        center: [f64; 3],
        radius: f64,
    },
    /// `over plane` — atom must lie above the plane.
    AbovePlane {
        normal: [f64; 3],
        distance: f64,
    },
    /// `below plane` — atom must lie below the plane.
    BelowPlane {
        normal: [f64; 3],
        distance: f64,
    },
    /// `inside cube x0 y0 z0 d` — atom inside an axis-aligned cube.
    InsideCube {
        origin: [f64; 3],
        side: f64,
    },
    /// `outside cube x0 y0 z0 d` — atom outside an axis-aligned cube.
    OutsideCube {
        origin: [f64; 3],
        side: f64,
    },
    /// `outside box xmin ymin zmin xmax ymax zmax` — atom outside a box.
    OutsideBox {
        min: [f64; 3],
        max: [f64; 3],
    },
    /// `inside ellipsoid a1 a2 a3 b1 b2 b3 c` — center, semi-axes, exponent.
    InsideEllipsoid {
        center: [f64; 3],
        axes: [f64; 3],
        exponent: f64,
    },
    /// `outside ellipsoid a1 a2 a3 b1 b2 b3 c` — center, semi-axes, exponent.
    OutsideEllipsoid {
        center: [f64; 3],
        axes: [f64; 3],
        exponent: f64,
    },
    /// `inside cylinder a1 a2 a3 d1 d2 d3 r l` — center, axis, radius, length.
    InsideCylinder {
        center: [f64; 3],
        axis: [f64; 3],
        radius: f64,
        length: f64,
    },
    /// `outside cylinder a1 a2 a3 d1 d2 d3 r l` — center, axis, radius, length.
    OutsideCylinder {
        center: [f64; 3],
        axis: [f64; 3],
        radius: f64,
        length: f64,
    },
    /// `profile <distribution> <geometry> <params> [density|histogram]` — bias
    /// the selected sites toward a target 1-D spatial distribution. Lowered to a
    /// [`ProfileRestraint`](crate::restraint::profile::ProfileRestraint).
    Profile {
        /// Target-distribution shape (gaussian / erf / tanh / exponential /
        /// tabulated) with its own numeric parameters.
        dist: ProfileDistribution,
        /// Reaction-coordinate geometry (plane / radial / cylinder).
        geometry: ProfileGeometry,
        /// Whether a tabulated input is a density or a count histogram
        /// (default histogram); ignored by the analytic shapes.
        input_kind: ProfileInputKind,
    },
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
    let mut nloop: Option<usize> = None;
    let mut avoid_overlap = true;
    let mut pbc: Option<PbcSpec> = None;
    let mut structures: Vec<Structure> = Vec::new();

    enum State {
        TopLevel,
        InStructure(Structure),
        InAtoms {
            structure: Structure,
            group: AtomGroup,
        },
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
                    nloop = Some(parse_usize(&tokens, 1, "nloop", lineno)?);
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
                "pbc" => {
                    pbc = Some(parse_pbc(&tokens, lineno)?);
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
                _ => {
                    return Err(unknown_keyword(lineno, &keyword, "top-level"));
                }
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
                "center" | "centerofmass" => {
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
                "profile" => {
                    let r = parse_profile(&tokens, lineno)?;
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
                _ => {
                    return Err(unknown_keyword(lineno, &keyword, "structure block"));
                }
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
                    State::InAtoms {
                        structure: s,
                        group,
                    }
                }
                "outside" => {
                    let r = parse_outside(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms {
                        structure: s,
                        group,
                    }
                }
                "over" | "above" => {
                    let r = parse_plane_above(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms {
                        structure: s,
                        group,
                    }
                }
                "below" => {
                    let r = parse_plane_below(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms {
                        structure: s,
                        group,
                    }
                }
                "profile" => {
                    let r = parse_profile(&tokens, lineno)?;
                    group.restraints.push(r);
                    State::InAtoms {
                        structure: s,
                        group,
                    }
                }
                _ => {
                    return Err(unknown_keyword(lineno, &keyword, "atoms block"));
                }
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
        // Packmol default (getinp.f90:537-539): unset `nloop` resolves to
        // 200 * ntype, where ntype is the number of structure types.
        nloop: nloop.unwrap_or(200 * structures.len()),
        avoid_overlap,
        pbc,
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
        "cube" => {
            let (origin, side) = parse_cube_args(tokens, lineno, "inside cube")?;
            Ok(RestraintSpec::InsideCube { origin, side })
        }
        "sphere" => {
            let center = parse_vec3(tokens, 2, "inside sphere center", lineno)?;
            let radius = parse_f64(tokens, 5, "inside sphere radius", lineno)?;
            Ok(RestraintSpec::InsideSphere { center, radius })
        }
        "ellipsoid" => {
            let (center, axes, exponent) =
                parse_ellipsoid_args(tokens, lineno, "inside ellipsoid")?;
            Ok(RestraintSpec::InsideEllipsoid {
                center,
                axes,
                exponent,
            })
        }
        "cylinder" => {
            let (center, axis, radius, length) =
                parse_cylinder_args(tokens, lineno, "inside cylinder")?;
            Ok(RestraintSpec::InsideCylinder {
                center,
                axis,
                radius,
                length,
            })
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
        "box" => {
            let min = parse_vec3(tokens, 2, "outside box min", lineno)?;
            let max = parse_vec3(tokens, 5, "outside box max", lineno)?;
            Ok(RestraintSpec::OutsideBox { min, max })
        }
        "cube" => {
            let (origin, side) = parse_cube_args(tokens, lineno, "outside cube")?;
            Ok(RestraintSpec::OutsideCube { origin, side })
        }
        "sphere" => {
            let center = parse_vec3(tokens, 2, "outside sphere center", lineno)?;
            let radius = parse_f64(tokens, 5, "outside sphere radius", lineno)?;
            Ok(RestraintSpec::OutsideSphere { center, radius })
        }
        "ellipsoid" => {
            let (center, axes, exponent) =
                parse_ellipsoid_args(tokens, lineno, "outside ellipsoid")?;
            Ok(RestraintSpec::OutsideEllipsoid {
                center,
                axes,
                exponent,
            })
        }
        "cylinder" => {
            let (center, axis, radius, length) =
                parse_cylinder_args(tokens, lineno, "outside cylinder")?;
            Ok(RestraintSpec::OutsideCylinder {
                center,
                axis,
                radius,
                length,
            })
        }
        _ => Err(parse_err(
            lineno,
            format!("unsupported `outside` shape `{shape}`"),
        )),
    }
}

/// `<kw> cube x0 y0 z0 d` — origin corner + side length.
fn parse_cube_args(
    tokens: &[&str],
    lineno: usize,
    ctx: &str,
) -> Result<([f64; 3], f64), ScriptError> {
    let origin = parse_vec3(tokens, 2, ctx, lineno)?;
    let side = parse_f64(tokens, 5, ctx, lineno)?;
    Ok((origin, side))
}

/// `<kw> ellipsoid a1 a2 a3 b1 b2 b3 c` — center, semi-axes, exponent.
fn parse_ellipsoid_args(
    tokens: &[&str],
    lineno: usize,
    ctx: &str,
) -> Result<([f64; 3], [f64; 3], f64), ScriptError> {
    let center = parse_vec3(tokens, 2, ctx, lineno)?;
    let axes = parse_vec3(tokens, 5, ctx, lineno)?;
    if axes.iter().any(|&a| a <= 0.0) {
        return Err(parse_err(
            lineno,
            format!("{ctx} semi-axes must be strictly positive, got {axes:?}"),
        ));
    }
    let exponent = parse_f64(tokens, 8, ctx, lineno)?;
    Ok((center, axes, exponent))
}

/// `<kw> cylinder a1 a2 a3 d1 d2 d3 r l` — center, axis, radius, length.
fn parse_cylinder_args(
    tokens: &[&str],
    lineno: usize,
    ctx: &str,
) -> Result<([f64; 3], [f64; 3], f64, f64), ScriptError> {
    let center = parse_vec3(tokens, 2, ctx, lineno)?;
    let axis = parse_vec3(tokens, 5, ctx, lineno)?;
    if axis == [0.0; 3] {
        return Err(parse_err(
            lineno,
            format!("{ctx} axis must be a non-zero direction"),
        ));
    }
    let radius = parse_f64(tokens, 8, ctx, lineno)?;
    let length = parse_f64(tokens, 9, ctx, lineno)?;
    Ok((center, axis, radius, length))
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

/// Parse `pbc` in either the 3-value (`pbc X Y Z`) or 6-value
/// (`pbc X0 Y0 Z0  X1 Y1 Z1`) form — matching packmol `getinp.f90`
/// lines 175-191.
fn parse_pbc(tokens: &[&str], lineno: usize) -> Result<PbcSpec, ScriptError> {
    // tokens[0] == "pbc"; the remaining tokens carry the numeric payload.
    match tokens.len() - 1 {
        3 => {
            let size = parse_vec3(tokens, 1, "pbc box size", lineno)?;
            if size.iter().any(|&v| v <= 0.0) {
                return Err(parse_err(
                    lineno,
                    format!("`pbc` box length must be positive on every axis; got {size:?}"),
                ));
            }
            Ok(PbcSpec {
                min: [0.0; 3],
                max: size,
            })
        }
        6 => {
            let min = parse_vec3(tokens, 1, "pbc min", lineno)?;
            let max = parse_vec3(tokens, 4, "pbc max", lineno)?;
            for k in 0..3 {
                if max[k] <= min[k] {
                    return Err(parse_err(
                        lineno,
                        format!(
                            "`pbc` max must exceed min on every axis; axis {k}: min={}, max={}",
                            min[k], max[k]
                        ),
                    ));
                }
            }
            Ok(PbcSpec { min, max })
        }
        n => Err(parse_err(
            lineno,
            format!("`pbc` expects 3 or 6 numeric values, got {n}"),
        )),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Numeric helpers
// ──────────────────────────────────────────────────────────────────────────────

pub(super) fn parse_f64(
    tokens: &[&str],
    idx: usize,
    ctx: &str,
    lineno: usize,
) -> Result<f64, ScriptError> {
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

fn parse_usize(
    tokens: &[&str],
    idx: usize,
    ctx: &str,
    lineno: usize,
) -> Result<usize, ScriptError> {
    let tok = tokens
        .get(idx)
        .ok_or_else(|| parse_err(lineno, format!("`{ctx}` — missing value at position {idx}")))?;
    tok.parse::<usize>()
        .map_err(|_| parse_err(lineno, format!("`{ctx}` — `{tok}` is not a valid integer")))
}

pub(super) fn parse_vec3(
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

pub(super) fn parse_err(lineno: usize, message: impl Into<String>) -> ScriptError {
    ScriptError::Parse {
        line: lineno,
        message: message.into(),
    }
}

fn unknown_keyword(lineno: usize, keyword: &str, context: &'static str) -> ScriptError {
    ScriptError::UnknownKeyword {
        line: lineno,
        keyword: keyword.to_string(),
        context,
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
    fn parse_centerofmass_alias() {
        // Official Packmol scripts (interface, solvprotein) spell the `center`
        // keyword as `centerofmass`; the parser must treat them identically.
        let src = "\
tolerance 2.0
filetype xyz
output interface.xyz

structure t3.xyz
  number 1
  centerofmass
  fixed 0. 20. 20. 1.57 1.57 1.57
end structure
";
        let inp = parse(src).expect("parse failed");
        assert!(inp.structures[0].center);
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
        assert!(matches!(
            s.mol_restraints[0],
            RestraintSpec::InsideSphere { .. }
        ));
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

    #[test]
    fn parse_pbc_three_values() {
        let src = "\
tolerance 2.0
output out.pdb
pbc 100.0 80.0 60.0

structure mol.pdb
  number 1
end structure
";
        let inp = parse(src).expect("parse failed");
        let pbc = inp.pbc.expect("pbc should be set");
        assert_eq!(pbc.min, [0.0, 0.0, 0.0]);
        assert_eq!(pbc.max, [100.0, 80.0, 60.0]);
    }

    #[test]
    fn parse_pbc_six_values() {
        let src = "\
tolerance 2.0
output out.pdb
pbc -50.0 -40.0 -30.0  50.0 40.0 30.0

structure mol.pdb
  number 1
end structure
";
        let inp = parse(src).expect("parse failed");
        let pbc = inp.pbc.expect("pbc should be set");
        assert_eq!(pbc.min, [-50.0, -40.0, -30.0]);
        assert_eq!(pbc.max, [50.0, 40.0, 30.0]);
    }

    #[test]
    fn parse_pbc_rejects_bad_arity() {
        let src = "\
output out.pdb
pbc 100.0 80.0

structure mol.pdb
  number 1
end structure
";
        let err = parse(src).expect_err("should fail");
        assert!(matches!(err, ScriptError::Parse { .. }));
    }

    #[test]
    fn parse_pbc_rejects_non_positive_extent() {
        let src = "\
output out.pdb
pbc 100.0 0.0 60.0

structure mol.pdb
  number 1
end structure
";
        let err = parse(src).expect_err("should fail");
        assert!(matches!(err, ScriptError::Parse { .. }));
    }

    #[test]
    fn parse_unknown_top_level_keyword_errors() {
        // Regression: parser used to silently drop unknown top-level
        // keywords, which caused `pbc` to be dropped and drove the cell
        // grid to ~10⁸ cells at pack time.
        let src = "\
output out.pdb
wibble 1 2 3

structure mol.pdb
  number 1
end structure
";
        let err = parse(src).expect_err("should fail");
        match err {
            ScriptError::UnknownKeyword {
                keyword, context, ..
            } => {
                assert_eq!(keyword, "wibble");
                assert_eq!(context, "top-level");
            }
            other => panic!("expected UnknownKeyword, got {other:?}"),
        }
    }

    #[test]
    fn parse_unknown_keyword_inside_structure_errors() {
        let src = "\
output out.pdb

structure mol.pdb
  number 1
  nonsense foo
end structure
";
        let err = parse(src).expect_err("should fail");
        match err {
            ScriptError::UnknownKeyword {
                keyword, context, ..
            } => {
                assert_eq!(keyword, "nonsense");
                assert_eq!(context, "structure block");
            }
            other => panic!("expected UnknownKeyword, got {other:?}"),
        }
    }

    /// Parse a single restraint line inside a structure block and return the
    /// parsed `RestraintSpec`.
    fn parse_one_restraint(line: &str) -> RestraintSpec {
        let src =
            format!("output out.pdb\n\nstructure mol.pdb\n  number 1\n  {line}\nend structure\n");
        let inp = parse(&src).expect("parse failed");
        inp.structures[0].mol_restraints[0].clone()
    }

    #[test]
    fn parse_inside_cube() {
        assert!(matches!(
            parse_one_restraint("inside cube 1. 2. 3. 10."),
            RestraintSpec::InsideCube { origin, side }
                if origin == [1.0, 2.0, 3.0] && side == 10.0
        ));
    }

    #[test]
    fn parse_outside_cube() {
        assert!(matches!(
            parse_one_restraint("outside cube 0. 0. 0. 5."),
            RestraintSpec::OutsideCube { side, .. } if side == 5.0
        ));
    }

    #[test]
    fn parse_outside_box() {
        assert!(matches!(
            parse_one_restraint("outside box 0. 0. 0. 4. 5. 6."),
            RestraintSpec::OutsideBox { min, max }
                if min == [0.0, 0.0, 0.0] && max == [4.0, 5.0, 6.0]
        ));
    }

    #[test]
    fn parse_inside_ellipsoid() {
        assert!(matches!(
            parse_one_restraint("inside ellipsoid 0. 0. 0. 3. 2. 1.5 1."),
            RestraintSpec::InsideEllipsoid { axes, exponent, .. }
                if axes == [3.0, 2.0, 1.5] && exponent == 1.0
        ));
    }

    #[test]
    fn parse_outside_ellipsoid() {
        assert!(matches!(
            parse_one_restraint("outside ellipsoid 1. 1. 1. 2. 2. 2. 1."),
            RestraintSpec::OutsideEllipsoid { center, .. } if center == [1.0, 1.0, 1.0]
        ));
    }

    #[test]
    fn parse_inside_cylinder() {
        assert!(matches!(
            parse_one_restraint("inside cylinder 0. 0. 0. 1. 0. 0. 2. 8."),
            RestraintSpec::InsideCylinder { axis, radius, length, .. }
                if axis == [1.0, 0.0, 0.0] && radius == 2.0 && length == 8.0
        ));
    }

    #[test]
    fn parse_outside_cylinder() {
        assert!(matches!(
            parse_one_restraint("outside cylinder 0. 0. 0. 0. 0. 1. 3. 5."),
            RestraintSpec::OutsideCylinder { axis, .. } if axis == [0.0, 0.0, 1.0]
        ));
    }

    #[test]
    fn parse_ellipsoid_rejects_zero_axis() {
        let src = "output out.pdb\n\nstructure mol.pdb\n  number 1\n  \
                   inside ellipsoid 0. 0. 0. 0. 2. 2. 1.\nend structure\n";
        assert!(parse(src).is_err(), "zero semi-axis must be rejected");
    }

    #[test]
    fn parse_cylinder_rejects_zero_axis() {
        let src = "output out.pdb\n\nstructure mol.pdb\n  number 1\n  \
                   inside cylinder 0. 0. 0. 0. 0. 0. 2. 8.\nend structure\n";
        assert!(parse(src).is_err(), "zero cylinder axis must be rejected");
    }

    // ── ac-001: `profile` keyword round-trip ─────────────────────────────────

    #[test]
    fn parse_profile_gaussian_plane_defaults_to_histogram() {
        let spec = parse_one_restraint("profile gaussian plane 0. 0. 1. 0. 0. 0. mu 10. sigma 2.");
        match spec {
            RestraintSpec::Profile {
                dist,
                geometry,
                input_kind,
            } => {
                assert_eq!(
                    geometry,
                    ProfileGeometry::Plane {
                        normal: [0.0, 0.0, 1.0],
                        point: [0.0, 0.0, 0.0],
                    }
                );
                assert_eq!(
                    dist,
                    ProfileDistribution::Gaussian {
                        mu: 10.0,
                        sigma: 2.0,
                    }
                );
                // §6.2 default-safe: unflagged input is a histogram.
                assert_eq!(input_kind, ProfileInputKind::Histogram);
            }
            other => panic!("expected Profile, got {other:?}"),
        }
    }

    #[test]
    fn parse_profile_density_flag_flips_input_kind() {
        let spec =
            parse_one_restraint("profile tabulated radial 0. 0. 0. density 1. 0.9 3. 0.5 5. 0.1");
        match spec {
            RestraintSpec::Profile {
                dist,
                geometry,
                input_kind,
            } => {
                assert_eq!(geometry, ProfileGeometry::Radial { center: [0.0; 3] });
                assert_eq!(input_kind, ProfileInputKind::Density);
                match dist {
                    ProfileDistribution::Tabulated { nodes } => {
                        assert_eq!(nodes, vec![(1.0, 0.9), (3.0, 0.5), (5.0, 0.1)]);
                    }
                    other => panic!("expected Tabulated, got {other:?}"),
                }
            }
            other => panic!("expected Profile, got {other:?}"),
        }
    }

    #[test]
    fn parse_profile_histogram_flag_is_explicit_and_case_insensitive() {
        let spec = parse_one_restraint("profile tabulated radial 0. 0. 0. HISTOGRAM 1. 0.9 3. 0.5");
        assert!(matches!(
            spec,
            RestraintSpec::Profile {
                input_kind: ProfileInputKind::Histogram,
                ..
            }
        ));
    }

    #[test]
    fn parse_profile_erf_plane_with_falling_side() {
        let spec =
            parse_one_restraint("profile erf plane 1. 0. 0. 5. 0. 0. xi0 5. width 1.5 falling");
        match spec {
            RestraintSpec::Profile { dist, .. } => assert_eq!(
                dist,
                ProfileDistribution::Erf {
                    xi0: 5.0,
                    w: 1.5,
                    rising: false,
                }
            ),
            other => panic!("expected Profile, got {other:?}"),
        }
    }

    #[test]
    fn parse_profile_exponential_radial() {
        let spec = parse_one_restraint("profile exponential radial 1. 2. 3. lambda 2.5");
        match spec {
            RestraintSpec::Profile { dist, geometry, .. } => {
                assert_eq!(
                    geometry,
                    ProfileGeometry::Radial {
                        center: [1.0, 2.0, 3.0],
                    }
                );
                assert_eq!(dist, ProfileDistribution::Exponential { lambda: 2.5 });
            }
            other => panic!("expected Profile, got {other:?}"),
        }
    }

    #[test]
    fn parse_profile_inside_atoms_block_scopes_indices() {
        let src = "\
output out.pdb

structure mol.pdb
  number 1
  atoms 3 7 9
    profile gaussian plane 0. 0. 1. 0. 0. 0. mu 4. sigma 1.
  end atoms
end structure
";
        let inp = parse(src).expect("parse failed");
        let s = &inp.structures[0];
        // The molecule-wide restraint list is untouched.
        assert!(s.mol_restraints.is_empty());
        assert_eq!(s.atom_groups.len(), 1);
        let group = &s.atom_groups[0];
        assert_eq!(group.atom_indices, vec![3, 7, 9]);
        assert_eq!(group.restraints.len(), 1);
        assert!(matches!(
            group.restraints[0],
            RestraintSpec::Profile {
                geometry: ProfileGeometry::Plane { .. },
                ..
            }
        ));
    }

    #[test]
    fn parse_profile_rejects_bad_token_with_lineno() {
        // sigma is not a number → Parse error carrying the offending line.
        let src = "\
output out.pdb

structure mol.pdb
  number 1
  profile gaussian plane 0. 0. 1. 0. 0. 0. mu 10. sigma oops
end structure
";
        let err = parse(src).expect_err("should fail");
        match err {
            ScriptError::Parse { line, .. } => assert_eq!(line, 5),
            other => panic!("expected Parse error, got {other:?}"),
        }
    }

    #[test]
    fn parse_profile_rejects_unknown_distribution() {
        let src = "\
output out.pdb

structure mol.pdb
  number 1
  profile wibble plane 0. 0. 1. 0. 0. 0.
end structure
";
        assert!(parse(src).is_err(), "unknown distribution must be rejected");
    }

    #[test]
    fn parse_profile_rejects_radial_histogram_node_at_origin() {
        // Radial + histogram + xi=0 is the shell singularity the spline layer
        // refuses; the parser must reject it up front.
        let src = "\
output out.pdb

structure mol.pdb
  number 1
  profile tabulated radial 0. 0. 0. 0. 1. 3. 0.5
end structure
";
        assert!(
            parse(src).is_err(),
            "radial histogram node at xi=0 must be rejected"
        );
    }

    #[test]
    fn parse_profile_rejects_non_positive_sigma() {
        let src = "\
output out.pdb

structure mol.pdb
  number 1
  profile gaussian plane 0. 0. 1. 0. 0. 0. mu 4. sigma 0.
end structure
";
        assert!(parse(src).is_err(), "non-positive sigma must be rejected");
    }
}
