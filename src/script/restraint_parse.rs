//! Sub-parsers for the restraint keywords (`inside`, `outside`,
//! `over plane`, `below plane`) that may appear inside a `structure`
//! block or an `atoms … end atoms` sub-block.
//!
//! Split out from [`super::parser`] so the restraint grammar lives in one
//! place — the natural home for adding the kinds that are defined in the
//! Rust API but not yet reachable from `.inp` (cube / ellipsoid / cylinder
//! / gaussian and the `outside-*` variants).

use super::error::ScriptError;
use super::parser::{RestraintSpec, parse_err, parse_f64, parse_vec3};

/// Parse `inside <shape> …` — currently `box` and `sphere`.
pub(super) fn parse_inside(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
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
            if radius <= 0.0 {
                return Err(parse_err(lineno, "inside sphere radius must be > 0"));
            }
            Ok(RestraintSpec::InsideSphere { center, radius })
        }
        _ => Err(parse_err(
            lineno,
            format!("unsupported `inside` shape `{shape}`"),
        )),
    }
}

/// Parse `outside <shape> …` — currently `sphere`.
pub(super) fn parse_outside(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let shape = tokens
        .get(1)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, "`outside` requires a shape keyword"))?;
    match shape.as_str() {
        "sphere" => {
            let center = parse_vec3(tokens, 2, "outside sphere center", lineno)?;
            let radius = parse_f64(tokens, 5, "outside sphere radius", lineno)?;
            if radius <= 0.0 {
                return Err(parse_err(lineno, "outside sphere radius must be > 0"));
            }
            Ok(RestraintSpec::OutsideSphere { center, radius })
        }
        _ => Err(parse_err(
            lineno,
            format!("unsupported `outside` shape `{shape}`"),
        )),
    }
}

/// Parse `over plane <nx> <ny> <nz> <d>` — atom must lie above the plane.
pub(super) fn parse_plane_above(
    tokens: &[&str],
    lineno: usize,
) -> Result<RestraintSpec, ScriptError> {
    let normal = parse_vec3(tokens, 2, "over plane normal", lineno)?;
    let distance = parse_f64(tokens, 5, "over plane distance", lineno)?;
    Ok(RestraintSpec::AbovePlane { normal, distance })
}

/// Parse `below plane <nx> <ny> <nz> <d>` — atom must lie below the plane.
pub(super) fn parse_plane_below(
    tokens: &[&str],
    lineno: usize,
) -> Result<RestraintSpec, ScriptError> {
    let normal = parse_vec3(tokens, 2, "below plane normal", lineno)?;
    let distance = parse_f64(tokens, 5, "below plane distance", lineno)?;
    Ok(RestraintSpec::BelowPlane { normal, distance })
}
