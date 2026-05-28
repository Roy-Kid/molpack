//! Sub-parsers for the restraint keywords (`inside`, `outside`,
//! `over`, `below`) that may appear inside a `structure` block or an
//! `atoms … end atoms` sub-block.
//!
//! This is the single source of truth mapping `.inp` keyword → shape →
//! [`RestraintSpec`]. Every one of the 14 concrete `*Restraint` kinds in
//! [`crate::restraint`] is reachable from here; [`super::build`] lowers each
//! `RestraintSpec` to its matching Rust constructor.

use super::error::ScriptError;
use super::parser::{RestraintSpec, parse_err, parse_f64, parse_vec3};

/// Reject a non-positive size/radius with a consistent message.
fn require_positive(value: f64, what: &str, lineno: usize) -> Result<f64, ScriptError> {
    if value <= 0.0 {
        return Err(parse_err(lineno, format!("{what} must be > 0")));
    }
    Ok(value)
}

/// Parse `inside <shape> …` — `box`, `cube`, `sphere`, `ellipsoid`, `cylinder`.
pub(super) fn parse_inside(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let shape = tokens
        .get(1)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, "`inside` requires a shape keyword"))?;
    match shape.as_str() {
        "box" => {
            let (min, max) = parse_box(tokens, lineno)?;
            Ok(RestraintSpec::InsideBox { min, max })
        }
        "cube" => {
            let (origin, side) = parse_cube(tokens, lineno)?;
            Ok(RestraintSpec::InsideCube { origin, side })
        }
        "sphere" => {
            let (center, radius) = parse_sphere(tokens, "inside sphere", lineno)?;
            Ok(RestraintSpec::InsideSphere { center, radius })
        }
        "ellipsoid" => {
            let (center, axes, exponent) = parse_ellipsoid(tokens, lineno)?;
            Ok(RestraintSpec::InsideEllipsoid {
                center,
                axes,
                exponent,
            })
        }
        "cylinder" => {
            let (center, axis, radius, length) = parse_cylinder(tokens, lineno)?;
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

/// Parse `outside <shape> …` — `box`, `cube`, `sphere`, `ellipsoid`, `cylinder`.
pub(super) fn parse_outside(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let shape = tokens
        .get(1)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, "`outside` requires a shape keyword"))?;
    match shape.as_str() {
        "box" => {
            let (min, max) = parse_box(tokens, lineno)?;
            Ok(RestraintSpec::OutsideBox { min, max })
        }
        "cube" => {
            let (origin, side) = parse_cube(tokens, lineno)?;
            Ok(RestraintSpec::OutsideCube { origin, side })
        }
        "sphere" => {
            let (center, radius) = parse_sphere(tokens, "outside sphere", lineno)?;
            Ok(RestraintSpec::OutsideSphere { center, radius })
        }
        "ellipsoid" => {
            let (center, axes, exponent) = parse_ellipsoid(tokens, lineno)?;
            Ok(RestraintSpec::OutsideEllipsoid {
                center,
                axes,
                exponent,
            })
        }
        "cylinder" => {
            let (center, axis, radius, length) = parse_cylinder(tokens, lineno)?;
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

/// Parse `over plane …` or `over gaussian …` (atom must lie above the surface).
pub(super) fn parse_above(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    match surface_kind(tokens, "over", lineno)?.as_str() {
        "plane" => {
            let (normal, distance) = parse_plane(tokens, "over plane", lineno)?;
            Ok(RestraintSpec::AbovePlane { normal, distance })
        }
        "gaussian" => {
            let (cx, cy, sx, sy, z0, height) = parse_gaussian(tokens, "over gaussian", lineno)?;
            Ok(RestraintSpec::AboveGaussian {
                cx,
                cy,
                sx,
                sy,
                z0,
                height,
            })
        }
        other => Err(parse_err(
            lineno,
            format!("unsupported `over` surface `{other}`"),
        )),
    }
}

/// Parse `below plane …` or `below gaussian …` (atom must lie below the surface).
pub(super) fn parse_below(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    match surface_kind(tokens, "below", lineno)?.as_str() {
        "plane" => {
            let (normal, distance) = parse_plane(tokens, "below plane", lineno)?;
            Ok(RestraintSpec::BelowPlane { normal, distance })
        }
        "gaussian" => {
            let (cx, cy, sx, sy, z0, height) = parse_gaussian(tokens, "below gaussian", lineno)?;
            Ok(RestraintSpec::BelowGaussian {
                cx,
                cy,
                sx,
                sy,
                z0,
                height,
            })
        }
        other => Err(parse_err(lineno, format!("unsupported surface `{other}`"))),
    }
}

// ── shared shape sub-parsers ────────────────────────────────────────────────

fn surface_kind(tokens: &[&str], kw: &str, lineno: usize) -> Result<String, ScriptError> {
    tokens
        .get(1)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, format!("`{kw}` requires a surface keyword")))
}

fn parse_box(tokens: &[&str], lineno: usize) -> Result<([f64; 3], [f64; 3]), ScriptError> {
    let min = parse_vec3(tokens, 2, "box min", lineno)?;
    let max = parse_vec3(tokens, 5, "box max", lineno)?;
    Ok((min, max))
}

fn parse_cube(tokens: &[&str], lineno: usize) -> Result<([f64; 3], f64), ScriptError> {
    let origin = parse_vec3(tokens, 2, "cube origin", lineno)?;
    let side = require_positive(
        parse_f64(tokens, 5, "cube side", lineno)?,
        "cube side",
        lineno,
    )?;
    Ok((origin, side))
}

fn parse_sphere(tokens: &[&str], ctx: &str, lineno: usize) -> Result<([f64; 3], f64), ScriptError> {
    let center = parse_vec3(tokens, 2, ctx, lineno)?;
    let radius = require_positive(
        parse_f64(tokens, 5, &format!("{ctx} radius"), lineno)?,
        &format!("{ctx} radius"),
        lineno,
    )?;
    Ok((center, radius))
}

fn parse_ellipsoid(
    tokens: &[&str],
    lineno: usize,
) -> Result<([f64; 3], [f64; 3], f64), ScriptError> {
    let center = parse_vec3(tokens, 2, "ellipsoid center", lineno)?;
    let axes = parse_vec3(tokens, 5, "ellipsoid axes", lineno)?;
    for a in axes {
        require_positive(a, "ellipsoid axis", lineno)?;
    }
    let exponent = require_positive(
        parse_f64(tokens, 8, "ellipsoid exponent", lineno)?,
        "ellipsoid exponent",
        lineno,
    )?;
    Ok((center, axes, exponent))
}

fn parse_cylinder(
    tokens: &[&str],
    lineno: usize,
) -> Result<([f64; 3], [f64; 3], f64, f64), ScriptError> {
    let center = parse_vec3(tokens, 2, "cylinder center", lineno)?;
    let axis = parse_vec3(tokens, 5, "cylinder axis", lineno)?;
    let radius = require_positive(
        parse_f64(tokens, 8, "cylinder radius", lineno)?,
        "cylinder radius",
        lineno,
    )?;
    let length = require_positive(
        parse_f64(tokens, 9, "cylinder length", lineno)?,
        "cylinder length",
        lineno,
    )?;
    Ok((center, axis, radius, length))
}

fn parse_plane(tokens: &[&str], ctx: &str, lineno: usize) -> Result<([f64; 3], f64), ScriptError> {
    let normal = parse_vec3(tokens, 2, &format!("{ctx} normal"), lineno)?;
    let distance = parse_f64(tokens, 5, &format!("{ctx} distance"), lineno)?;
    Ok((normal, distance))
}

#[allow(clippy::type_complexity)]
fn parse_gaussian(
    tokens: &[&str],
    ctx: &str,
    lineno: usize,
) -> Result<(f64, f64, f64, f64, f64, f64), ScriptError> {
    let cx = parse_f64(tokens, 2, ctx, lineno)?;
    let cy = parse_f64(tokens, 3, ctx, lineno)?;
    let sx = require_positive(
        parse_f64(tokens, 4, ctx, lineno)?,
        &format!("{ctx} sx"),
        lineno,
    )?;
    let sy = require_positive(
        parse_f64(tokens, 5, ctx, lineno)?,
        &format!("{ctx} sy"),
        lineno,
    )?;
    let z0 = parse_f64(tokens, 6, ctx, lineno)?;
    let height = parse_f64(tokens, 7, ctx, lineno)?;
    Ok((cx, cy, sx, sy, z0, height))
}
