//! Parsing for the `profile` keyword — the script-layer surface of the
//! profile-distribution restraint family.
//!
//! Grammar (Packmol-style: lowercase keyword, whitespace-delimited tokens):
//!
//! ```text
//! profile gaussian     <geometry> mu <μ> sigma <σ>
//! profile erf          <geometry> xi0 <ξ0> width <w> [rising|falling]
//! profile tanh         <geometry> xi0 <ξ0> width <w> [rising|falling]
//! profile exponential  <geometry> lambda <λ>
//! profile tabulated    <geometry> [density|histogram] <ξ1> <ρ1> <ξ2> <ρ2> …
//! ```
//!
//! where `<geometry>` is one of:
//!
//! ```text
//! plane    <nx ny nz> <x0 y0 z0>      (signed distance to a plane)
//! radial   <cx cy cz>                 (distance from a centre)
//! cylinder <cx cy cz> <dx dy dz> <L>  (perpendicular distance to an axis)
//! ```
//!
//! The trailing `density|histogram` flag is OPTIONAL and case-insensitive; the
//! default — matching the domain-safe §6.2 rule — is `histogram`, so a count
//! histogram has its shell-volume Jacobian divided out before Boltzmann
//! inversion. The flag is meaningful only for radial/cylindrical geometries
//! (planar has a unit Jacobian) and only the `tabulated` distribution reads
//! inline nodes after it; for the analytic shapes it simply records intent.
//!
//! This module owns the profile-specific spec sub-types and the line parser; it
//! is re-exported from [`super::parser`] so the public AST surface
//! ([`RestraintSpec::Profile`](super::parser::RestraintSpec)) stays in one
//! place. Keeping it separate holds `parser.rs` under its 800-LOC budget.

use super::error::ScriptError;
use super::parser::{RestraintSpec, parse_err, parse_f64, parse_vec3};

/// Whether a tabulated profile's values are a volumetric density or a raw count
/// histogram. Mirrors the restraint layer's `InputKind` without depending on it,
/// keeping the parser AST free of restraint-side types. Default is
/// [`Histogram`](ProfileInputKind::Histogram) — the §6.2 domain-safe choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ProfileInputKind {
    /// Already a volumetric density ρ*(ξ); invert directly.
    Density,
    /// A count histogram n(ξ); divide by the shell Jacobian before inverting.
    #[default]
    Histogram,
}

/// Reaction-coordinate geometry parsed from a `profile` line. The numeric
/// payload mirrors the [`Coordinate`](crate::restraint::profile::Coordinate)
/// constructors the lowering step feeds it into.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProfileGeometry {
    /// `plane <nx ny nz> <x0 y0 z0>` — signed distance to a plane.
    Plane {
        /// Plane normal (need not be unit; normalized at lowering).
        normal: [f64; 3],
        /// A point on the plane.
        point: [f64; 3],
    },
    /// `radial <cx cy cz>` — distance from a centre.
    Radial {
        /// Sphere centre.
        center: [f64; 3],
    },
    /// `cylinder <cx cy cz> <dx dy dz> <L>` — perpendicular distance to an axis.
    Cylinder {
        /// A point on the axis.
        origin: [f64; 3],
        /// Axis direction (need not be unit; normalized at lowering).
        axis: [f64; 3],
        /// Cylinder length L (Å), used by the shell Jacobian.
        length: f64,
    },
}

/// Target-distribution shape parsed from a `profile` line, with its own numeric
/// parameters. Mirrors the analytic
/// [`Distribution`](crate::restraint::profile::Distribution) plus a tabulated
/// node list.
#[derive(Debug, Clone, PartialEq)]
pub enum ProfileDistribution {
    /// `gaussian … mu <μ> sigma <σ>` — harmonic well at μ, width σ.
    Gaussian { mu: f64, sigma: f64 },
    /// `erf … xi0 <ξ0> width <w> [rising|falling]` — error-function step.
    Erf { xi0: f64, w: f64, rising: bool },
    /// `tanh … xi0 <ξ0> width <w> [rising|falling]` — hyperbolic-tangent step.
    Tanh { xi0: f64, w: f64, rising: bool },
    /// `exponential … lambda <λ>` — exponential decay, constant inward force.
    Exponential { lambda: f64 },
    /// `tabulated … [density|histogram] <ξ1> <ρ1> …` — inline node table.
    Tabulated { nodes: Vec<(f64, f64)> },
}

/// Parse a whole `profile …` line into a [`RestraintSpec::Profile`].
///
/// `tokens[0]` is the `profile` keyword; the remainder carry the distribution
/// name, geometry, and parameters. Token counts and keyword spelling are
/// validated through the shared [`parse_f64`] / [`parse_vec3`] / [`parse_err`]
/// helpers so error messages match the rest of the grammar (line number + a
/// human-readable cause).
pub fn parse_profile(tokens: &[&str], lineno: usize) -> Result<RestraintSpec, ScriptError> {
    let dist_kw = tokens
        .get(1)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, "`profile` requires a distribution keyword"))?;

    // The geometry block starts at token 2 and consumes a known number of
    // tokens; parameters follow. `parse_geometry` returns the index of the
    // first post-geometry token.
    let (geometry, after_geom) = parse_geometry(tokens, 2, lineno)?;

    let (dist, input_kind) = match dist_kw.as_str() {
        "gaussian" => {
            analytic_with_flag(parse_gaussian(tokens, after_geom, lineno)?, tokens, lineno)?
        }
        "erf" => analytic_with_flag(
            parse_step(tokens, after_geom, lineno, true)?,
            tokens,
            lineno,
        )?,
        "tanh" => analytic_with_flag(
            parse_step(tokens, after_geom, lineno, false)?,
            tokens,
            lineno,
        )?,
        "exponential" => analytic_with_flag(
            parse_exponential(tokens, after_geom, lineno)?,
            tokens,
            lineno,
        )?,
        "tabulated" => parse_tabulated(tokens, after_geom, lineno)?,
        other => {
            return Err(parse_err(
                lineno,
                format!("unknown `profile` distribution `{other}`"),
            ));
        }
    };

    // Fail fast at the boundary: a malformed tabulated grid must be rejected here
    // so the (public, infallible) lowering step never has to.
    if let ProfileDistribution::Tabulated { nodes } = &dist {
        validate_table(nodes, &geometry, input_kind, lineno)?;
    }

    Ok(RestraintSpec::Profile {
        dist,
        geometry,
        input_kind,
    })
}

/// Read the optional trailing `density|histogram` flag for an analytic
/// distribution, given the parser's `(distribution, end_index)`. The flag is
/// case-insensitive and defaults to histogram (§6.2). Any token at `end` that is
/// neither flag is a hard error — a stray trailing token must not be ignored.
fn analytic_with_flag(
    parsed: (ProfileDistribution, usize),
    tokens: &[&str],
    lineno: usize,
) -> Result<(ProfileDistribution, ProfileInputKind), ScriptError> {
    let (dist, end) = parsed;
    let input_kind = match tokens.get(end).map(|s| s.to_ascii_lowercase()).as_deref() {
        None => ProfileInputKind::default(),
        Some("density") => ProfileInputKind::Density,
        Some("histogram") => ProfileInputKind::Histogram,
        Some(other) => {
            return Err(parse_err(
                lineno,
                format!("unexpected trailing token `{other}`; expected `density` or `histogram`"),
            ));
        }
    };
    Ok((dist, input_kind))
}

/// Reject a tabulated grid the spline layer would later refuse: < 2 nodes,
/// non-finite values, negative density, a non-strictly-increasing ξ grid, or —
/// for a radial/cylindrical histogram — a node at ξ ≤ 0 where the shell factor
/// `dV/dξ` vanishes. Validating here keeps [`super::build`] lowering total.
fn validate_table(
    nodes: &[(f64, f64)],
    geometry: &ProfileGeometry,
    input_kind: ProfileInputKind,
    lineno: usize,
) -> Result<(), ScriptError> {
    if nodes.len() < 2 {
        return Err(parse_err(
            lineno,
            format!(
                "`profile tabulated` needs >= 2 node pairs, got {}",
                nodes.len()
            ),
        ));
    }
    let shell_geometry = !matches!(geometry, ProfileGeometry::Plane { .. });
    let is_histogram = input_kind == ProfileInputKind::Histogram;
    for (i, &(xi, rho)) in nodes.iter().enumerate() {
        if !xi.is_finite() || !rho.is_finite() {
            return Err(parse_err(
                lineno,
                format!("`profile tabulated` node {i} is non-finite"),
            ));
        }
        if rho < 0.0 {
            return Err(parse_err(
                lineno,
                format!("`profile tabulated` density at node {i} is negative ({rho})"),
            ));
        }
        if i > 0 && xi <= nodes[i - 1].0 {
            return Err(parse_err(
                lineno,
                format!("`profile tabulated` xi grid is not strictly increasing at node {i}"),
            ));
        }
        if shell_geometry && is_histogram && xi <= 0.0 {
            return Err(parse_err(
                lineno,
                format!(
                    "`profile tabulated` radial/cylindrical histogram node {i} has xi <= 0 \
                     where the shell factor vanishes; use density input or xi > 0"
                ),
            ));
        }
    }
    Ok(())
}

/// Parse the geometry block starting at `start`, returning the parsed geometry
/// and the index of the first token after it.
fn parse_geometry(
    tokens: &[&str],
    start: usize,
    lineno: usize,
) -> Result<(ProfileGeometry, usize), ScriptError> {
    let geom_kw = tokens
        .get(start)
        .map(|s| s.to_ascii_lowercase())
        .ok_or_else(|| parse_err(lineno, "`profile` requires a geometry keyword"))?;
    match geom_kw.as_str() {
        "plane" => {
            let normal = parse_vec3(tokens, start + 1, "profile plane normal", lineno)?;
            let point = parse_vec3(tokens, start + 4, "profile plane point", lineno)?;
            if normal == [0.0; 3] {
                return Err(parse_err(
                    lineno,
                    "profile plane normal must be a non-zero direction",
                ));
            }
            Ok((ProfileGeometry::Plane { normal, point }, start + 7))
        }
        "radial" => {
            let center = parse_vec3(tokens, start + 1, "profile radial center", lineno)?;
            Ok((ProfileGeometry::Radial { center }, start + 4))
        }
        "cylinder" => {
            let origin = parse_vec3(tokens, start + 1, "profile cylinder origin", lineno)?;
            let axis = parse_vec3(tokens, start + 4, "profile cylinder axis", lineno)?;
            if axis == [0.0; 3] {
                return Err(parse_err(
                    lineno,
                    "profile cylinder axis must be a non-zero direction",
                ));
            }
            let length = parse_f64(tokens, start + 7, "profile cylinder length", lineno)?;
            Ok((
                ProfileGeometry::Cylinder {
                    origin,
                    axis,
                    length,
                },
                start + 8,
            ))
        }
        other => Err(parse_err(
            lineno,
            format!("unknown `profile` geometry `{other}`"),
        )),
    }
}

/// `… mu <μ> sigma <σ>` — keyword-tagged so order is checked, not positional.
/// Returns the distribution and the index of the first token after `σ`.
fn parse_gaussian(
    tokens: &[&str],
    start: usize,
    lineno: usize,
) -> Result<(ProfileDistribution, usize), ScriptError> {
    expect_keyword(tokens, start, "mu", lineno)?;
    let mu = parse_f64(tokens, start + 1, "profile gaussian mu", lineno)?;
    expect_keyword(tokens, start + 2, "sigma", lineno)?;
    let sigma = parse_f64(tokens, start + 3, "profile gaussian sigma", lineno)?;
    if sigma <= 0.0 {
        return Err(parse_err(
            lineno,
            format!("profile gaussian sigma must be positive, got {sigma}"),
        ));
    }
    Ok((ProfileDistribution::Gaussian { mu, sigma }, start + 4))
}

/// `… xi0 <ξ0> width <w> [rising|falling]` shared by erf and tanh. Returns the
/// distribution and the index of the first token after the (optional) side flag.
fn parse_step(
    tokens: &[&str],
    start: usize,
    lineno: usize,
    is_erf: bool,
) -> Result<(ProfileDistribution, usize), ScriptError> {
    expect_keyword(tokens, start, "xi0", lineno)?;
    let xi0 = parse_f64(tokens, start + 1, "profile step xi0", lineno)?;
    expect_keyword(tokens, start + 2, "width", lineno)?;
    let w = parse_f64(tokens, start + 3, "profile step width", lineno)?;
    if w <= 0.0 {
        return Err(parse_err(
            lineno,
            format!("profile step width must be positive, got {w}"),
        ));
    }
    // Optional rising|falling flag (default rising). A trailing density|histogram
    // flag is NOT a side keyword, so consume the side flag only when present.
    let (rising, end) = match tokens
        .get(start + 4)
        .map(|s| s.to_ascii_lowercase())
        .as_deref()
    {
        Some("rising") => (true, start + 5),
        Some("falling") => (false, start + 5),
        _ => (true, start + 4),
    };
    let dist = if is_erf {
        ProfileDistribution::Erf { xi0, w, rising }
    } else {
        ProfileDistribution::Tanh { xi0, w, rising }
    };
    Ok((dist, end))
}

/// `… lambda <λ>`. Returns the distribution and the index after `λ`.
fn parse_exponential(
    tokens: &[&str],
    start: usize,
    lineno: usize,
) -> Result<(ProfileDistribution, usize), ScriptError> {
    expect_keyword(tokens, start, "lambda", lineno)?;
    let lambda = parse_f64(tokens, start + 1, "profile exponential lambda", lineno)?;
    if lambda <= 0.0 {
        return Err(parse_err(
            lineno,
            format!("profile exponential lambda must be positive, got {lambda}"),
        ));
    }
    Ok((ProfileDistribution::Exponential { lambda }, start + 2))
}

/// `… [density|histogram] <ξ1> <ρ1> <ξ2> <ρ2> …`. The optional flag selects the
/// input kind; the remaining tokens are read as (ξ, ρ) pairs.
fn parse_tabulated(
    tokens: &[&str],
    start: usize,
    lineno: usize,
) -> Result<(ProfileDistribution, ProfileInputKind), ScriptError> {
    // Optional leading density|histogram flag.
    let (input_kind, nodes_start) =
        match tokens.get(start).map(|s| s.to_ascii_lowercase()).as_deref() {
            Some("density") => (ProfileInputKind::Density, start + 1),
            Some("histogram") => (ProfileInputKind::Histogram, start + 1),
            _ => (ProfileInputKind::default(), start),
        };

    let rest = &tokens[nodes_start.min(tokens.len())..];
    if rest.is_empty() {
        return Err(parse_err(
            lineno,
            "`profile tabulated` requires at least one (xi rho) node pair",
        ));
    }
    if rest.len() % 2 != 0 {
        return Err(parse_err(
            lineno,
            format!(
                "`profile tabulated` needs an even count of node values (xi rho …), got {}",
                rest.len()
            ),
        ));
    }
    let mut nodes = Vec::with_capacity(rest.len() / 2);
    for pair in rest.chunks_exact(2) {
        let xi = pair[0]
            .parse::<f64>()
            .map_err(|_| parse_err(lineno, format!("invalid tabulated xi `{}`", pair[0])))?;
        let rho = pair[1]
            .parse::<f64>()
            .map_err(|_| parse_err(lineno, format!("invalid tabulated rho `{}`", pair[1])))?;
        nodes.push((xi, rho));
    }
    Ok((ProfileDistribution::Tabulated { nodes }, input_kind))
}

/// Assert the token at `idx` is exactly `kw` (case-insensitive), else error.
fn expect_keyword(tokens: &[&str], idx: usize, kw: &str, lineno: usize) -> Result<(), ScriptError> {
    match tokens.get(idx).map(|s| s.to_ascii_lowercase()) {
        Some(t) if t == kw => Ok(()),
        Some(t) => Err(parse_err(
            lineno,
            format!("expected `{kw}` keyword, got `{t}`"),
        )),
        None => Err(parse_err(lineno, format!("missing `{kw}` keyword"))),
    }
}
