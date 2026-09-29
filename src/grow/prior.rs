//! Geometric conformer priors for the growth solver.
//!
//! A prior is **data the caller supplies** — torsion state weights, spreads,
//! or (with the Task 4 calibration helpers) a target characteristic ratio.
//! Never a force-field object: molpack is purely geometric, and a caller who
//! wants force-field-derived weights computes them outside and passes the
//! numbers in.

use molrs::types::F;
use rand::Rng;

use crate::grow::internal::wrap_pi;
use crate::random::uniform01;

const PI: F = std::f64::consts::PI as F;
const TWO_PI: F = std::f64::consts::TAU as F;

/// Prior over a free torsion variable, sampled once per growth trial.
///
/// The choice is load-bearing, not cosmetic: uniform sampling produces
/// freely-rotating-chain statistics (C∞ = 2.00 for a tetrahedral backbone),
/// which undershoots melt chain dimensions by ~40% for PEO — and crowding at
/// melt density does not repair it. See the spec's Domain basis §5.1.
#[derive(Debug, Clone)]
pub enum TorsionPrior {
    /// Uniform on `(-π, π]`. Negative-control / regression-baseline use only.
    Uniform,
    /// von Mises-like spread of concentration `kappa` around the template's
    /// own torsion value. `kappa = 0` degenerates to [`Uniform`][Self::Uniform].
    Template { kappa: F },
    /// RIS-style discrete states as `(angle_rad, weight)` pairs. Weights are
    /// normalized at use; they need not sum to 1.
    States(Vec<(F, F)>),
}

impl TorsionPrior {
    /// Draw one torsion value.
    ///
    /// `template_value` is the template's own dihedral for this variable —
    /// consumed only by [`Template`][Self::Template]; the other priors are
    /// absolute. All angles share [`dihedral`](super::internal)'s convention:
    /// values in `(-π, π]`, trans ≈ ±π for a backbone.
    pub fn sample(&self, template_value: F, rng: &mut impl Rng) -> F {
        match self {
            TorsionPrior::Uniform => -PI + TWO_PI * uniform01(rng),
            TorsionPrior::Template { kappa } => {
                if *kappa <= 0.0 {
                    return -PI + TWO_PI * uniform01(rng);
                }
                // Wrapped-normal approximation of the von Mises distribution
                // (σ = κ^{-1/2}): indistinguishable for the concentrations a
                // conformer prior uses, and cheap. Box–Muller from two
                // uniforms keeps the RNG stream consumption fixed per draw.
                let u1 = uniform01(rng).max(F::MIN_POSITIVE);
                let u2 = uniform01(rng);
                let gauss = (-2.0 * u1.ln()).sqrt() * (TWO_PI * u2).cos();
                wrap_pi(template_value + gauss / kappa.sqrt())
            }
            TorsionPrior::States(states) => {
                let total: F = states.iter().map(|&(_, w)| w.max(0.0)).sum();
                if total <= 0.0 || states.is_empty() {
                    return -PI + TWO_PI * uniform01(rng);
                }
                let mut ticket = uniform01(rng) * total;
                for &(angle, w) in states {
                    let w = w.max(0.0);
                    if ticket < w {
                        return angle;
                    }
                    ticket -= w;
                }
                states[states.len() - 1].0
            }
        }
    }

    /// Three-state trans/gauche± prior whose trans fraction is solved from a
    /// target characteristic ratio — the single-scalar calibration of the
    /// spec's Domain basis §5.1.
    ///
    /// With valence angle `theta` (radians), the freely-rotating baseline is
    /// `C_FRC = (1 − cos θ)/(1 + cos θ)`, and independent hindered rotation
    /// gives `C∞ = C_FRC · (1 + ⟨cos φ′⟩)/(1 − ⟨cos φ′⟩)` with `φ′` measured
    /// from trans (Flory). Solving for the trans fraction of states at
    /// trans = π and gauche± = ±π/3 (absolute convention, matching
    /// [`sample`][Self::sample]):
    /// `x = (r − 1)/(r + 1)`, `r = C∞/C_FRC`, `p_t = (2x + 1)/3`.
    ///
    /// PEO (C∞ = 5.5, tetrahedral θ): `p_t ≈ 0.645`.
    pub fn three_state_from_c_inf(c_inf: F, theta: F) -> Self {
        let cos_t = theta.cos();
        let c_frc = (1.0 - cos_t) / (1.0 + cos_t);
        let r = c_inf / c_frc;
        let x = (r - 1.0) / (r + 1.0);
        let p_t = ((2.0 * x + 1.0) / 3.0).clamp(0.0, 1.0);
        let p_g = (1.0 - p_t) / 2.0;
        TorsionPrior::States(vec![(PI, p_t), (PI / 3.0, p_g), (-PI / 3.0, p_g)])
    }
}

/// Prior over the placement (bond) angles of grown sites.
///
/// All-atom templates keep their angles verbatim ([`Template`][Self::Template],
/// the default): bond angles are stiff coordinates there. Coarse-grained
/// templates are different — CG angle potentials are fitted to reproduce a
/// persistence length, so the angle is a soft, statistically distributed
/// coordinate, and for a torsion-free CG chain the angle prior is the *only*
/// control over chain stiffness (spec Domain basis §5.5).
#[derive(Debug, Clone)]
pub enum AnglePrior {
    /// Copy every site's angle from the template (the all-atom behavior).
    Template,
    /// Discrete worm-like chain: bond-deflection angles θ′ sampled from the
    /// tilted density `p(cos θ′) ∝ exp(κ cos θ′)`.
    Wlc { kappa: F },
}

impl AnglePrior {
    /// Calibrate the WLC tilt from a target characteristic ratio:
    /// `⟨cos θ′⟩ = (C∞ − 1)/(C∞ + 1)`, then κ from the inverse Langevin
    /// function (Cohen's approximation `κ ≈ x(3 − x²)/(1 − x²)`).
    ///
    /// The relation assumes **uniform torsions** — combining a Wlc angle
    /// prior with a non-uniform [`TorsionPrior`] double-counts stiffness
    /// and needs an empirical recalibration.
    ///
    /// Kremer–Grest melts: `wlc_from_c_inf(1.76)` (spec §5.5).
    pub fn wlc_from_c_inf(c_inf: F) -> Self {
        let x = ((c_inf - 1.0) / (c_inf + 1.0)).clamp(0.0, 0.999_999);
        let kappa = x * (3.0 - x * x) / (1.0 - x * x);
        AnglePrior::Wlc { kappa }
    }

    /// Draw one *interior* placement angle (`π − θ′`), which is what NeRF
    /// consumes. [`Template`][Self::Template] returns the template's value
    /// and consumes no randomness.
    pub(crate) fn sample_interior(&self, template_angle: F, rng: &mut impl Rng) -> F {
        match self {
            AnglePrior::Template => template_angle,
            AnglePrior::Wlc { kappa } => {
                let k = kappa.max(1e-9);
                let u = uniform01(rng);
                // Inverse CDF of the tilted density on cos θ′ ∈ [−1, 1].
                let c = (1.0 + (u + (1.0 - u) * (-2.0 * k).exp()).ln() / k).clamp(-1.0, 1.0);
                PI - c.acos()
            }
        }
    }
}
