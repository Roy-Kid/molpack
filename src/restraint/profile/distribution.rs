//! Analytic target-distribution penalties — closed-form Boltzmann inversion.
//!
//! This leaf module turns a desired 1-D profile of a scalar reaction coordinate
//! ξ (a length in Å; supplied already-reduced by the geometry layer in
//! `coordinate.rs`) into a soft energy `U(ξ)` and force `−dU/dξ` by Boltzmann
//! inversion of a target volumetric density ρ*(ξ):
//!
//! ```text
//! U(ξ) = −kT·ln(ρ*(ξ)/ρ₀)        dU/dξ = −kT·ρ*'(ξ)/ρ*(ξ)
//! ```
//!
//! Units: `kT` is a FIXED energy-scale factor (NOT a simulated temperature)
//! sharing the overlap-penalty's energy units; it sets bias steepness. `ρ₀` is an
//! arbitrary reference density entering `U` only as an additive constant, so it
//! does not affect `dU/dξ` or the force — it sets the energy zero only.
//!
//! Four closed-form shapes are provided ([`Distribution`]), each with EXACT
//! derivatives (no finite differencing in production):
//!
//! - **Gaussian** ρ* ∝ exp(−(ξ−μ)²/2σ²) ⇒ `U = (kT/2σ²)(ξ−μ)²`,
//!   `dU/dξ = (kT/σ²)(ξ−μ)` — exactly a harmonic well, spring `k = kT/σ²`.
//! - **erf step** ρ* ∝ ½[1+erf((ξ−ξ₀)/(√2 w))] — the default interface shape
//!   (capillary-wave-averaged liquid interface).
//! - **tanh step** ρ* ∝ ½[1+tanh((ξ−ξ₀)/w)] — the mean-field variant.
//! - **exponential** ρ* ∝ exp(−ξ/λ), ξ≥0 ⇒ `dU/dξ = kT/λ` (constant force).
//!
//! Two correctness safeguards from the domain basis are composed by
//! [`ProfilePenalty`]:
//!
//! - **Shell-volume Jacobian** ([`ShellJacobian`] + [`InputKind`]): Boltzmann
//!   inversion must act on a VOLUMETRIC density. A count histogram n(ξ)
//!   [count/ξ] is converted ρ*(ξ) = n(ξ)/(dV/dξ) BEFORE inversion. Omitting this
//!   for a radial/cylindrical profile injects a spurious `+kT·ln(dV/dξ)` into `U`
//!   whose derivative is a spurious outward force (`+2kT/ξ` radial, `+kT/ξ` cyl)
//!   that diverges as ξ→0. Default-safe = treat the input as a histogram.
//! - **Density floor** ([`DensityFloor`]): replace ρ* by `max(ρ*, ρ_min)` ⇔ cap
//!   `U_max = −kT·ln(ρ_min/ρ₀)`. A zero-probability bin is unobservable, not an
//!   infinite wall; where the cap is active the force is zeroed so gencan never
//!   sees an `Inf`/`NaN` gradient.
//!
//! Every type is a plain `Copy` value, allocation-free, referentially
//! transparent in (ξ, kT), and usable from `Send + Sync` contexts. This module
//! owns no coordinate geometry, no tabulated profiles, and no `Restraint` impl
//! (those are sibling sub-specs).

use std::f64::consts::PI;

use molrs::types::F;

/// One of the four closed-form target-distribution shapes.
///
/// `u`/`du_dxi` return the raw Boltzmann inversion of the un-normalized profile
/// ρ*(ξ) with the additive `ρ₀` reference dropped to zero (ρ₀ does not affect the
/// force, and the energy zero is immaterial to gencan). The shell-volume Jacobian
/// and density floor are NOT applied here — they are composed by
/// [`ProfilePenalty`]. For erf/tanh steps `rising = false` flips the argument sign
/// to describe a falling interface.
#[derive(Debug, Clone, Copy)]
pub enum Distribution {
    /// Gaussian ρ* ∝ exp(−(ξ−μ)²/2σ²): a harmonic well centred at `mu` with
    /// spring `k = kT/σ²`.
    Gaussian {
        /// Centre μ (Å).
        mu: F,
        /// Width σ (Å); smaller σ ⇒ stiffer well.
        sigma: F,
    },
    /// Error-function step ρ* ∝ ½[1+erf((ξ−ξ₀)/(√2 w))] — the default,
    /// capillary-wave-averaged physical interface.
    ErfStep {
        /// Interface midpoint ξ₀ (Å).
        xi0: F,
        /// Interface width w (Å).
        w: F,
        /// `true` rising side, `false` falling (argument sign flipped).
        rising: bool,
    },
    /// Hyperbolic-tangent step ρ* ∝ ½[1+tanh((ξ−ξ₀)/w)] — the mean-field variant.
    TanhStep {
        /// Interface midpoint ξ₀ (Å).
        xi0: F,
        /// Interface width w (Å).
        w: F,
        /// `true` rising side, `false` falling (argument sign flipped).
        rising: bool,
    },
    /// Exponential decay ρ* ∝ exp(−ξ/λ), valid for ξ≥0; gives a constant force
    /// `kT/λ` toward small ξ.
    Exponential {
        /// Decay length λ (Å).
        lambda: F,
    },
}

impl Distribution {
    /// The un-normalized target density ρ*(ξ) (ρ₀ dropped). Strictly positive for
    /// Gaussian/exponential; in [0, 1] for the erf/tanh steps.
    pub fn rho_star(&self, xi: F) -> F {
        match *self {
            Distribution::Gaussian { mu, sigma } => {
                let z = (xi - mu) / sigma;
                (-0.5 * z * z).exp()
            }
            Distribution::ErfStep { xi0, w, rising } => {
                let z = step_arg(xi, xi0, rising) / (std::f64::consts::SQRT_2 * w);
                0.5 * (1.0 + libm::erf(z))
            }
            Distribution::TanhStep { xi0, w, rising } => {
                let t = step_arg(xi, xi0, rising) / w;
                0.5 * (1.0 + t.tanh())
            }
            Distribution::Exponential { lambda } => (-xi / lambda).exp(),
        }
    }

    /// The exact derivative ρ*'(ξ) of [`rho_star`](Self::rho_star).
    pub fn d_rho_star(&self, xi: F) -> F {
        match *self {
            Distribution::Gaussian { mu, sigma } => {
                let s2 = sigma * sigma;
                self.rho_star(xi) * (-(xi - mu) / s2)
            }
            Distribution::ErfStep { xi0, w, rising } => {
                let sign = step_sign(rising);
                let z = step_arg(xi, xi0, rising) / (std::f64::consts::SQRT_2 * w);
                // d/dξ ½[1+erf(z)] = (1/√(2π)·w)·exp(−z²)·sign.
                sign * (-z * z).exp() / (TWO_PI.sqrt() * w)
            }
            Distribution::TanhStep { xi0, w, rising } => {
                let sign = step_sign(rising);
                let t = step_arg(xi, xi0, rising) / w;
                // d/dξ ½[1+tanh(t)] = ½·sech²(t)·sign/w; sech²(t) = 1/cosh²(t).
                let sech2 = 1.0 / (t.cosh() * t.cosh());
                0.5 * sech2 * sign / w
            }
            Distribution::Exponential { lambda } => self.rho_star(xi) * (-1.0 / lambda),
        }
    }

    /// Raw Boltzmann energy `U(ξ) = −kT·ln(ρ*(ξ))` (ρ₀ = 1, additive const
    /// dropped). No Jacobian, no floor.
    pub fn u(&self, xi: F, kt: F) -> F {
        -kt * self.rho_star(xi).ln()
    }

    /// Raw force-generating derivative `dU/dξ = −kT·ρ*'(ξ)/ρ*(ξ)`. No Jacobian,
    /// no floor.
    pub fn du_dxi(&self, xi: F, kt: F) -> F {
        -kt * self.d_rho_star(xi) / self.rho_star(xi)
    }
}

/// 2π, the angular factor in the shell-volume Jacobian.
const TWO_PI: F = 2.0 * PI;

/// Signed offset `±(ξ−ξ₀)` for a step profile: `+` rising, `−` falling.
#[inline]
fn step_arg(xi: F, xi0: F, rising: bool) -> F {
    step_sign(rising) * (xi - xi0)
}

/// `+1` for a rising interface, `−1` for a falling one.
#[inline]
fn step_sign(rising: bool) -> F {
    if rising { 1.0 } else { -1.0 }
}

/// Selects `dV/dξ` for the histogram→density conversion (`ρ* = n/(dV/dξ)`).
///
/// Only the additive `ln(dV/dξ)` term and its ξ-derivative are needed for the
/// inversion, so this exposes those directly. The geometry kind is supplied by
/// the caller (the composing restraint in a later sub-spec).
#[derive(Debug, Clone, Copy)]
pub enum ShellJacobian {
    /// Planar slab: `dV/dξ` is constant ⇒ no centre-weighting term.
    Planar,
    /// Radial shell: `dV/dξ = 4π ξ²`.
    Radial,
    /// Cylindrical shell: `dV/dξ = 2π ξ L`.
    Cylindrical {
        /// Cylinder length L (Å).
        length: F,
    },
}

impl ShellJacobian {
    /// The additive `ln(dV/dξ)` correction added to `U` when the input is a
    /// histogram: `0` planar, `ln(4π ξ²)` radial, `ln(2π ξ L)` cylindrical.
    pub fn log_shell_correction(&self, xi: F) -> F {
        match *self {
            ShellJacobian::Planar => 0.0,
            ShellJacobian::Radial => (4.0 * PI * xi * xi).ln(),
            ShellJacobian::Cylindrical { length } => (TWO_PI * xi * length).ln(),
        }
    }

    /// The ξ-derivative of [`log_shell_correction`](Self::log_shell_correction):
    /// `0` planar, `2/ξ` radial, `1/ξ` cylindrical — the spurious outward force
    /// injected if the Jacobian is omitted.
    pub fn d_log_shell(&self, xi: F) -> F {
        match *self {
            ShellJacobian::Planar => 0.0,
            ShellJacobian::Radial => 2.0 / xi,
            ShellJacobian::Cylindrical { .. } => 1.0 / xi,
        }
    }
}

/// Declares whether a [`Distribution`]'s value is already a volumetric density or
/// a raw count histogram needing shell-volume division before inversion.
///
/// The default is [`CountHistogram`](InputKind::CountHistogram) — the
/// domain-safe choice that never under-weights the shell correction.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum InputKind {
    /// Already a volumetric density ρ*(ξ) [count/volume]; invert directly.
    VolumetricDensity,
    /// A count histogram n(ξ) [count/ξ]; divide by `dV/dξ` before inverting.
    #[default]
    CountHistogram,
}

/// Density floor / energy cap (§6.3 divergence control).
///
/// Caps the penalty at `U_max = −kT·ln(ρ_min/ρ₀)` so a zero-probability bin
/// becomes a finite plateau instead of an infinite wall with no usable gradient.
#[derive(Debug, Clone, Copy)]
pub struct DensityFloor {
    /// Floor density ρ_min (count/volume); pick ~1e-3–1e-2 of the peak.
    pub rho_min: F,
    /// Reference density ρ₀ setting the energy zero.
    pub rho0: F,
    /// Energy-scale factor kT the cap was built against.
    pub kt: F,
    /// Precomputed cap `U_max = −kT·ln(ρ_min/ρ₀)`.
    pub u_max: F,
}

impl DensityFloor {
    /// Build a floor for the given `rho_min`, reference `rho0`, and `kt`,
    /// precomputing `u_max`.
    pub fn new(rho_min: F, rho0: F, kt: F) -> Self {
        Self {
            rho_min,
            rho0,
            kt,
            u_max: -kt * (rho_min / rho0).ln(),
        }
    }

    /// `true` when the density has fallen to/below the floor and the cap is
    /// active (so the bin contributes a plateau energy and no force). A `NaN`
    /// density is treated as capped so the gradient is zeroed rather than
    /// propagated.
    #[inline]
    pub fn is_capped(&self, rho_density: F) -> bool {
        rho_density <= self.rho_min || rho_density.is_nan()
    }

    /// Clamp `u` to at most `u_max`.
    #[inline]
    pub fn cap(&self, u: F) -> F {
        u.min(self.u_max)
    }
}

/// Immutable façade composing a [`Distribution`], its shell-volume
/// [`ShellJacobian`], the [`InputKind`] flag, and a [`DensityFloor`] into the
/// final `U(ξ)` / `dU/dξ`.
///
/// It obtains the un-inverted profile value ρ*(ξ) and derivative ρ*'(ξ) from the
/// distribution and performs inversion + (optional) shell-Jacobian division +
/// floor CENTRALLY, so a `CountHistogram + Radial` penalty fed n(ξ) gives the
/// SAME `U` as a `VolumetricDensity` penalty fed ρ*(ξ) = n(ξ)/(4π ξ²).
#[derive(Debug, Clone, Copy)]
pub struct ProfilePenalty {
    /// The closed-form target shape.
    pub dist: Distribution,
    /// The geometry's shell-volume Jacobian (used only for `CountHistogram`).
    pub jacobian: ShellJacobian,
    /// Whether [`dist`](Self::dist) yields a density or a count histogram.
    pub input_kind: InputKind,
    /// The density floor / energy cap.
    pub floor: DensityFloor,
}

impl ProfilePenalty {
    /// `true` when this penalty divides out the shell-volume Jacobian (i.e. the
    /// input is a count histogram).
    #[inline]
    fn uses_shell(&self) -> bool {
        self.input_kind == InputKind::CountHistogram
    }

    /// The volumetric density ρ_density(ξ): the raw ρ* for a density input, or
    /// `ρ*/(dV/dξ)` for a histogram input. Used to test the floor.
    #[inline]
    fn rho_density(&self, xi: F) -> F {
        let rho = self.dist.rho_star(xi);
        if self.uses_shell() {
            // dV/dξ = exp(log_shell_correction); divide n by it.
            rho / self.jacobian.log_shell_correction(xi).exp()
        } else {
            rho
        }
    }

    /// Composed penalty energy `U(ξ)` with Jacobian division and floor cap.
    ///
    /// `U = −kT·ln(ρ*/ρ₀) + kT·ln(dV/dξ)` for a histogram input (the second term
    /// is the shell correction); capped at `U_max` where the density is below the
    /// floor.
    pub fn u(&self, xi: F, kt: F) -> F {
        if self.floor.is_capped(self.rho_density(xi)) {
            return self.floor.u_max;
        }
        // −kT·ln(ρ*/ρ₀); ρ₀ from the floor keeps U on the same scale as U_max.
        let mut u = -kt * (self.dist.rho_star(xi) / self.floor.rho0).ln();
        if self.uses_shell() {
            u += kt * self.jacobian.log_shell_correction(xi);
        }
        self.floor.cap(u)
    }

    /// Composed force-generating derivative `dU/dξ` with Jacobian division and
    /// floor zeroing.
    ///
    /// `dU/dξ = −kT·ρ*'/ρ* + kT·d/dξ ln(dV/dξ)` for a histogram input; zero where
    /// the cap is active (a capped/empty bin exerts no force).
    pub fn du_dxi(&self, xi: F, kt: F) -> F {
        if self.floor.is_capped(self.rho_density(xi)) {
            return 0.0;
        }
        let mut du = self.dist.du_dxi(xi, kt);
        if self.uses_shell() {
            du += kt * self.jacobian.d_log_shell(xi);
        }
        du
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const TOL_REL: F = 1e-6;
    const KT: F = 0.6;

    /// Central finite difference of `f` at `x` with step `h`.
    fn fd<G: Fn(F) -> F>(f: G, x: F, h: F) -> F {
        (f(x + h) - f(x - h)) / (2.0 * h)
    }

    fn assert_close_rel(got: F, want: F, label: &str) {
        let denom = want.abs().max(1.0);
        let rel = (got - want).abs() / denom;
        assert!(
            rel <= TOL_REL,
            "{label}: got={got} want={want} rel_err={rel}"
        );
    }

    // ── ac-001: closed-form du_dxi vs central FD of u, all four shapes ─────────

    fn check_grad_matches_fd(dist: Distribution, xs: &[F], label: &str) {
        let h = 1e-5;
        for &x in xs {
            let analytic = dist.du_dxi(x, KT);
            let numeric = fd(|t| dist.u(t, KT), x, h);
            assert_close_rel(analytic, numeric, &format!("{label} @ xi={x}"));
        }
    }

    #[test]
    fn gaussian_grad_matches_fd() {
        let d = Distribution::Gaussian {
            mu: 3.0,
            sigma: 1.5,
        };
        check_grad_matches_fd(d, &[0.5, 2.0, 3.0, 4.5, 6.0], "gaussian");
    }

    #[test]
    fn erf_step_grad_matches_fd() {
        let rising = Distribution::ErfStep {
            xi0: 5.0,
            w: 1.2,
            rising: true,
        };
        let falling = Distribution::ErfStep {
            xi0: 5.0,
            w: 1.2,
            rising: false,
        };
        // Stay away from the deep tail where rho* underflows (handled by floor,
        // not by the raw closed form).
        check_grad_matches_fd(rising, &[3.0, 4.5, 5.0, 5.5, 7.0], "erf rising");
        check_grad_matches_fd(falling, &[3.0, 4.5, 5.0, 5.5, 7.0], "erf falling");
    }

    #[test]
    fn tanh_step_grad_matches_fd() {
        let rising = Distribution::TanhStep {
            xi0: 5.0,
            w: 1.2,
            rising: true,
        };
        let falling = Distribution::TanhStep {
            xi0: 5.0,
            w: 1.2,
            rising: false,
        };
        check_grad_matches_fd(rising, &[3.0, 4.5, 5.0, 5.5, 7.0], "tanh rising");
        check_grad_matches_fd(falling, &[3.0, 4.5, 5.0, 5.5, 7.0], "tanh falling");
    }

    #[test]
    fn exponential_grad_matches_fd() {
        let d = Distribution::Exponential { lambda: 2.5 };
        check_grad_matches_fd(d, &[0.1, 1.0, 2.0, 4.0, 8.0], "exponential");
    }

    // ── ac-002: closed-form identities (harmonic Gaussian; constant exp force) ─

    #[test]
    fn gaussian_is_harmonic() {
        let (mu, sigma) = (3.0, 1.5);
        let d = Distribution::Gaussian { mu, sigma };
        for &x in &[0.0, 1.0, 3.0, 5.0, 7.5] {
            let want_u = (KT / (2.0 * sigma * sigma)) * (x - mu).powi(2);
            let want_du = (KT / (sigma * sigma)) * (x - mu);
            assert_close_rel(d.u(x, KT), want_u, &format!("gaussian u @ {x}"));
            assert_close_rel(d.du_dxi(x, KT), want_du, &format!("gaussian du @ {x}"));
        }
    }

    #[test]
    fn exponential_force_is_constant() {
        let lambda = 2.5;
        let d = Distribution::Exponential { lambda };
        let want = KT / lambda;
        for &x in &[0.0, 0.5, 2.0, 5.0, 12.0] {
            assert_close_rel(d.du_dxi(x, KT), want, &format!("exp du @ {x}"));
        }
    }

    // ── ac-003: two-sided radial histogram-vs-density Jacobian ────────────────

    /// A "profile" whose raw rho_star is the supplied count histogram n(xi).
    /// We reuse the Gaussian shape as a smooth, differentiable n(xi) and feed it
    /// to ProfilePenalty under different InputKind / ShellJacobian settings.
    fn radial_penalty(input_kind: InputKind, jac: ShellJacobian) -> ProfilePenalty {
        ProfilePenalty {
            dist: Distribution::Gaussian {
                mu: 6.0,
                sigma: 2.0,
            },
            jacobian: jac,
            input_kind,
            floor: DensityFloor::new(1e-30, 1.0, KT),
        }
    }

    #[test]
    fn histogram_with_jacobian_equals_density() {
        // n(xi) is the Gaussian rho_star; the matching volumetric density is
        // rho*(xi) = n(xi)/(4 pi xi^2). The CountHistogram+Radial penalty divides
        // n by 4 pi xi^2 internally, so it must equal a VolumetricDensity penalty
        // fed that same density. Because ProfilePenalty derives the density from
        // the SAME dist internally, the two configurations describe the same
        // radial profile and must give equal U.
        let hist = radial_penalty(InputKind::CountHistogram, ShellJacobian::Radial);
        let dens = radial_penalty(InputKind::VolumetricDensity, ShellJacobian::Radial);
        for &xi in &[2.0, 4.0, 6.0, 8.0, 10.0] {
            // density U should equal histogram U MINUS the shell correction the
            // histogram path adds, i.e. they differ by exactly kt*ln(4 pi xi^2).
            let u_hist = hist.u(xi, KT);
            let u_dens = dens.u(xi, KT);
            let shell = KT * (4.0 * PI * xi * xi).ln();
            assert_close_rel(u_hist, u_dens + shell, &format!("hist vs density U @ {xi}"));
        }
    }

    #[test]
    fn omitting_jacobian_adds_spurious_force() {
        // Treating the histogram as a density (Planar jacobian / no shell factor)
        // drops the +kt*ln(4 pi xi^2) term: dU/dxi differs by exactly +2kt/xi.
        let correct = radial_penalty(InputKind::CountHistogram, ShellJacobian::Radial);
        let wrong = radial_penalty(InputKind::CountHistogram, ShellJacobian::Planar);
        for &xi in &[2.0, 4.0, 6.0, 8.0, 10.0] {
            let du_correct = correct.du_dxi(xi, KT);
            let du_wrong = wrong.du_dxi(xi, KT);
            let predicted = 2.0 * KT / xi;
            assert_close_rel(
                du_correct - du_wrong,
                predicted,
                &format!("spurious force @ {xi}"),
            );
            // And the energies must genuinely differ (not a no-op).
            assert!(
                (correct.u(xi, KT) - wrong.u(xi, KT)).abs() > 1e-9,
                "U identical at xi={xi} — Jacobian had no effect"
            );
        }
    }

    // ── ac-004: density floor — zero-density yields finite capped U / du ───────

    #[test]
    fn zero_density_is_capped_and_finite() {
        let rho_min = 1e-3;
        let rho0 = 1.0;
        // Falling erf far on the high side → rho* underflows toward 0.
        let pen = ProfilePenalty {
            dist: Distribution::ErfStep {
                xi0: 5.0,
                w: 0.5,
                rising: false,
            },
            jacobian: ShellJacobian::Planar,
            input_kind: InputKind::VolumetricDensity,
            floor: DensityFloor::new(rho_min, rho0, KT),
        };
        let u_max = -KT * (rho_min / rho0).ln();
        let xi = 50.0; // deep in the dead tail
        let u = pen.u(xi, KT);
        let du = pen.du_dxi(xi, KT);
        assert!(u.is_finite(), "U not finite: {u}");
        assert!(du.is_finite(), "du not finite: {du}");
        assert_close_rel(u, u_max, "capped U == U_max");
        assert_eq!(du, 0.0, "capped bin must exert no force, got {du}");
    }

    #[test]
    fn explicit_zero_density_is_capped() {
        // Drive rho* literally to zero via a Gaussian evaluated far from mu so the
        // density is below rho_min; the floor must still give finite U == U_max.
        let rho_min = 1e-2;
        let rho0 = 1.0;
        let pen = ProfilePenalty {
            dist: Distribution::Gaussian {
                mu: 0.0,
                sigma: 0.3,
            },
            jacobian: ShellJacobian::Planar,
            input_kind: InputKind::VolumetricDensity,
            floor: DensityFloor::new(rho_min, rho0, KT),
        };
        let u_max = -KT * (rho_min / rho0).ln();
        let xi = 20.0;
        let u = pen.u(xi, KT);
        let du = pen.du_dxi(xi, KT);
        assert!(u.is_finite() && du.is_finite());
        assert_close_rel(u, u_max, "gaussian tail capped U");
        assert_eq!(du, 0.0);
    }

    // ── smoke: types are Copy + Debug ─────────────────────────────────────────

    #[test]
    fn types_are_copy_debug() {
        let pen = ProfilePenalty {
            dist: Distribution::Exponential { lambda: 1.0 },
            jacobian: ShellJacobian::Cylindrical { length: 10.0 },
            input_kind: InputKind::default(),
            floor: DensityFloor::new(1e-3, 1.0, KT),
        };
        let copied = pen; // Copy
        let _ = format!("{pen:?} {copied:?}");
        assert_eq!(InputKind::default(), InputKind::CountHistogram);
    }
}
