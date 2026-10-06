//! Truncated-Newton line search (`tnls` in Packmol `gencan.f`).

use molrs::op::types::F;

use crate::Objective;
use crate::eval::EvalMode;

/// Truncated-Newton line search result.
pub(super) struct LsResult {
    pub f: F,
    pub fcnt: usize,
    pub gcnt: usize,
    pub inform: i32,
}

/// Reusable buffers for TN line search (`tnls` in Packmol `gencan.f`).
pub(super) struct TnLsScratch {
    pub(super) xret: Vec<F>,
    pub(super) gret: Vec<F>,
    pub(super) xplus: Vec<F>,
    pub(super) xtmp: Vec<F>,
    pub(super) gplus: Vec<F>,
}

impl TnLsScratch {
    pub(super) fn new(n: usize) -> Self {
        Self {
            xret: vec![0.0; n],
            gret: vec![0.0; n],
            xplus: vec![0.0; n],
            xtmp: vec![0.0; n],
            gplus: vec![0.0; n],
        }
    }

    pub(super) fn ensure_len(&mut self, n: usize) {
        if self.xret.len() != n {
            self.xret.resize(n, 0.0);
        }
        if self.gret.len() != n {
            self.gret.resize(n, 0.0);
        }
        if self.xplus.len() != n {
            self.xplus.resize(n, 0.0);
        }
        if self.xtmp.len() != n {
            self.xtmp.resize(n, 0.0);
        }
        if self.gplus.len() != n {
            self.gplus.resize(n, 0.0);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn tn_linesearch(
    nind: usize,
    ind: &[usize],
    n: usize,
    x: &[F],
    g: &[F],
    d: &[F],
    l: &[F],
    u: &[F],
    f0: F,
    amax: F,
    rbdtype: i32,
    rbdind: usize,
    nint: F,
    next: F,
    mininterp: usize,
    maxextrap: usize,
    fmin: F,
    maxfc: usize,
    mut fcnt: usize,
    mut gcnt: usize,
    gamma: F,
    beta: F,
    _sterel: F,
    _steabs: F,
    sigma1: F,
    sigma2: F,
    epsrel: F,
    epsabs: F,
    scratch: &mut TnLsScratch,
    obj: &mut dyn Objective,
) -> LsResult {
    scratch.ensure_len(n);
    let nind = nind.min(ind.len());
    let xret = &mut scratch.xret;
    let gret = &mut scratch.gret;
    let xplus = &mut scratch.xplus;
    let xtmp = &mut scratch.xtmp;
    let gplus = &mut scratch.gplus;
    xret.copy_from_slice(x);
    gret.copy_from_slice(g);
    let mut fret = f0;

    let mut gplus_valid = false;

    // gtd = <g,d> in the free-variable subspace.
    let mut gtd = 0.0 as F;
    for &ii in &ind[..nind] {
        gtd += g[ii] * d[ii];
    }

    // First trial alpha = min(1, amax)
    let mut alpha = (1.0 as F).min(amax);
    xplus.copy_from_slice(x);
    for &ii in &ind[..nind] {
        xplus[ii] = x[ii] + alpha * d[ii];
    }
    if alpha == amax && rbdtype != 0 {
        if rbdtype == 1 {
            xplus[rbdind] = l[rbdind];
        } else {
            xplus[rbdind] = u[rbdind];
        }
    }

    let mut fplus = obj.evaluate(xplus, EvalMode::FOnly, None).f_total;
    fcnt += 1;

    let mut do_extrap = false;

    // Decide between extrapolation and interpolation.
    if amax > 1.0 {
        if fplus <= f0 + gamma * alpha * gtd {
            obj.evaluate(xplus, EvalMode::GradientOnly, Some(gplus));
            gcnt += 1;
            gplus_valid = true;

            let mut gptd = 0.0 as F;
            for &ii in &ind[..nind] {
                gptd += gplus[ii] * d[ii];
            }

            if gptd < beta * gtd {
                do_extrap = true;
            } else {
                xret.copy_from_slice(xplus);
                gret.copy_from_slice(gplus);
                fret = fplus;
                return LsResult {
                    f: fret,
                    fcnt,
                    gcnt,
                    inform: 0,
                };
            }
        }
    } else if fplus < f0 {
        do_extrap = true;
    }

    // ------------------------------------------------------------------
    // Extrapolation
    // ------------------------------------------------------------------
    if do_extrap {
        let mut extrap = 0usize;

        loop {
            if fplus <= fmin {
                xret.copy_from_slice(xplus);
                fret = fplus;
                if extrap != 0 || amax <= 1.0 {
                    obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
                    gcnt += 1;
                } else if gplus_valid {
                    gret.copy_from_slice(gplus);
                }
                return LsResult {
                    f: fret,
                    fcnt,
                    gcnt,
                    inform: 4,
                };
            }

            if fcnt >= maxfc {
                xret.copy_from_slice(xplus);
                fret = fplus;
                if extrap != 0 || amax <= 1.0 {
                    obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
                    gcnt += 1;
                } else if gplus_valid {
                    gret.copy_from_slice(gplus);
                }
                return LsResult {
                    f: fret,
                    fcnt,
                    gcnt,
                    inform: 8,
                };
            }

            if extrap >= maxextrap {
                xret.copy_from_slice(xplus);
                fret = fplus;
                if extrap != 0 || amax <= 1.0 {
                    obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
                    gcnt += 1;
                } else if gplus_valid {
                    gret.copy_from_slice(gplus);
                }
                return LsResult {
                    f: fret,
                    fcnt,
                    gcnt,
                    inform: 7,
                };
            }

            let atmp = if alpha < amax && next * alpha > amax {
                amax
            } else {
                next * alpha
            };

            xtmp.copy_from_slice(x);
            for &ii in &ind[..nind] {
                xtmp[ii] = x[ii] + atmp * d[ii];
            }
            if atmp == amax && rbdtype != 0 {
                if rbdtype == 1 {
                    xtmp[rbdind] = l[rbdind];
                } else {
                    xtmp[rbdind] = u[rbdind];
                }
            }
            if atmp > amax {
                for &ii in &ind[..nind] {
                    xtmp[ii] = xtmp[ii].clamp(l[ii], u[ii]);
                }
            }

            if alpha > amax {
                let mut samep = true;
                for &ii in &ind[..nind] {
                    if (xtmp[ii] - xplus[ii]).abs() > (epsrel * xplus[ii].abs()).max(epsabs) {
                        samep = false;
                        break;
                    }
                }

                if samep {
                    xret.copy_from_slice(xplus);
                    fret = fplus;
                    if extrap != 0 || amax <= 1.0 {
                        obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
                        gcnt += 1;
                    } else if gplus_valid {
                        gret.copy_from_slice(gplus);
                    }
                    return LsResult {
                        f: fret,
                        fcnt,
                        gcnt,
                        inform: 0,
                    };
                }
            }

            let ftmp = obj.evaluate(xtmp, EvalMode::FOnly, None).f_total;
            fcnt += 1;

            if ftmp < fplus {
                alpha = atmp;
                fplus = ftmp;
                xplus.copy_from_slice(xtmp);
                gplus_valid = false;
                extrap += 1;
                continue;
            }

            xret.copy_from_slice(xplus);
            fret = fplus;
            if extrap != 0 || amax <= 1.0 {
                obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
                gcnt += 1;
            } else if gplus_valid {
                gret.copy_from_slice(gplus);
            } else {
                obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
                gcnt += 1;
            }
            return LsResult {
                f: fret,
                fcnt,
                gcnt,
                inform: 0,
            };
        }
    }

    // ------------------------------------------------------------------
    // Interpolation
    // ------------------------------------------------------------------
    let mut interp = 0usize;
    loop {
        if fplus <= fmin {
            xret.copy_from_slice(xplus);
            fret = fplus;
            obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
            gcnt += 1;
            return LsResult {
                f: fret,
                fcnt,
                gcnt,
                inform: 4,
            };
        }

        if fcnt >= maxfc {
            if fplus < f0 {
                xret.copy_from_slice(xplus);
                fret = fplus;
                obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
                gcnt += 1;
            }
            return LsResult {
                f: fret,
                fcnt,
                gcnt,
                inform: 8,
            };
        }

        if fplus <= f0 + gamma * alpha * gtd {
            xret.copy_from_slice(xplus);
            fret = fplus;
            obj.evaluate(xret, EvalMode::GradientOnly, Some(gret));
            gcnt += 1;
            return LsResult {
                f: fret,
                fcnt,
                gcnt,
                inform: 0,
            };
        }

        interp += 1;
        if alpha < sigma1 {
            alpha /= nint;
        } else {
            let denom = 2.0 * (fplus - f0 - alpha * gtd);
            let atmp = if denom != 0.0 {
                (-gtd * alpha * alpha) / denom
            } else {
                alpha / nint
            };
            if atmp < sigma1 || atmp > sigma2 * alpha {
                alpha /= nint;
            } else {
                alpha = atmp;
            }
        }

        xplus.copy_from_slice(x);
        for &ii in &ind[..nind] {
            xplus[ii] = x[ii] + alpha * d[ii];
        }

        fplus = obj.evaluate(xplus, EvalMode::FOnly, None).f_total;
        fcnt += 1;

        let mut samep = true;
        for &ii in &ind[..nind] {
            if (alpha * d[ii]).abs() > (epsrel * x[ii].abs()).max(epsabs) {
                samep = false;
                break;
            }
        }
        if interp >= mininterp && samep {
            return LsResult {
                f: fret,
                fcnt,
                gcnt,
                inform: 6,
            };
        }
    }
}
