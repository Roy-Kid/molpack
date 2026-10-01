//! GENCAN loop and the projected-gradient helpers it owns.

use molrs::types::F;

use super::linesearch::tn_linesearch;
use super::{GencanParams, GencanResult, GencanWorkspace, cg, spg};
use crate::eval::EvalMode;
use crate::numerics::{numeric_controls, positive_norm_floor};
use crate::objective::Objective;

/// Main GENCAN loop.
/// Port of `gencan.f` subroutine with Packmol-specific additions.
pub fn gencan(
    x: &mut [F],
    l: &[F],
    u: &[F],
    obj: &mut dyn Objective,
    params: &GencanParams,
    precision: F,
    workspace: &mut GencanWorkspace,
) -> GencanResult {
    let n = x.len();
    workspace.ensure_len(n);

    // Constants (from easygencan parameters section)
    const INFREL: F = 1.0e20;
    const INFABS: F = F::MAX;
    const BETA_LS: F = 0.5;
    const GAMMA: F = 1.0e-4;
    const THETA: F = 1.0e-6;
    const SIGMA1: F = 0.1;
    const SIGMA2: F = 0.9;
    const MAXEXTRAP: usize = 100;
    const MININTERP: usize = 4;
    const NINT: F = 2.0;
    const NEXT: F = 2.0;
    const ETA: F = 0.9;
    const LSPGMA: F = 1.0e10;
    const LSPGMI: F = 1.0e-10;
    const FMIN: F = 1.0e-5;
    const EPSGPEN: F = 0.0;
    let numeric = numeric_controls();

    let cgepsi = 0.1 as F;
    let cgepsf = 1.0e-5 as F;
    let cggpnf = F::max(1.0e-4, params.epsgpsn);
    let epsnqmp = 1.0e-4 as F;
    let maxitnqmp = 5usize;
    let epsnfp = 0.0 as F;
    let maxitnfp = params.maxit;
    let maxitngp = 1000usize;

    // Project initial point
    for i in 0..n {
        x[i] = x[i].clamp(l[i], u[i]);
    }

    // Initial function value + gradient.
    let mut fcnt = 0usize;
    let g = &mut workspace.g;
    g.fill(0.0);
    let mut f = obj.evaluate(x, EvalMode::FAndGradient, Some(g)).f_total;

    // Packmol behavior: check convergence before counting this first eval.
    if converged(obj, precision) {
        return GencanResult {
            f,
            gpsupn: 0.0,
            iter: 0,
            fcnt,
            gcnt: 0,
            cgcnt: 0,
            inform: 0,
        };
    }
    fcnt += 1;
    let mut gcnt = 1usize;
    let mut cgcnt = 0usize;

    // Compute xnorm
    let mut xnorm = x.iter().map(|xi| xi * xi).sum::<F>().sqrt();

    let ind = &mut workspace.ind;
    let cg_scratch = &mut workspace.cg_scratch;
    let spg_scratch = &mut workspace.spg_scratch;
    let tnls_scratch = &mut workspace.tnls_scratch;

    // Compute projected gradient
    let (mut gpsupn, mut gpeucn2, mut gieucn2, mut nind) =
        projected_gradient_info(n, x, g.as_slice(), l, u, ind);

    // CG epsilon scaling
    let (acgeps, bcgeps) = gp_ieee_signal(gpsupn, cgepsf, cgepsi, cggpnf);

    // Track initial projected gradient for kappa computation in cgmaxit
    // (Packmol gp_ieee_signal2: cgscre=2, uses sup-norm)
    let gpsupn0 = gpsupn;

    let mut iter = 0usize;
    let mut inform = 7i32; // default: max iterations

    // Trust radius — computed fresh at the start of each TN iteration (see below).
    // Fortran: iter==1 → max(delmin, 0.1*xnorm); iter>1 → max(delmin, 10*sqrt(sts)).
    let mut delta: F;

    // BB spectral step
    let mut sts = 0.0 as F;
    let mut sty = 0.0 as F;
    let ometa2 = (1.0 - ETA).powi(2);

    // No-progress tracking
    let mut fprev = INFABS;
    let mut bestprog = 0.0 as F;
    let mut itnfp = 0usize;
    let mut lastgpns = vec![INFABS; maxitngp];

    // Working vectors
    let d = &mut workspace.d;
    let s = &mut workspace.s;
    let y = &mut workspace.y;

    // Main loop
    loop {
        // Packmol behavior: recompute precision test with computef at each iteration.
        if converged(obj, precision) {
            break;
        }

        if gpeucn2 <= EPSGPEN * EPSGPEN {
            inform = 0;
            break;
        }
        // Check convergence: sup-norm of projected gradient
        if gpsupn <= params.epsgpsn {
            inform = 1;
            break;
        }

        // No function progress
        let currprog = fprev - f;
        bestprog = bestprog.max(currprog);
        if currprog <= epsnfp * bestprog {
            itnfp += 1;
            if itnfp >= maxitnfp {
                inform = 2;
                break;
            }
        } else {
            itnfp = 0;
        }

        // No gradient progress
        let gpnmax = lastgpns.iter().copied().fold(0.0 as F, F::max);
        lastgpns[iter % maxitngp] = gpeucn2;
        if gpeucn2 >= gpnmax {
            inform = 3;
            break;
        }

        if f <= FMIN {
            inform = 4;
            break;
        }
        if iter >= params.maxit {
            inform = 7;
            break;
        }
        if fcnt >= params.maxfc {
            inform = 8;
            break;
        }

        // New iteration
        iter += 1;
        fprev = f;

        // Save x → s, g → y
        s.copy_from_slice(x);
        y.copy_from_slice(g.as_slice());

        if gieucn2 <= ometa2 * gpeucn2 {
            // SPG iteration: abandon current face
            let lamspg = if iter == 1 || sty <= 0.0 {
                F::max(1.0, xnorm) / gpeucn2.sqrt().max(positive_norm_floor())
            } else {
                sts / sty
            };
            let lamspg = lamspg.clamp(LSPGMI, LSPGMA);

            let spg_res = spg::spgls(
                n,
                x,
                g.as_slice(),
                l,
                u,
                lamspg,
                f,
                NINT,
                MININTERP,
                FMIN,
                params.maxfc,
                fcnt,
                GAMMA,
                SIGMA1,
                SIGMA2,
                numeric.sterel,
                numeric.steabs,
                numeric.epsrel,
                numeric.epsabs,
                spg_scratch,
                obj,
            );
            f = spg_res.f;
            fcnt = spg_res.fcnt;
            x.copy_from_slice(&spg_scratch.xtrial);

            if spg_res.inform < 0 {
                inform = spg_res.inform;
                break;
            }

            obj.evaluate(x, EvalMode::GradientOnly, Some(g.as_mut_slice()));
            gcnt += 1;
        } else {
            // TN iteration: compute Newton direction via CG

            // Compute trust-region radius (Fortran gencan.f lines 2120-2128):
            //   iter==1: delta = max(delmin, 0.1 * max(1, xnorm))
            //   iter>1:  delta = max(delmin, 10 * sqrt(sts))
            delta = if iter == 1 {
                F::max(params.delmin, 0.1 * F::max(1.0, xnorm))
            } else {
                F::max(params.delmin, 10.0 * sts.sqrt())
            };

            let cgeps = compute_cgeps(gpsupn, acgeps, bcgeps, cgepsf, cgepsi);

            // Packmol gp_ieee_signal2 formula (cgscre=2, nearlyq=false):
            //   kappa = clamp(log10(gpsupn/gpsupn0) / log10(epsgpsn/gpsupn0), 0, 1)
            //   cgmaxit = min(20, (1-kappa)*max(1, 10*log10(nind)) + kappa*nind)
            let cgmaxit = {
                let mut kappa = (gpsupn / gpsupn0).log10() / (params.epsgpsn / gpsupn0).log10();
                kappa = F::max(0.0, F::min(1.0, kappa));

                let nind_f = nind as F;
                let base = (1.0 - kappa) * F::max(1.0, 10.0 * nind_f.log10()) + kappa * nind_f;
                usize::min(20, base as usize)
            };

            let cg_res = cg::cg_solve(
                nind,
                ind,
                n,
                x,
                g.as_slice(),
                delta,
                l,
                u,
                cgeps,
                epsnqmp,
                maxitnqmp,
                cgmaxit,
                false, // nearlyq = .false. in packmol easygencan defaults
                1,     // trtype=1 (sup-norm)
                THETA,
                numeric.sterel,
                numeric.steabs,
                numeric.epsrel,
                numeric.epsabs,
                INFREL,
                INFABS,
                d,
                cg_scratch,
                obj,
            );
            cgcnt += cg_res.iter;

            // Compute maximum feasible step along d (packmol gencan.f lines 2204-2225).
            let mut amax = INFABS;
            let mut rbdtype = 0i32;
            let mut rbdind = if nind > 0 { ind[0] } else { 0 };

            if cg_res.inform == 2 {
                amax = 1.0;
                rbdtype = cg_res.rbdtype;
                rbdind = cg_res.rbdind.unwrap_or(rbdind);
            } else {
                for &ii in &ind[..nind] {
                    if d[ii] > 0.0 {
                        let amaxx = (u[ii] - x[ii]) / d[ii];
                        if amaxx < amax {
                            amax = amaxx;
                            rbdind = ii;
                            rbdtype = 2;
                        }
                    } else if d[ii] < 0.0 {
                        let amaxx = (l[ii] - x[ii]) / d[ii];
                        if amaxx < amax {
                            amax = amaxx;
                            rbdind = ii;
                            rbdtype = 1;
                        }
                    }
                }
            }

            // TN line search (full port of tnls behavior).
            let ls_res = tn_linesearch(
                nind,
                ind,
                n,
                x,
                g.as_slice(),
                d,
                l,
                u,
                f,
                amax,
                rbdtype,
                rbdind,
                NINT,
                NEXT,
                MININTERP,
                MAXEXTRAP,
                FMIN,
                params.maxfc,
                fcnt,
                gcnt,
                GAMMA,
                BETA_LS,
                numeric.sterel,
                numeric.steabs,
                SIGMA1,
                SIGMA2,
                numeric.epsrel,
                numeric.epsabs,
                tnls_scratch,
                obj,
            );
            f = ls_res.f;
            fcnt = ls_res.fcnt;
            gcnt = ls_res.gcnt;
            x.copy_from_slice(&tnls_scratch.xret);
            g.copy_from_slice(&tnls_scratch.gret);

            if ls_res.inform < 0 {
                inform = ls_res.inform;
                break;
            }
            inform = ls_res.inform;

            // packmol behavior: if tnls stops with inform=6, discard TN step and force SPG.
            if ls_res.inform == 6 {
                let lamspg = if iter == 1 || sty <= 0.0 {
                    F::max(1.0, xnorm) / gpeucn2.sqrt().max(positive_norm_floor())
                } else {
                    sts / sty
                };
                let lamspg = lamspg.clamp(LSPGMI, LSPGMA);

                let spg_res = spg::spgls(
                    n,
                    x,
                    g.as_slice(),
                    l,
                    u,
                    lamspg,
                    f,
                    NINT,
                    MININTERP,
                    FMIN,
                    params.maxfc,
                    fcnt,
                    GAMMA,
                    SIGMA1,
                    SIGMA2,
                    numeric.sterel,
                    numeric.steabs,
                    numeric.epsrel,
                    numeric.epsabs,
                    spg_scratch,
                    obj,
                );
                f = spg_res.f;
                fcnt = spg_res.fcnt;
                x.copy_from_slice(&spg_scratch.xtrial);

                if spg_res.inform < 0 {
                    inform = spg_res.inform;
                    break;
                }

                let infotmp = spg_res.inform;
                obj.evaluate(x, EvalMode::GradientOnly, Some(g.as_mut_slice()));
                gcnt += 1;
                inform = infotmp;
            }
        }

        // Adjust to bounds near machine precision (packmol gencan.f lines 2363-2371).
        for i in 0..n {
            if x[i] <= l[i] + (numeric.epsrel * l[i].abs()).max(numeric.epsabs) {
                x[i] = l[i];
            } else if x[i] >= u[i] - (numeric.epsrel * u[i].abs()).max(numeric.epsabs) {
                x[i] = u[i];
            }
        }

        // Update x norm.
        xnorm = x.iter().map(|xi| xi * xi).sum::<F>().sqrt();

        // Update BB steplength: sts = (x-s)^T(x-s), sty = (x-s)^T(g-y)
        sts = 0.0;
        sty = 0.0;
        for i in 0..n {
            let ds = x[i] - s[i];
            let dy = g[i] - y[i];
            sts += ds * ds;
            sty += ds * dy;
        }

        // Update projected gradient info (reuses pre-allocated ind buffer)
        let pg = projected_gradient_info(n, x, g.as_slice(), l, u, ind);
        gpsupn = pg.0;
        gpeucn2 = pg.1;
        gieucn2 = pg.2;
        nind = pg.3;
    }

    // No extra compute_f here — packer.rs calls compute_f with unscaled radii
    // immediately after pgencan returns, so this would be redundant.

    GencanResult {
        f,
        gpsupn,
        iter,
        fcnt,
        gcnt,
        cgcnt,
        inform,
    }
}

/// Compute projected gradient info into pre-allocated `ind` buffer.
/// Returns (gpsupn, gpeucn2, gieucn2, nind).
fn projected_gradient_info(
    n: usize,
    x: &[F],
    g: &[F],
    l: &[F],
    u: &[F],
    ind: &mut Vec<usize>,
) -> (F, F, F, usize) {
    let mut gpsupn = 0.0 as F;
    let mut gpeucn2 = 0.0 as F;
    let mut gieucn2 = 0.0 as F;

    ind.clear();
    for i in 0..n {
        let gpi = (x[i] - g[i]).clamp(l[i], u[i]) - x[i];
        gpsupn = gpsupn.max(gpi.abs());
        gpeucn2 += gpi * gpi;
        if x[i] > l[i] && x[i] < u[i] {
            gieucn2 += gpi * gpi;
            ind.push(i);
        }
    }

    (gpsupn, gpeucn2, gieucn2, ind.len())
}

/// Packmol precision check (`packmolprecision` in `pgencan.f90`):
/// recompute objective-side violations and test `fdist/frest`.
fn converged(obj: &dyn Objective, precision: F) -> bool {
    obj.fdist() < precision && obj.frest() < precision
}

/// Compute scaling for CG epsilon.
fn gp_ieee_signal(gpsupn: F, cgepsf: F, cgepsi: F, cggpnf: F) -> (F, F) {
    if gpsupn > 0.0 {
        let acgeps = (cgepsf / cgepsi).log10() / (cggpnf / gpsupn).log10();
        let bcgeps = cgepsi.log10() - acgeps * gpsupn.log10();
        (acgeps, bcgeps)
    } else {
        (0.0, cgepsf)
    }
}

fn compute_cgeps(gpsupn: F, acgeps: F, bcgeps: F, cgepsf: F, cgepsi: F) -> F {
    let cgeps = (10.0 as F).powf(acgeps * gpsupn.log10() + bcgeps);
    cgeps.clamp(cgepsf, cgepsi)
}
