//! GENCAN optimizer — faithful Rust port of `gencan.f` and `pgencan.f90`.
//!
//! Reference: Birgin & Martinez, Comp.Opt.Appl. 23:101-125, 2002.

use molrs::types::F;
pub mod cg;
pub mod entry;
pub mod phases;
pub mod solver;
pub mod spg;

use crate::constraints::EvalMode;
use crate::numerics::{numeric_controls, positive_norm_floor};
use crate::objective::Objective;

/// Parameters for the GENCAN call (matches `easygencan` defaults from `pgencan.f90`).
pub struct GencanParams {
    pub epsgpsn: F,
    pub maxit: usize,
    pub maxfc: usize,
    pub delmin: F,
    pub iprint: i32,
    pub ncomp: usize,
}

impl Default for GencanParams {
    fn default() -> Self {
        Self {
            epsgpsn: 1.0e-6,
            maxit: 20,
            maxfc: 200,     // 10 * maxit
            delmin: 1.0e-2, // Packmol easygencan default (gencan.f: delmin = 1.d-2)
            iprint: 0,
            ncomp: 50,
        }
    }
}

/// Result of a GENCAN run.
pub struct GencanResult {
    pub f: F,
    pub gpsupn: F,
    pub iter: usize,
    pub fcnt: usize,
    pub gcnt: usize,
    pub cgcnt: usize,
    /// 0=converged(eucl), 1=converged(sup), 2=noFprogress, 3=noGprogress,
    /// 4=fSmall, 7=maxIter, 8=maxFeval, <0=error
    pub inform: i32,
}

/// Reusable work buffers for repeated `pgencan` calls.
pub struct GencanWorkspace {
    g: Vec<F>,
    ind: Vec<usize>,
    d: Vec<F>,
    s: Vec<F>,
    y: Vec<F>,
    cg_scratch: cg::CgScratch,
    spg_scratch: spg::SpgScratch,
    tnls_scratch: TnLsScratch,
}

impl GencanWorkspace {
    pub fn new() -> Self {
        Self {
            g: Vec::new(),
            ind: Vec::new(),
            d: Vec::new(),
            s: Vec::new(),
            y: Vec::new(),
            cg_scratch: cg::CgScratch::new(0),
            spg_scratch: spg::SpgScratch::new(0),
            tnls_scratch: TnLsScratch::new(0),
        }
    }

    fn ensure_len(&mut self, n: usize) {
        if self.g.len() != n {
            self.g.resize(n, 0.0);
        }
        if self.d.len() != n {
            self.d.resize(n, 0.0);
        }
        if self.s.len() != n {
            self.s.resize(n, 0.0);
        }
        if self.y.len() != n {
            self.y.resize(n, 0.0);
        }
        if self.ind.capacity() < n {
            self.ind.reserve(n - self.ind.capacity());
        }
    }
}

impl Default for GencanWorkspace {
    fn default() -> Self {
        Self::new()
    }
}

/// Entry point — mirrors `pgencan.f90` → `easygencan` → `gencan`.
///
/// The variable layout is:
///   x[0..3N]   = COM positions (free molecules)
///   x[3N..6N]  = Euler angles (free molecules)
///
/// Bounds: COM variables are unbounded; Euler angles may be bounded by
/// `constrain_rotation` constraints (Packmol pgencan.f90).
pub fn pgencan(
    x: &mut [F],
    obj: &mut dyn Objective,
    params: &GencanParams,
    precision: F,
    workspace: &mut GencanWorkspace,
) -> GencanResult {
    let n = x.len();

    let mut l = vec![0.0 as F; n];
    let mut u = vec![0.0 as F; n];
    obj.bounds(&mut l, &mut u);

    gencan(x, &l, &u, obj, params, precision, workspace)
}

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

/// Truncated-Newton line search result.
struct LsResult {
    pub f: F,
    pub fcnt: usize,
    pub gcnt: usize,
    pub inform: i32,
}

/// Reusable buffers for TN line search (`tnls` in Packmol `gencan.f`).
struct TnLsScratch {
    xret: Vec<F>,
    gret: Vec<F>,
    xplus: Vec<F>,
    xtmp: Vec<F>,
    gplus: Vec<F>,
}

impl TnLsScratch {
    fn new(n: usize) -> Self {
        Self {
            xret: vec![0.0; n],
            gret: vec![0.0; n],
            xplus: vec![0.0; n],
            xtmp: vec![0.0; n],
            gplus: vec![0.0; n],
        }
    }

    fn ensure_len(&mut self, n: usize) {
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
fn tn_linesearch(
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

#[cfg(test)]
mod tests {
    //! Unit-level guards for the GENCAN optimizer (`src/gencan/`).
    //!
    //! These tests pin the optimizer's behaviour on synthetic
    //! quadratics with known analytic minima — fast (< 50 ms each) and
    //! decoupled from the rest of the pack pipeline.
    //!
    //! The `Objective` trait is the only contract `gencan` depends on, so we
    //! provide a hand-rolled `Quadratic` impl rather than poking at
    //! `PackContext` internals. `fdist`/`frest` are reported as 0 so the
    //! Packmol-style early-exit check (`fdist < precision && frest < precision`)
    //! never fires when `precision = 0.0`; gencan then has to converge on its
    //! own gpsupn / maxit criterion.

    use crate::constraints::{EvalMode, EvalOutput};
    use crate::gencan::{GencanParams, GencanWorkspace, gencan, pgencan};
    use crate::objective::Objective;
    use molrs::types::F;

    /// f(x) = 0.5 · Σ (xᵢ − μᵢ)²; ∇f = (x − μ); minimum at x = μ, f = 0.
    struct Quadratic {
        mu: Vec<F>,
        ncf: usize,
        ncg: usize,
    }

    impl Quadratic {
        fn new(mu: Vec<F>) -> Self {
            Self { mu, ncf: 0, ncg: 0 }
        }
    }

    impl Objective for Quadratic {
        fn evaluate(&mut self, x: &[F], mode: EvalMode, gradient: Option<&mut [F]>) -> EvalOutput {
            debug_assert_eq!(x.len(), self.mu.len());

            let mut f = 0.0;
            for (xi, mu_i) in x.iter().zip(self.mu.iter()) {
                let dx = xi - mu_i;
                f += 0.5 * dx * dx;
            }

            match mode {
                EvalMode::FOnly => {
                    self.ncf += 1;
                }
                EvalMode::FAndGradient | EvalMode::GradientOnly | EvalMode::RestMol => {
                    self.ncf += 1;
                    self.ncg += 1;
                    if let Some(g) = gradient {
                        for ((gi, xi), mu_i) in g.iter_mut().zip(x.iter()).zip(self.mu.iter()) {
                            *gi = xi - mu_i;
                        }
                    }
                }
            }

            EvalOutput {
                f_total: f,
                fdist_max: 0.0,
                frest_max: 0.0,
            }
        }

        fn fdist(&self) -> F {
            0.0
        }

        fn frest(&self) -> F {
            0.0
        }

        fn ncf(&self) -> usize {
            self.ncf
        }

        fn ncg(&self) -> usize {
            self.ncg
        }

        fn reset_eval_counters(&mut self) {
            self.ncf = 0;
            self.ncg = 0;
        }
    }

    // ── unconstrained CG path ──────────────────────────────────────────────────

    /// Drive `pgencan` on a 3-D positive-definite quadratic with no bounds.
    /// The CG inner solver should hit the analytic minimum (μ) within the
    /// gpsupn tolerance after a handful of outer iterations. If a future
    /// refactor breaks `cg_solve`'s descent direction or the outer
    /// projected-gradient bookkeeping, x will not converge to μ.
    #[test]
    fn pgencan_unconstrained_quadratic_converges_to_minimum() {
        let mu = vec![3.0, -2.0, 1.5];
        let mut obj = Quadratic::new(mu.clone());
        let mut x = vec![0.0; mu.len()];
        let params = GencanParams::default();
        let mut ws = GencanWorkspace::new();

        // precision = 0.0 → packmolprecision (`fdist < 0 && frest < 0`) is
        // never satisfied; gencan terminates only on gpsupn/maxit.
        let result = pgencan(&mut x, &mut obj, &params, 0.0, &mut ws);

        assert!(
            result.inform == 0 || result.inform == 1,
            "gencan did not converge: inform={}, gpsupn={}, iter={}",
            result.inform,
            result.gpsupn,
            result.iter
        );
        for i in 0..x.len() {
            let err = (x[i] - mu[i]).abs();
            assert!(
                err < 1e-4,
                "x[{i}] = {} expected ≈ {} (err={err}, gpsupn={})",
                x[i],
                mu[i],
                result.gpsupn
            );
        }
        assert!(result.f < 1e-8, "f = {} expected ≈ 0", result.f);
    }

    // ── bound-constrained SPG path ─────────────────────────────────────────────

    /// Drive `gencan` directly with explicit bounds so the projected
    /// gradient hits a wall: minimum of f(x) = 0.5(x − 5)² over x ∈ [0, 1]
    /// is at x = 1 (clamped, projected gradient = 0 at the upper bound).
    /// Exercises the SPG / projection branch of gencan that the
    /// unconstrained test above never reaches.
    #[test]
    fn gencan_bound_constrained_quadratic_clamps_to_boundary() {
        let mu = vec![5.0];
        let mut obj = Quadratic::new(mu.clone());
        let mut x = vec![0.0];
        let l = vec![0.0];
        let u = vec![1.0];
        let params = GencanParams::default();
        let mut ws = GencanWorkspace::new();

        let result = gencan(&mut x, &l, &u, &mut obj, &params, 0.0, &mut ws);

        assert!(
            result.inform >= 0,
            "gencan errored: inform={}",
            result.inform
        );
        let err = (x[0] - 1.0).abs();
        assert!(
            err < 1e-6,
            "expected x clamped to upper bound 1.0, got {} (err={err}, gpsupn={})",
            x[0],
            result.gpsupn
        );
        let expected_f = 0.5 * (1.0 - 5.0_f64).powi(2);
        assert!(
            (result.f - expected_f).abs() < 1e-6,
            "expected f = {expected_f}, got {}",
            result.f
        );
    }

    /// Bilateral clamp: minimum at μ = −3, search bounded on [0, 4] so the
    /// unique constrained minimum is at the lower bound x = 0. Mirrors the
    /// `_clamps_to_boundary` test on the opposite face.
    #[test]
    fn gencan_clamps_to_lower_bound() {
        let mut obj = Quadratic::new(vec![-3.0]);
        let mut x = vec![2.0];
        let l = vec![0.0];
        let u = vec![4.0];
        let params = GencanParams::default();
        let mut ws = GencanWorkspace::new();

        let result = gencan(&mut x, &l, &u, &mut obj, &params, 0.0, &mut ws);

        assert!(result.inform >= 0);
        assert!(
            (x[0] - 0.0).abs() < 1e-6,
            "expected x clamped to 0.0, got {}",
            x[0]
        );
    }

    // ── the seam markers the GENCAN stage declares ─────────────────────────────

    /// Owner-side half of acceptance ac-008: what `GencanStage` declares on the
    /// stage seam belongs here, not in `stage::tests` (which knows only fakes).
    ///
    /// `Placed::None` because the stage seeds its own placements with `initial()`
    /// — it needs nothing placed on entry; `Placed::All` because it returns with
    /// every free molecule placed. The two markers are what a pipeline checks
    /// before chaining anything after this stage.
    #[test]
    fn gencan_stage_requires_none_guarantees_all() {
        use crate::gencan::solver::{GencanSettings, GencanStage};
        use crate::{Placed, Stage};

        // The declarations are construction-time constants: an empty system is
        // enough to read them, and using one keeps this test off the algorithm.
        let stage = GencanStage::new(GencanSettings::default(), Vec::new(), None, 0, 0);

        assert_eq!(
            stage.requires().placed,
            Placed::None,
            "GENCAN places the molecules itself, so it requires nothing placed"
        );
        assert_eq!(
            stage.guarantees().placed,
            Placed::All,
            "GENCAN returns with every free molecule placed"
        );
    }
}
