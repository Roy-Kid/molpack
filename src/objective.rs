//! Objective function and gradient computation.
//! Exact port of `computef.f90`, `computeg.f90`, `fparc.f90`, `gparc.f90`.

use crate::cell::{index_cell, setcell};
use crate::constraints::{EvalMode, EvalOutput};
use crate::context::{ATOM_FLAG_FIXED, ATOM_FLAG_SHORT, NONE_IDX, PackContext};
use crate::euler::{compcart, eulerrmat, eulerrmat_derivatives};
use molrs::types::F;
#[cfg(feature = "rayon")]
use rayon::prelude::*;

#[derive(Clone, Copy)]
enum ExpandMode {
    F,
    G,
    FG,
}

/// Precomputed PBC constants so the inner pair loops can wrap without
/// a per-axis division. `inv_length[k] = 1/length[k]` when PBC is active on
/// axis k, else 0. `any_active` short-circuits the entire wrap in the common
/// non-PBC workload.
#[derive(Clone, Copy)]
struct PbcConstants {
    length: [F; 3],
    inv_length: [F; 3],
    any_active: bool,
}

#[inline(always)]
fn pbc_constants(sys: &PackContext) -> PbcConstants {
    let length = sys.pbc_length;
    let mut inv_length = [0.0 as F; 3];
    let mut any_active = false;
    for k in 0..3 {
        if sys.pbc_periodic[k] && length[k] > 0.0 {
            inv_length[k] = 1.0 / length[k];
            any_active = true;
        }
    }
    PbcConstants {
        length,
        inv_length,
        any_active,
    }
}

#[inline(always)]
fn pbc_wrap_delta(dx: F, dy: F, dz: F, pbc: &PbcConstants) -> (F, F, F) {
    if !pbc.any_active {
        return (dx, dy, dz);
    }
    (
        dx - (dx * pbc.inv_length[0]).round() * pbc.length[0],
        dy - (dy * pbc.inv_length[1]).round() * pbc.length[1],
        dz - (dz * pbc.inv_length[2]).round() * pbc.length[2],
    )
}

// Parallel pair evaluation is user-selected via
// `Molpack::parallel_pair_eval(true)` → `PackContext::parallel_pair_eval`.
// The library does not attempt to auto-detect when rayon pays off: the
// crossover is workload-shaped (call frequency, work-per-call,
// movebad proportion) rather than something inferrable from a single
// per-pack metric like `active_cells.len()`.

/// Per-atom hot-path state pulled out once before the inner `jcart` loop
/// inside `fparc` / `gparc` / `fgparc` (and the rayon variants `fparc_stats`
/// and `fgparc_into`).
///
/// Before this type, each of those four kernels had the same 15-line
/// prologue: read `atom_props[icart]`, derive `fixed_i`, `use_short_i`,
/// `shrad_i`, `shscl_i`, and cache `any_fixed_atoms` / `any_short_radius`
/// as `has_fixed` / `has_short`. Drift between the four copies would
/// silently change kernel behaviour — centralising it here means one
/// edit site. The struct is `#[inline(always)]`-constructed so the call
/// compiles to the same loads the inlined prologue produced (verified
/// by the `pair_kernel` bench — see the commit message for Fix 2).
#[derive(Clone, Copy)]
struct AtomHotState {
    props: crate::context::AtomProps,
    /// Cached `sys.any_short_radius` — guards the cold short-radius fetch.
    has_short: bool,
    /// Atom `i` itself is fixed — used only together with the `j`-fixed
    /// check to short-circuit the pair.
    fixed_i: bool,
    /// Atom `i` itself opts into the short-radius penalty — one half of
    /// the `use_short_i || use_short_j` gate.
    use_short_i: bool,
    shrad_i: F,
    shscl_i: F,
}

impl AtomHotState {
    #[inline(always)]
    fn load(icart: usize, sys: &PackContext) -> Self {
        let props = sys.atom_props[icart];
        let has_short = sys.any_short_radius;
        let has_fixed = sys.any_fixed_atoms;
        let fixed_i = has_fixed && (props.flags & ATOM_FLAG_FIXED != 0);
        let use_short_i = has_short && (props.flags & ATOM_FLAG_SHORT != 0);
        // The short-radius `Vec`s are cold when `any_short_radius` is
        // false, so skip the load entirely — the branch below folds away
        // when `has_short` is a compile-time constant in the caller
        // (rare) or branch-predicts perfectly otherwise.
        let shrad_i = if has_short {
            sys.short_radius[icart]
        } else {
            0.0
        };
        let shscl_i = if has_short {
            sys.short_radius_scale[icart]
        } else {
            0.0
        };
        Self {
            props,
            has_short,
            fixed_i,
            use_short_i,
            shrad_i,
            shscl_i,
        }
    }
}

/// Compute objective function value and update `sys.xcart`.
/// Port of `computef.f90`.
///
/// `x` layout: [COM₀..COMₙ (3N)] ++ [euler₀..eulerₙ (3N)]
///
/// Takes the geometry-cache fast path when `x` and the cell grid are identical
/// to the last expansion: the Cartesian rebuild and linked-cell rebuild are
/// skipped and only the constraint + pair kernels re-run on the stored state.
/// The packer hits this path on every outer iteration when it re-evaluates at
/// unscaled radii after `pgencan` returns, since radius mutation does not
/// invalidate the cache key.
pub fn compute_f(x: &[F], sys: &mut PackContext) -> F {
    sys.debug_assert_atom_props_sync();
    sys.increment_ncf();
    sys.fdist = 0.0;
    sys.frest = 0.0;

    if matches_cached_geometry(x, sys) {
        let mut f = accumulate_constraint_values_from_xcart(sys);
        if !sys.init1 {
            f += accumulate_pair_f(sys);
            f += accumulate_collective_f(sys);
        }
        return f;
    }

    if !sys.init1 {
        sys.resetcells();
    }

    let mut f = expand_molecules(x, sys, ExpandMode::F);

    if sys.init1 {
        update_cached_geometry(x, sys);
        return f;
    }

    f += accumulate_pair_f(sys);
    f += accumulate_collective_f(sys);
    update_cached_geometry(x, sys);

    f
}

/// Compute objective function value and gradient in one pass over geometry/state.
/// This avoids rebuilding Cartesian coordinates and cell lists twice.
///
/// Takes the geometry-cache fast path when `x` and the cell grid match the
/// last expansion — the Cartesian rebuild and linked-cell rebuild are skipped.
pub fn compute_fg(x: &[F], sys: &mut PackContext, g: &mut [F]) -> F {
    sys.debug_assert_atom_props_sync();
    sys.increment_ncf();
    sys.increment_ncg();
    sys.fdist = 0.0;
    sys.frest = 0.0;
    sys.work.gxcar.fill([0.0; 3]);

    if matches_cached_geometry(x, sys) {
        let mut f = accumulate_constraint_values_and_gradients_from_xcart(sys);
        if !sys.init1 {
            f += accumulate_pair_fg(sys);
            f += accumulate_collective_fg(sys);
        }
        project_cartesian_gradient(x, sys, g);
        return f;
    }

    if !sys.init1 {
        sys.resetcells();
    }

    let mut f = expand_molecules(x, sys, ExpandMode::FG);

    if !sys.init1 {
        f += accumulate_pair_fg(sys);
        f += accumulate_collective_fg(sys);
    }

    update_cached_geometry(x, sys);
    project_cartesian_gradient(x, sys, g);

    f
}

/// Contribution of one *i*–*j* atom pair to the packing objective.
///
/// `grad` is the force to **add** to atom *i* and **subtract** from atom *j*
/// (Newton's third law); it folds the main and short-range terms together and
/// stays `[0; 3]` whenever the pair doesn't overlap or `GRAD` is `false`.
/// `overlap` is `datom < tol` — callers use it to skip the gradient write for
/// the common non-overlapping pair. `violation` is `tol_ini - datom`, only
/// meaningful (and only computed) when `VIOLATION` is `true`.
#[derive(Clone, Copy)]
struct PairContribution {
    energy: F,
    grad: [F; 3],
    violation: F,
    overlap: bool,
}

/// The single source of truth for per-pair packing physics: the distance
/// penalty, its gradient, the short-range correction, and the overlap
/// violation. Every neighbor-walk variant (serial/parallel, value-only/fused,
/// move-tracking or not) calls this — they differ only in *where* they apply
/// the result, never in the formula.
///
/// Returns `None` when the pair is skipped (same molecule, or both fixed). The
/// two const generics let each caller drop the work it doesn't need at zero
/// runtime cost: `GRAD` removes the force math for value-only passes,
/// `VIOLATION` removes the `fdist` term for the gradient-only pass.
///
/// Per-atom reads go through `sys.atom_props` — an AoS mirror packing the hot
/// fields into one cache line per atom (see [`AtomProps`]).
#[inline(always)]
fn pair_term<const GRAD: bool, const VIOLATION: bool>(
    hot: &AtomHotState,
    xi: [F; 3],
    jcart: usize,
    sys: &PackContext,
    pbc: &PbcConstants,
) -> Option<PairContribution> {
    let props_j = sys.atom_props[jcart];
    // Skip same molecule.
    if hot.props.ibmol == props_j.ibmol && hot.props.ibtype == props_j.ibtype {
        return None;
    }
    // Skip two fixed atoms (only matters if the system has any at all).
    if hot.fixed_i && (props_j.flags & ATOM_FLAG_FIXED != 0) {
        return None;
    }

    let xj = sys.xcart[jcart];
    let (dx, dy, dz) = pbc_wrap_delta(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2], pbc);
    let datom = dx * dx + dy * dy + dz * dz;
    let rsum = hot.props.radius + props_j.radius;
    let tol = rsum * rsum;

    let mut energy: F = 0.0;
    let mut grad = [0.0; 3];
    let overlap = datom < tol;

    if overlap {
        let penalty = datom - tol;
        let scale = hot.props.fscale * props_j.fscale;
        energy += scale * penalty * penalty;
        if GRAD {
            let dtemp = scale * 4.0 * penalty;
            grad[0] += dtemp * dx;
            grad[1] += dtemp * dy;
            grad[2] += dtemp * dz;
        }
        if hot.has_short && (hot.use_short_i || (props_j.flags & ATOM_FLAG_SHORT != 0)) {
            let short_rsum = hot.shrad_i + sys.short_radius[jcart];
            let short_tol = short_rsum * short_rsum;
            if datom < short_tol {
                let short_penalty = datom - short_tol;
                let mut sr_scale = (hot.shscl_i * sys.short_radius_scale[jcart]).sqrt();
                sr_scale *= (tol * tol) / (short_tol * short_tol);
                let sr_pair_scale = scale * sr_scale;
                energy += sr_pair_scale * short_penalty * short_penalty;
                if GRAD {
                    let dtemp2 = sr_pair_scale * 4.0 * short_penalty;
                    grad[0] += dtemp2 * dx;
                    grad[1] += dtemp2 * dy;
                    grad[2] += dtemp2 * dz;
                }
            }
        }
    }

    let violation = if VIOLATION {
        let rsum_ini = hot.props.radius_ini + props_j.radius_ini;
        let tol_ini = rsum_ini * rsum_ini;
        tol_ini - datom
    } else {
        0.0
    };

    Some(PairContribution {
        energy,
        grad,
        violation,
        overlap,
    })
}

/// Atom-pair distance penalty function.
/// Port of `fparc.f90`.
///
/// Returns `(penalty_sum, fdist_max)`. The caller aggregates `fdist_max`
/// across pair traversal and updates `sys.fdist` once at the end, so the
/// inner loop keeps a data dependency on a local register instead of an
/// `&mut PackContext` field.
///
/// The per-atom reads go through `sys.atom_props` — an AoS mirror that
/// packs the ten or so hot fields into one cache line per atom (see
/// [`AtomProps`]).
#[inline(always)]
fn fparc(icart: usize, first_jcart: u32, sys: &mut PackContext, pbc: &PbcConstants) -> (F, F) {
    let mut result = 0.0;
    let mut local_fdist: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);
    let move_flag = sys.move_flag;

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        if let Some(c) = pair_term::<false, true>(&hot, xi, jcart, sys, pbc) {
            result += c.energy;
            if c.violation > local_fdist {
                local_fdist = c.violation;
            }
            if move_flag {
                if c.violation > sys.fdist_atom[icart] {
                    sys.fdist_atom[icart] = c.violation;
                }
                if c.violation > sys.fdist_atom[jcart] {
                    sys.fdist_atom[jcart] = c.violation;
                }
            }
        }
        jcart_id = next;
    }

    (result, local_fdist)
}

/// Compute gradient `g` from current system state.
/// Port of `computeg.f90`.
pub fn compute_g(x: &[F], sys: &mut PackContext, g: &mut [F]) {
    sys.debug_assert_atom_props_sync();
    sys.increment_ncg();
    // Zero Cartesian gradient
    sys.work.gxcar.fill([0.0; 3]);

    if matches_cached_geometry(x, sys) {
        accumulate_constraint_gradient_from_xcart(sys);
        if !sys.init1 {
            accumulate_pair_g(sys);
            let _ = accumulate_collective_fg(sys);
        }
        project_cartesian_gradient(x, sys, g);
        return;
    }

    if !sys.init1 {
        sys.resetcells();
    }

    expand_molecules(x, sys, ExpandMode::G);

    if !sys.init1 {
        accumulate_pair_g(sys);
        let _ = accumulate_collective_fg(sys);
    }

    update_cached_geometry(x, sys);
    project_cartesian_gradient(x, sys, g);
}

/// Expand all active molecules into `sys.xcart`, applying either constraint
/// values (`computef`) or constraint gradients (`computeg`), and rebuilding the
/// cell list.
///
/// One phase-structured implementation shared by the serial and rayon paths.
/// **Phase A** rebuilds each molecule's Cartesian coordinates (`eulerrmat` +
/// `compcart`) and evaluates its box/restraint constraints, writing **disjoint**
/// per-atom `xcart` / `gxcar` (and, under `move_flag`, `frest_atom`) slots — so it
/// has no cross-task write conflict and parallelises cleanly. **Phase B** replays
/// the linked-cell insertion serially in the identical molecule/atom order (that
/// insertion mutates a shared list and is inherently serial). The *only*
/// difference between a serial and a parallel run is whether Phase A iterates
/// with `iter` or `par_iter` — see [`expand_reduce`]. There is no second
/// algorithm, and the resulting cell lists are bit-identical either way.
fn expand_molecules(x: &[F], sys: &mut PackContext, mode: ExpandMode) -> F {
    // Cheap serial pass: one descriptor per active molecule (pure index
    // arithmetic, no trig), reusing the persistent workspace buffer.
    let mut descs = std::mem::take(&mut sys.work.mol_descs);
    fill_active_mol_descs(&mut descs, sys);

    // Move the per-atom write targets out of `sys` so Phase A can share
    // `&PackContext` immutably while each molecule writes its own disjoint slots
    // through raw pointers (constraint reads — coor, restraints, iratom — stay on
    // the shared borrow). Sound for the serial `iter` path too: it just runs the
    // per-molecule closures sequentially.
    let mut xcart = std::mem::take(&mut sys.xcart);
    let mut gxcar = std::mem::take(&mut sys.work.gxcar);
    let mut frest_atom = std::mem::take(&mut sys.frest_atom);
    let scale = sys.scale;
    let scale2 = sys.scale2;
    let move_flag = sys.move_flag;
    // Per-atom `frest_atom` bookkeeping is movebad-only; that path always runs
    // serially, so it never races even though the parallel branch could write it.
    let parallel = sys.parallel_pair_eval && !move_flag;

    #[derive(Clone, Copy)]
    struct Slots {
        xcart: *mut [F; 3],
        gxcar: *mut [F; 3],
        frest_atom: *mut F,
    }
    // SAFETY: each molecule owns a disjoint, contiguous `icart` range, so no two
    // tasks ever write the same slot.
    unsafe impl Send for Slots {}
    unsafe impl Sync for Slots {}
    impl Slots {
        // Accessors take `self` so a closure that calls them captures the whole
        // (`Sync`) wrapper, not its bare `*mut` fields (which would make the
        // closure non-`Sync` and break `par_iter`).
        /// # Safety
        /// `i` must be in bounds and owned solely by the calling task.
        #[inline(always)]
        unsafe fn xcart_at(self, i: usize) -> *mut [F; 3] {
            unsafe { self.xcart.add(i) }
        }
        /// # Safety
        /// `i` must be in bounds and owned solely by the calling task.
        #[inline(always)]
        unsafe fn gxcar_at<'a>(self, i: usize) -> &'a mut [F; 3] {
            unsafe { &mut *self.gxcar.add(i) }
        }
        /// # Safety
        /// `i` must be in bounds and owned solely by the calling task.
        #[inline(always)]
        unsafe fn frest_atom_at(self, i: usize) -> *mut F {
            unsafe { self.frest_atom.add(i) }
        }
    }
    let slots = Slots {
        xcart: xcart.as_mut_ptr(),
        gxcar: gxcar.as_mut_ptr(),
        frest_atom: frest_atom.as_mut_ptr(),
    };

    // Run Phase A in an inner scope so the shared `&PackContext` borrow (held by
    // the closure) is released before we move the buffers back into `sys`.
    let (f_total, frest_max) = {
        let sys_ro: &PackContext = sys;
        // Phase A, per molecule — independent, so it runs under `iter`/`par_iter`.
        let body = |&(itype, icart0, ilubar, ilugan): &(usize, usize, usize, usize)| -> (F, F) {
            let (v1, v2, v3) = eulerrmat(x[ilugan], x[ilugan + 1], x[ilugan + 2]);
            let xcm = [x[ilubar], x[ilubar + 1], x[ilubar + 2]];
            let idbase = sys_ro.idfirst[itype];
            let na = sys_ro.natoms[itype];
            let mut f_local: F = 0.0;
            let mut frest_local: F = 0.0;
            for iatom in 0..na {
                let icart = icart0 + iatom;
                let pos = compcart(&xcm, &sys_ro.coor[idbase + iatom], &v1, &v2, &v3);
                // SAFETY: `icart` is owned by this molecule alone.
                unsafe {
                    *slots.xcart_at(icart) = pos;
                }
                let start = sys_ro.iratom_offsets[icart];
                let end = sys_ro.iratom_offsets[icart + 1];
                if start == end {
                    continue;
                }
                // Value (F / FG): same `.f` call order as the legacy serial loop.
                if matches!(mode, ExpandMode::F | ExpandMode::FG) {
                    let mut fplus = 0.0;
                    for &irest in &sys_ro.iratom_data[start..end] {
                        fplus += sys_ro.restraints[irest].f(&pos, scale, scale2);
                    }
                    f_local += fplus;
                    if fplus > frest_local {
                        frest_local = fplus;
                    }
                    if move_flag {
                        // SAFETY: disjoint per-atom slot.
                        unsafe {
                            *slots.frest_atom_at(icart) += fplus;
                        }
                    }
                }
                // Gradient (G / FG): accumulate into this atom's own slot.
                if matches!(mode, ExpandMode::G | ExpandMode::FG) {
                    // SAFETY: disjoint slot, as above.
                    let gc = unsafe { slots.gxcar_at(icart) };
                    for &irest in &sys_ro.iratom_data[start..end] {
                        let _ = sys_ro.restraints[irest].fg(&pos, scale, scale2, gc);
                    }
                }
            }
            (f_local, frest_local)
        };
        expand_reduce(&descs, &body, parallel)
    };

    // Restore the moved-out buffers before the (serial) cell-insertion phase.
    sys.xcart = xcart;
    sys.work.gxcar = gxcar;
    sys.frest_atom = frest_atom;
    if frest_max > sys.frest {
        sys.frest = frest_max;
    }

    // Phase B: serial linked-cell insertion in the identical order, so the lists
    // are bit-identical no matter how Phase A iterated.
    if !sys.init1 {
        for &(itype, icart0, _, _) in &descs {
            let na = sys.natoms[itype];
            for iatom in 0..na {
                let icart = icart0 + iatom;
                let pos = sys.xcart[icart];
                insert_atom_in_cell(icart, &pos, sys);
            }
        }
    }

    sys.work.mol_descs = descs;
    f_total
}

/// Run the independent per-molecule Phase-A `body` over `descs` and reduce to
/// `(f_total, frest_max)`. This is the single place where the serial and
/// parallel expand paths diverge: `par_iter` when `parallel`, plain `iter`
/// otherwise — the only `iter`-vs-`par_iter` switch in the expand pass.
#[cfg(feature = "rayon")]
fn expand_reduce<B>(descs: &[(usize, usize, usize, usize)], body: &B, parallel: bool) -> (F, F)
where
    B: Fn(&(usize, usize, usize, usize)) -> (F, F) + Sync,
{
    use rayon::prelude::*;
    if parallel {
        descs
            .par_iter()
            .map(body)
            .reduce(|| (0.0 as F, 0.0 as F), |a, b| (a.0 + b.0, a.1.max(b.1)))
    } else {
        descs
            .iter()
            .map(body)
            .fold((0.0 as F, 0.0 as F), |a, b| (a.0 + b.0, a.1.max(b.1)))
    }
}

/// Non-rayon build: the expand pass is always serial.
#[cfg(not(feature = "rayon"))]
fn expand_reduce<B>(descs: &[(usize, usize, usize, usize)], body: &B, _parallel: bool) -> (F, F)
where
    B: Fn(&(usize, usize, usize, usize)) -> (F, F),
{
    descs
        .iter()
        .map(body)
        .fold((0.0 as F, 0.0 as F), |a, b| (a.0 + b.0, a.1.max(b.1)))
}

/// Fill `descs` with one `(itype, icart0, ilubar, ilugan)` descriptor per
/// **active** molecule (clearing any previous contents first), mirroring the
/// expand/project loops' `comptype` skipping and `ilubar`/`ilugan` compaction.
/// `ilugan` (the Euler-angle offset) starts at `ntotmol*3`, which is correct for
/// both full and per-type compact `x` (`SwapState::set_type` shrinks `ntotmol`
/// to the phase's molecule count).
fn fill_active_mol_descs(descs: &mut Vec<(usize, usize, usize, usize)>, sys: &PackContext) {
    descs.clear();
    let mut ilubar = 0usize;
    let mut ilugan = sys.ntotmol * 3;
    let mut icart = 0usize;
    for itype in 0..sys.ntype {
        if !sys.comptype[itype] {
            icart += sys.nmols[itype] * sys.natoms[itype];
            continue;
        }
        for _ in 0..sys.nmols[itype] {
            descs.push((itype, icart, ilubar, ilugan));
            icart += sys.natoms[itype];
            ilubar += 3;
            ilugan += 3;
        }
    }
}

#[inline(always)]
fn accumulate_constraint_value(icart: usize, pos: &[F; 3], sys: &mut PackContext) -> F {
    let mut fplus = 0.0;
    let start = sys.iratom_offsets[icart];
    let end = sys.iratom_offsets[icart + 1];
    for &irest in &sys.iratom_data[start..end] {
        fplus += sys.restraints[irest].f(pos, sys.scale, sys.scale2);
    }
    if fplus > sys.frest {
        sys.frest = fplus;
    }
    if sys.move_flag {
        sys.frest_atom[icart] += fplus;
    }
    fplus
}

#[inline(always)]
fn accumulate_constraint_gradient(icart: usize, pos: &[F; 3], sys: &mut PackContext) {
    let start = sys.iratom_offsets[icart];
    let end = sys.iratom_offsets[icart + 1];
    let scale = sys.scale;
    let scale2 = sys.scale2;
    let gc = &mut sys.work.gxcar[icart];
    for &irest in &sys.iratom_data[start..end] {
        // fg returns the penalty value too; discard it — only gradient accumulation matters here
        let _ = sys.restraints[irest].fg(pos, scale, scale2, gc);
    }
}

#[inline]
fn accumulate_constraint_gradient_from_xcart(sys: &mut PackContext) {
    let mut icart = 0usize;

    for itype in 0..sys.ntype {
        if !sys.comptype[itype] {
            icart += sys.nmols[itype] * sys.natoms[itype];
            continue;
        }

        for _imol in 0..sys.nmols[itype] {
            for _iatom in 0..sys.natoms[itype] {
                let pos = sys.xcart[icart];
                accumulate_constraint_gradient(icart, &pos, sys);
                icart += 1;
            }
        }
    }
}

/// F-only counterpart to [`accumulate_constraint_gradient_from_xcart`]. Walks
/// `sys.xcart` (assumed current from a prior expansion) and re-applies each
/// restraint's function value, returning the accumulated penalty. Used by the
/// `compute_f` cache fast path when the Cartesian expansion can be skipped.
#[inline]
fn accumulate_constraint_values_from_xcart(sys: &mut PackContext) -> F {
    let mut f = 0.0;
    let mut icart = 0usize;

    for itype in 0..sys.ntype {
        if !sys.comptype[itype] {
            icart += sys.nmols[itype] * sys.natoms[itype];
            continue;
        }

        for _imol in 0..sys.nmols[itype] {
            for _iatom in 0..sys.natoms[itype] {
                let pos = sys.xcart[icart];
                f += accumulate_constraint_value(icart, &pos, sys);
                icart += 1;
            }
        }
    }

    f
}

/// Combined F+G counterpart to [`accumulate_constraint_gradient_from_xcart`].
/// Used by the `compute_fg` cache fast path.
#[inline]
fn accumulate_constraint_values_and_gradients_from_xcart(sys: &mut PackContext) -> F {
    let mut f = 0.0;
    let mut icart = 0usize;

    for itype in 0..sys.ntype {
        if !sys.comptype[itype] {
            icart += sys.nmols[itype] * sys.natoms[itype];
            continue;
        }

        for _imol in 0..sys.nmols[itype] {
            for _iatom in 0..sys.natoms[itype] {
                let pos = sys.xcart[icart];
                f += accumulate_constraint_value(icart, &pos, sys);
                accumulate_constraint_gradient(icart, &pos, sys);
                icart += 1;
            }
        }
    }

    f
}

/// Starting `icart` of free type `itype` — prefix sum of `nmols·natoms` over
/// preceding types. Cheap (`ntype` is small) and used only by collective terms.
#[inline]
fn type_icart_start(sys: &PackContext, itype: usize) -> usize {
    let mut start = 0usize;
    for t in 0..itype {
        start += sys.nmols[t] * sys.natoms[t];
    }
    start
}

/// Sum of all collective (group-level) restraint penalties. Each restraint is
/// evaluated once over **all copies** of its type (gathered from `sys.xcart`).
/// Read-only — used by the value-only `compute_f` paths.
fn accumulate_collective_f(sys: &PackContext) -> F {
    if sys.collective.is_empty() {
        return 0.0;
    }
    let (scale, scale2) = (sys.scale, sys.scale2);
    let mut total = 0.0;
    for (itype, r) in &sys.collective {
        let itype = *itype;
        if !sys.comptype[itype] {
            continue;
        }
        let start = type_icart_start(sys, itype);
        let len = sys.nmols[itype] * sys.natoms[itype];
        if len == 0 {
            continue;
        }
        total += r.f(&sys.xcart[start..start + len], scale, scale2);
    }
    total
}

/// Value + gradient of all collective restraints. The coupled gradient is
/// scattered into `sys.work.gxcar` so the subsequent `project_cartesian_gradient`
/// maps it onto each molecule's COM/Euler DOF (monatomic species → identity on z).
/// Returns the summed penalty value. Used by `compute_fg` (value) and
/// `compute_g` (gradient; value discarded).
fn accumulate_collective_fg(sys: &mut PackContext) -> F {
    if sys.collective.is_empty() {
        return 0.0;
    }
    let (scale, scale2) = (sys.scale, sys.scale2);
    // Move the list out so `sys.xcart` (read) and `sys.work.gxcar` (write) are
    // borrowed as disjoint fields without aliasing the whole `sys`.
    let collective = std::mem::take(&mut sys.collective);
    let mut total = 0.0;
    for (itype, r) in &collective {
        let itype = *itype;
        if !sys.comptype[itype] {
            continue;
        }
        let start = type_icart_start(sys, itype);
        let len = sys.nmols[itype] * sys.natoms[itype];
        if len == 0 {
            continue;
        }
        let coords: Vec<[F; 3]> = sys.xcart[start..start + len].to_vec();
        let mut grads = vec![[0.0 as F; 3]; len];
        total += r.fg(&coords, scale, scale2, &mut grads);
        for (k, g) in grads.iter().enumerate() {
            let gc = &mut sys.work.gxcar[start + k];
            gc[0] += g[0];
            gc[1] += g[1];
            gc[2] += g[2];
        }
    }
    sys.collective = collective;
    total
}

#[inline(always)]
fn insert_atom_in_cell(icart: usize, pos: &[F; 3], sys: &mut PackContext) {
    let cell = setcell(
        pos,
        &sys.pbc_min,
        &sys.pbc_length,
        &sys.cell_length,
        &sys.ncells,
        &sys.pbc_periodic,
    );
    let icell = index_cell(&cell, &sys.ncells);
    sys.latomnext[icart] = sys.latomfirst[icell];
    sys.latomfirst[icell] = icart as u32;

    if sys.empty_cell[icell] {
        sys.empty_cell[icell] = false;
        sys.active_cells.push(icell);
        sys.lcellnext[icell] = sys.lcellfirst;
        sys.lcellfirst = icell as u32;
    }

    // `ibtype` / `ibmol` used to be written here on every eval; they are
    // now set once in `Molpack::pack` and mirrored into `atom_props`.
}

#[inline(always)]
fn accumulate_pair_f(sys: &mut PackContext) -> F {
    #[cfg(feature = "rayon")]
    if sys.parallel_pair_eval && !sys.move_flag {
        let (f, fdist_max) = accumulate_pair_f_parallel(sys);
        sys.fdist = sys.fdist.max(fdist_max);
        return f;
    }

    let pbc = pbc_constants(sys);
    let mut f = 0.0;
    let mut fdist_local: F = 0.0;
    let mut icell_id = sys.lcellfirst;
    while icell_id != NONE_IDX {
        let icell = icell_id as usize;
        let neighbors = sys.neighbor_cells_f[icell];

        let mut icart_id = sys.latomfirst[icell];
        while icart_id != NONE_IDX {
            let icart = icart_id as usize;
            let (df, dfd) = fparc(icart, sys.latomnext[icart], sys, &pbc);
            f += df;
            if dfd > fdist_local {
                fdist_local = dfd;
            }
            for &ncell in &neighbors {
                let (df, dfd) = fparc(icart, sys.latomfirst[ncell], sys, &pbc);
                f += df;
                if dfd > fdist_local {
                    fdist_local = dfd;
                }
            }

            icart_id = sys.latomnext[icart];
        }

        icell_id = sys.lcellnext[icell];
    }

    if fdist_local > sys.fdist {
        sys.fdist = fdist_local;
    }
    f
}

#[cfg(feature = "rayon")]
fn accumulate_pair_f_parallel(sys: &PackContext) -> (F, F) {
    let pbc = pbc_constants(sys);
    sys.active_cells
        .par_iter()
        .map(|&icell| {
            let mut f: F = 0.0;
            let mut fdist_max: F = 0.0;
            let neighbors = sys.neighbor_cells_f[icell];

            let mut icart_id = sys.latomfirst[icell];
            while icart_id != NONE_IDX {
                let icart = icart_id as usize;
                let (f_same, fdist_same) = fparc_stats(icart, sys.latomnext[icart], sys, &pbc);
                f += f_same;
                fdist_max = fdist_max.max(fdist_same);

                for &ncell in &neighbors {
                    let (f_neigh, fdist_neigh) =
                        fparc_stats(icart, sys.latomfirst[ncell], sys, &pbc);
                    f += f_neigh;
                    fdist_max = fdist_max.max(fdist_neigh);
                }

                icart_id = sys.latomnext[icart];
            }

            (f, fdist_max)
        })
        .reduce(
            || (0.0 as F, 0.0 as F),
            |(f1, d1), (f2, d2)| (f1 + f2, d1.max(d2)),
        )
}

#[inline(always)]
fn accumulate_pair_g(sys: &mut PackContext) {
    let pbc = pbc_constants(sys);
    let mut icell_id = sys.lcellfirst;
    while icell_id != NONE_IDX {
        let icell = icell_id as usize;
        let neighbors = sys.neighbor_cells_g[icell];

        let mut icart_id = sys.latomfirst[icell];
        while icart_id != NONE_IDX {
            let icart = icart_id as usize;
            gparc(icart, sys.latomnext[icart], sys, &pbc);
            for &ncell in &neighbors {
                gparc(icart, sys.latomfirst[ncell], sys, &pbc);
            }

            icart_id = sys.latomnext[icart];
        }

        icell_id = sys.lcellnext[icell];
    }
}

#[inline(always)]
fn accumulate_pair_fg(sys: &mut PackContext) -> F {
    #[cfg(feature = "rayon")]
    if sys.parallel_pair_eval && !sys.move_flag {
        let (f, fdist_max) = accumulate_pair_fg_parallel(sys);
        if fdist_max > sys.fdist {
            sys.fdist = fdist_max;
        }
        return f;
    }

    let pbc = pbc_constants(sys);
    let mut f = 0.0;
    let mut fdist_local: F = 0.0;
    let mut icell_id = sys.lcellfirst;
    while icell_id != NONE_IDX {
        let icell = icell_id as usize;
        let neighbors = sys.neighbor_cells_g[icell];

        let mut icart_id = sys.latomfirst[icell];
        while icart_id != NONE_IDX {
            let icart = icart_id as usize;
            let (df, dfd) = fgparc(icart, sys.latomnext[icart], sys, &pbc);
            f += df;
            if dfd > fdist_local {
                fdist_local = dfd;
            }
            for &ncell in &neighbors {
                let (df, dfd) = fgparc(icart, sys.latomfirst[ncell], sys, &pbc);
                f += df;
                if dfd > fdist_local {
                    fdist_local = dfd;
                }
            }

            icart_id = sys.latomnext[icart];
        }

        icell_id = sys.lcellnext[icell];
    }

    if fdist_local > sys.fdist {
        sys.fdist = fdist_local;
    }
    f
}

/// Base pointer into the per-worker scratch gradient buffer
/// ([`crate::context::WorkBuffers::grad_partials`]). `*mut [F; 3]` is not
/// `Send`/`Sync`; this wrapper asserts each worker only ever dereferences its
/// own `[t*ntotat .. (t+1)*ntotat)` region, which is disjoint across the
/// concurrently-running tasks (keyed by the unique rayon pool thread index).
#[cfg(feature = "rayon")]
#[derive(Clone, Copy)]
struct PartialPtr {
    base: *mut [F; 3],
    ntotat: usize,
}
// SAFETY: see type doc — regions are keyed by the unique pool thread index, so
// no two concurrently-running tasks alias the same slot.
#[cfg(feature = "rayon")]
unsafe impl Send for PartialPtr {}
#[cfg(feature = "rayon")]
unsafe impl Sync for PartialPtr {}
#[cfg(feature = "rayon")]
impl PartialPtr {
    /// # Safety
    /// `t` must be the calling worker's unique pool index and `i < ntotat`.
    #[inline(always)]
    unsafe fn slot<'a>(self, t: usize, i: usize) -> &'a mut [F; 3] {
        unsafe { &mut *self.base.add(t * self.ntotat + i) }
    }
}

/// Parallel counterpart to [`accumulate_pair_fg`]. rayon work-steals over
/// `active_cells` using the **same 13-neighbor half-stencil the serial path
/// walks** ([`PackContext::neighbor_cells_g`]), so each unordered pair is
/// visited exactly once — none of the ~2× redundant distance work an
/// atom-centric full-stencil pass incurs.
///
/// Race freedom *without* that redundancy comes from per-worker scratch buffers
/// ([`crate::context::WorkBuffers::grad_partials`]): worker `t` accumulates
/// every `gi += d` / `gj -= d` half-stencil write into its own private region
/// `[t*ntotat .. (t+1)*ntotat)`, so no two concurrently-running tasks touch the
/// same slot even though a half-stencil writes into neighbor cells. After the
/// sweep the regions are reduced into `sys.work.gxcar` (which already holds the
/// constraint gradient from `expand_molecules`). The zeroing and the reduction
/// are themselves parallel and cost O(ntotat) wall-clock when `nthreads ≈
/// cores`.
///
/// The contribution set is identical to the serial path (same pairs, same
/// signs); only the summation order differs, within the tolerance the
/// `parallel_equivalence` tests pin. Gated on `!move_flag` by the caller (the
/// per-atom `fdist_atom` bookkeeping is intentionally skipped).
#[cfg(feature = "rayon")]
fn accumulate_pair_fg_parallel(sys: &mut PackContext) -> (F, F) {
    let ntotat = sys.ntotat;
    let nthreads = rayon::current_num_threads().max(1);

    // Per-worker scratch: take it out of `sys` so the context can be shared
    // immutably across tasks while each worker writes its own region through a
    // raw pointer. Sized to `nthreads * ntotat` and zeroed in parallel (each
    // worker clears its own region).
    let mut partials = std::mem::take(&mut sys.work.grad_partials);
    if partials.len() != nthreads * ntotat {
        partials.resize(nthreads * ntotat, [0.0; 3]);
    }
    partials
        .par_chunks_mut(ntotat.max(1))
        .for_each(|region| region.fill([0.0; 3]));

    let pbc = pbc_constants(sys);
    let sys_ro: &PackContext = sys;
    let pptr = PartialPtr {
        base: partials.as_mut_ptr(),
        ntotat,
    };

    let (f_total, fdist_max) = sys_ro
        .active_cells
        .par_iter()
        .map(|&icell| {
            // The worker's private region index. Inside a rayon parallel
            // closure this is always `Some(0..nthreads)`.
            let t = rayon::current_thread_index().unwrap_or(0);
            let neighbors = &sys_ro.neighbor_cells_g[icell];
            let mut f_local: F = 0.0;
            let mut fdist_local: F = 0.0;
            let mut icart_id = sys_ro.latomfirst[icell];
            while icart_id != NONE_IDX {
                let icart = icart_id as usize;
                // Same cell: forward atoms only (latomnext), each pair once.
                let (df, dfd) = fgparc_into(icart, sys_ro.latomnext[icart], sys_ro, pptr, t, &pbc);
                f_local += df;
                if dfd > fdist_local {
                    fdist_local = dfd;
                }
                // The 13 forward neighbor cells.
                for &ncell in neighbors {
                    let (df, dfd) =
                        fgparc_into(icart, sys_ro.latomfirst[ncell], sys_ro, pptr, t, &pbc);
                    f_local += df;
                    if dfd > fdist_local {
                        fdist_local = dfd;
                    }
                }
                icart_id = sys_ro.latomnext[icart];
            }
            (f_local, fdist_local)
        })
        .reduce(|| (0.0 as F, 0.0 as F), |a, b| (a.0 + b.0, a.1.max(b.1)));

    // Reduce the per-worker regions into the constraint gradient already in
    // gxcar. Parallel over atom chunks; within a chunk each worker region is
    // read as a contiguous slice (cache-friendly).
    let mut grad = std::mem::take(&mut sys.work.gxcar);
    let chunk = ntotat.div_ceil(nthreads).max(1);
    {
        let partials_ref: &[[F; 3]] = &partials;
        grad.par_chunks_mut(chunk)
            .enumerate()
            .for_each(|(ci, gchunk)| {
                let start = ci * chunk;
                let len = gchunk.len();
                for t in 0..nthreads {
                    let base = t * ntotat + start;
                    let region = &partials_ref[base..base + len];
                    for (g, p) in gchunk.iter_mut().zip(region) {
                        g[0] += p[0];
                        g[1] += p[1];
                        g[2] += p[2];
                    }
                }
            });
    }

    sys.work.gxcar = grad;
    sys.work.grad_partials = partials;
    (f_total, fdist_max)
}

/// [`fgparc`] for the parallel half-stencil path: identical pair math and the
/// same `gi += d` / `gj -= d` writes (each unordered pair once), but the writes
/// land in worker `t`'s private scratch region via `pptr` instead of
/// `sys.work.gxcar`, and the `move_flag` per-atom bookkeeping is omitted (the
/// caller gates on `!move_flag`). Reads a shared `&PackContext`; returns
/// `(penalty_sum, fdist_max)`.
///
/// `gi` (slot `icart`) and `gj` (slot `jcart`) are always distinct atoms — the
/// same-molecule skip rules out `icart == jcart` — and each `&mut` is scoped so
/// the two never coexist, so the disjoint raw writes carry no aliasing hazard.
#[cfg(feature = "rayon")]
#[inline(always)]
fn fgparc_into(
    icart: usize,
    first_jcart: u32,
    sys: &PackContext,
    pptr: PartialPtr,
    t: usize,
    pbc: &PbcConstants,
) -> (F, F) {
    let mut result: F = 0.0;
    let mut local_fdist: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        if let Some(c) = pair_term::<true, true>(&hot, xi, jcart, sys, pbc) {
            result += c.energy;
            if c.overlap {
                // SAFETY: `t` is this worker's unique region; `icart`/`jcart`
                // are distinct in-bounds atoms; each `&mut` is dropped before
                // the next.
                {
                    let gi = unsafe { pptr.slot(t, icart) };
                    gi[0] += c.grad[0];
                    gi[1] += c.grad[1];
                    gi[2] += c.grad[2];
                }
                {
                    let gj = unsafe { pptr.slot(t, jcart) };
                    gj[0] -= c.grad[0];
                    gj[1] -= c.grad[1];
                    gj[2] -= c.grad[2];
                }
            }
            if c.violation > local_fdist {
                local_fdist = c.violation;
            }
        }
        jcart_id = next;
    }

    (result, local_fdist)
}

/// Atom-pair gradient accumulation into `sys.work.gxcar`.
/// Port of `gparc.f90`.
#[inline(always)]
fn gparc(icart: usize, first_jcart: u32, sys: &mut PackContext, pbc: &PbcConstants) {
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        if let Some(c) = pair_term::<true, false>(&hot, xi, jcart, sys, pbc) {
            if c.overlap {
                let gi = &mut sys.work.gxcar[icart];
                gi[0] += c.grad[0];
                gi[1] += c.grad[1];
                gi[2] += c.grad[2];
                let gj = &mut sys.work.gxcar[jcart];
                gj[0] -= c.grad[0];
                gj[1] -= c.grad[1];
                gj[2] -= c.grad[2];
            }
        }
        jcart_id = next;
    }
}

/// Atom-pair function and gradient accumulation into `sys.work.gxcar`.
///
/// Returns `(penalty_sum, fdist_max)`. Caller reduces `fdist_max` locally
/// and writes `sys.fdist` once after the cell walk completes.
#[inline(always)]
fn fgparc(icart: usize, first_jcart: u32, sys: &mut PackContext, pbc: &PbcConstants) -> (F, F) {
    let mut result = 0.0;
    let mut local_fdist: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);
    let move_flag = sys.move_flag;

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        if let Some(c) = pair_term::<true, true>(&hot, xi, jcart, sys, pbc) {
            result += c.energy;
            if c.overlap {
                let gi = &mut sys.work.gxcar[icart];
                gi[0] += c.grad[0];
                gi[1] += c.grad[1];
                gi[2] += c.grad[2];
                let gj = &mut sys.work.gxcar[jcart];
                gj[0] -= c.grad[0];
                gj[1] -= c.grad[1];
                gj[2] -= c.grad[2];
            }
            if c.violation > local_fdist {
                local_fdist = c.violation;
            }
            if move_flag {
                if c.violation > sys.fdist_atom[icart] {
                    sys.fdist_atom[icart] = c.violation;
                }
                if c.violation > sys.fdist_atom[jcart] {
                    sys.fdist_atom[jcart] = c.violation;
                }
            }
        }
        jcart_id = next;
    }

    (result, local_fdist)
}

#[cfg(feature = "rayon")]
#[inline(always)]
fn fparc_stats(icart: usize, first_jcart: u32, sys: &PackContext, pbc: &PbcConstants) -> (F, F) {
    let mut result: F = 0.0;
    let mut fdist_max: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        if let Some(c) = pair_term::<false, true>(&hot, xi, jcart, sys, pbc) {
            result += c.energy;
            if c.violation > fdist_max {
                fdist_max = c.violation;
            }
        }
        jcart_id = next;
    }

    (result, fdist_max)
}

/// Project each active molecule's atom Cartesian gradient onto its own 6 DOF
/// (COM + Euler), writing only its own **disjoint** `g[ilubar..]` / `g[ilugan..]`
/// slots. Because no slot is summed across molecules and each molecule's
/// accumulation order matches the legacy serial loop's (same atoms, same axis
/// order), the result is **bit-identical** whether Phase A iterates with `iter`
/// or `par_iter` — the only difference between the serial and parallel paths
/// (see [`project_for_each`]). The dominant cost is `eulerrmat_derivatives`
/// (trig) per molecule.
fn project_cartesian_gradient(x: &[F], sys: &mut PackContext, g: &mut [F]) {
    g.iter_mut().for_each(|v| *v = 0.0);

    let mut descs = std::mem::take(&mut sys.work.mol_descs);
    fill_active_mol_descs(&mut descs, sys);
    let move_flag = sys.move_flag;
    let parallel = sys.parallel_pair_eval && !move_flag;

    #[derive(Clone, Copy)]
    struct GPtr(*mut F);
    // SAFETY: per-molecule `ilubar`/`ilugan` ranges are disjoint, so no two tasks
    // write the same `g` index.
    unsafe impl Send for GPtr {}
    unsafe impl Sync for GPtr {}
    impl GPtr {
        // Takes `self` so a closure that calls it captures the whole (`Sync`)
        // wrapper, not the bare `*mut` field.
        /// # Safety
        /// `i` must index a slot owned solely by the calling task.
        #[inline(always)]
        unsafe fn write(self, i: usize, v: F) {
            unsafe {
                *self.0.add(i) = v;
            }
        }
    }
    let gptr = GPtr(g.as_mut_ptr());

    // Inner scope so the shared `&PackContext` borrow ends before `sys` is reused.
    {
        let sys_ro: &PackContext = sys;
        let body = |&(itype, icart0, ilubar, ilugan): &(usize, usize, usize, usize)| {
            let beta = x[ilugan];
            let gama = x[ilugan + 1];
            let teta = x[ilugan + 2];
            let (dv1beta, dv1gama, dv1teta, dv2beta, dv2gama, dv2teta, dv3beta, dv3gama, dv3teta) =
                eulerrmat_derivatives(beta, gama, teta);

            let idbase = sys_ro.idfirst[itype];
            let na = sys_ro.natoms[itype];
            let mut gcom = [0.0 as F; 3];
            let mut gang = [0.0 as F; 3];
            for iatom in 0..na {
                let gx = sys_ro.work.gxcar[icart0 + iatom];
                let cr = sys_ro.coor[idbase + iatom];
                for k in 0..3 {
                    gcom[k] += gx[k];
                }
                for k in 0..3 {
                    gang[0] +=
                        (cr[0] * dv1beta[k] + cr[1] * dv2beta[k] + cr[2] * dv3beta[k]) * gx[k];
                    gang[1] +=
                        (cr[0] * dv1gama[k] + cr[1] * dv2gama[k] + cr[2] * dv3gama[k]) * gx[k];
                    gang[2] +=
                        (cr[0] * dv1teta[k] + cr[1] * dv2teta[k] + cr[2] * dv3teta[k]) * gx[k];
                }
            }
            // SAFETY: `ilubar`/`ilugan` index this molecule's own disjoint slots.
            unsafe {
                for k in 0..3 {
                    gptr.write(ilubar + k, gcom[k]);
                    gptr.write(ilugan + k, gang[k]);
                }
            }
        };
        project_for_each(&descs, &body, parallel);
    }

    sys.work.mol_descs = descs;
}

/// Run the independent per-molecule projection `body` over `descs`. The single
/// place the serial and parallel projection paths diverge: `par_iter` when
/// `parallel`, plain `iter` otherwise. No reduction — each task writes disjoint
/// `g` slots.
#[cfg(feature = "rayon")]
fn project_for_each<B>(descs: &[(usize, usize, usize, usize)], body: &B, parallel: bool)
where
    B: Fn(&(usize, usize, usize, usize)) + Sync,
{
    use rayon::prelude::*;
    if parallel {
        descs.par_iter().for_each(body);
    } else {
        descs.iter().for_each(body);
    }
}

/// Non-rayon build: projection is always serial.
#[cfg(not(feature = "rayon"))]
fn project_for_each<B>(descs: &[(usize, usize, usize, usize)], body: &B, _parallel: bool)
where
    B: Fn(&(usize, usize, usize, usize)),
{
    descs.iter().for_each(body);
}

#[inline]
fn matches_cached_geometry(x: &[F], sys: &PackContext) -> bool {
    sys.work.matches_cached_geometry(
        x,
        &sys.comptype,
        sys.init1,
        sys.ncells,
        sys.cell_length,
        sys.pbc_min,
        sys.pbc_length,
        sys.pbc_periodic,
    )
}

#[inline]
fn update_cached_geometry(x: &[F], sys: &mut PackContext) {
    sys.work.update_cached_geometry(
        x,
        &sys.comptype,
        sys.init1,
        sys.ncells,
        sys.cell_length,
        sys.pbc_min,
        sys.pbc_length,
        sys.pbc_periodic,
    );
}

// ── Phase A.5 — Objective trait ────────────────────────────────────────────
//
// `Objective` is the abstraction the packer's GENCAN loop will talk to. At
// this checkpoint the trait is defined and implemented for `PackContext` but
// no call site has been rewired yet (`pgencan` still takes `&mut PackContext`
// directly). Phase A.6 swaps `pgencan`'s signature to `&mut dyn Objective`
// and is gated by an explicit extract-bench-loop commit so the dyn-dispatch
// cost lands with a measurement attached.
//
// The trait is intentionally shaped to match what the GENCAN loop reads and
// writes today — no speculative extra methods.

/// Abstracts the packer's objective function so the optimizer can talk to
/// any `(f, g)` oracle, not just `PackContext`.
///
/// Implementors are responsible for:
/// - Returning `f`, worst-atom distance violation (`fdist`), worst-molecule
///   restraint violation (`frest`) in one pass, via [`evaluate`].
/// - Exposing the cumulative function / gradient counters that the caller
///   resets between phases and reads for logging.
///
/// The trait does **not** own bounds (`l`, `u`): those are problem-specific
/// and the caller (e.g. `pgencan::build_bounds`) builds them from the
/// concrete context it has in hand.
///
/// See spec `.claude/specs/molrs-pack-plugin-arch.md` §9 Phase A step 5 for
/// the rationale; Phase B will expose this trait publicly as the extension
/// hook for custom objectives.
///
/// [`evaluate`]: Self::evaluate
pub trait Objective {
    /// Unified evaluation entry point. `mode` selects between `f` only,
    /// gradient only, and both. When a gradient is requested, `gradient`
    /// must be `Some(buf)` with `buf.len() == x.len()`; when it is not,
    /// `gradient` is ignored.
    ///
    /// On return, the implementor's internal `fdist` / `frest` state
    /// reflects this call (so a subsequent `self.fdist()` / `self.frest()`
    /// returns the same values as the `EvalOutput`'s `fdist_max` /
    /// `frest_max`).
    fn evaluate(&mut self, x: &[F], mode: EvalMode, gradient: Option<&mut [F]>) -> EvalOutput;

    /// Worst-atom distance violation from the most recent `evaluate` call.
    fn fdist(&self) -> F;

    /// Worst-molecule restraint violation from the most recent `evaluate`
    /// call.
    fn frest(&self) -> F;

    /// Cumulative count of function-value evaluations since the last
    /// `reset_eval_counters` call.
    fn ncf(&self) -> usize;

    /// Cumulative count of gradient evaluations since the last
    /// `reset_eval_counters` call.
    fn ncg(&self) -> usize;

    /// Zero the function / gradient counters. GENCAN calls this at the
    /// start of each outer iteration.
    fn reset_eval_counters(&mut self);

    /// Fill `l` and `u` (each of length `x.len()` at the `pgencan` entry)
    /// with per-variable bounds. The default implementation sets every
    /// variable to `[-1e20, +1e20]` (effectively unbounded); `PackContext`
    /// overrides it to add the Euler-angle bounds implied by
    /// `constrain_rotation`.
    ///
    /// Landed in A.6 so `pgencan` no longer needs `&mut PackContext` for
    /// anything beyond evaluation — bounds construction is now behind the
    /// trait too.
    fn bounds(&self, l: &mut [F], u: &mut [F]) {
        debug_assert_eq!(l.len(), u.len(), "bounds: l/u length mismatch");
        l.fill(-1.0e20);
        u.fill(1.0e20);
    }
}

impl Objective for PackContext {
    #[inline]
    fn evaluate(&mut self, x: &[F], mode: EvalMode, gradient: Option<&mut [F]>) -> EvalOutput {
        PackContext::evaluate(self, x, mode, gradient)
    }

    #[inline]
    fn fdist(&self) -> F {
        self.fdist
    }

    #[inline]
    fn frest(&self) -> F {
        self.frest
    }

    #[inline]
    fn ncf(&self) -> usize {
        PackContext::ncf(self)
    }

    #[inline]
    fn ncg(&self) -> usize {
        PackContext::ncg(self)
    }

    #[inline]
    fn reset_eval_counters(&mut self) {
        PackContext::reset_eval_counters(self);
    }

    /// Port of `gencan::build_bounds` (pre-A.6). COM variables (first `n/2`)
    /// stay `[-1e20, +1e20]`; Euler variables (last `n/2`) inherit
    /// `rot_bound` when `constrain_rot` is set for the owning type /
    /// molecule / axis.
    fn bounds(&self, l: &mut [F], u: &mut [F]) {
        debug_assert_eq!(l.len(), u.len(), "bounds: l/u length mismatch");
        let n = l.len();
        l.fill(-1.0e20);
        u.fill(1.0e20);
        let mut i = n / 2;
        for itype in 0..self.ntype {
            if !self.comptype[itype] {
                continue;
            }
            for _imol in 0..self.nmols[itype] {
                for axis in 0..3 {
                    if self.constrain_rot[itype][axis] {
                        let center = self.rot_bound[itype][axis][0];
                        let half_width = self.rot_bound[itype][axis][1].abs();
                        l[i] = center - half_width;
                        u[i] = center + half_width;
                    }
                    i += 1;
                }
            }
        }
        debug_assert_eq!(i, n);
    }
}

#[cfg(test)]
mod objective_trait_tests {
    use super::*;
    use crate::PackContext;

    /// Pins `<PackContext as Objective>::evaluate` against the inherent
    /// `PackContext::evaluate`: calling through a `&mut dyn Objective` must
    /// return byte-identical `EvalOutput` and leave identical `fdist` /
    /// `frest` state when fed the same empty-molecule context.
    ///
    /// With `ntotmol=0`, `x` is empty; `evaluate(FOnly)` runs the full
    /// constraints-container dispatch on empty state and returns
    /// `(0, 0, 0)` deterministically. The test still exercises the trait
    /// dispatch path so a future behavioural divergence (e.g. the impl
    /// forwarding to the wrong method) would fail loudly.
    #[test]
    fn dyn_objective_matches_inherent_evaluate() {
        fn build(ntotat: usize) -> PackContext {
            let mut sys = PackContext::new(ntotat, 0, 0);
            sys.radius.fill(0.75);
            sys.radius_ini.fill(1.5);
            sys.work.radiuswork.resize(ntotat, 0.0);
            sys.sync_atom_props();
            sys
        }
        let x: Vec<F> = Vec::new();

        let mut via_inherent = build(4);
        let out_inherent = PackContext::evaluate(&mut via_inherent, &x, EvalMode::FOnly, None);

        let mut via_trait_owner = build(4);
        let out_trait = {
            let obj: &mut dyn Objective = &mut via_trait_owner;
            obj.evaluate(&x, EvalMode::FOnly, None)
        };

        assert_eq!(out_inherent.f_total, out_trait.f_total, "f_total mismatch");
        assert_eq!(
            out_inherent.fdist_max, out_trait.fdist_max,
            "fdist_max mismatch"
        );
        assert_eq!(
            out_inherent.frest_max, out_trait.frest_max,
            "frest_max mismatch"
        );
        assert_eq!(
            via_inherent.fdist, via_trait_owner.fdist,
            "post-call fdist drift"
        );
        assert_eq!(
            via_inherent.frest, via_trait_owner.frest,
            "post-call frest drift"
        );

        // Counter + reset contract
        let obj: &mut dyn Objective = &mut via_trait_owner;
        assert_eq!(obj.ncf(), 1, "ncf should be 1 after one FOnly evaluate");
        assert_eq!(obj.ncg(), 0, "ncg should stay 0 under FOnly");
        obj.reset_eval_counters();
        assert_eq!(obj.ncf(), 0, "ncf should be 0 after reset");
        assert_eq!(obj.ncg(), 0, "ncg should be 0 after reset");
    }
}
