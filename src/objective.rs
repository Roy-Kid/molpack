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
        if length[k] > 0.0 {
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
/// inside `fparc` / `gparc` / `fgparc` (and the `*_stats` rayon variants).
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
pub fn compute_f(x: &[F], sys: &mut PackContext) -> F {
    sys.debug_assert_atom_props_sync();
    sys.increment_ncf();
    sys.fdist = 0.0;
    sys.frest = 0.0;

    if !sys.init1 {
        sys.resetcells();
    }

    let mut f = expand_molecules(x, sys, ExpandMode::F);

    if sys.init1 {
        update_cached_geometry(x, sys);
        return f;
    }

    f += accumulate_pair_f(sys);
    update_cached_geometry(x, sys);

    f
}

/// Compute objective function value and gradient in one pass over geometry/state.
/// This avoids rebuilding Cartesian coordinates and cell lists twice.
pub fn compute_fg(x: &[F], sys: &mut PackContext, g: &mut [F]) -> F {
    sys.debug_assert_atom_props_sync();
    sys.increment_ncf();
    sys.increment_ncg();
    sys.fdist = 0.0;
    sys.frest = 0.0;
    sys.work.gxcar.fill([0.0; 3]);

    if !sys.init1 {
        sys.resetcells();
    }

    let mut f = expand_molecules(x, sys, ExpandMode::FG);

    if !sys.init1 {
        f += accumulate_pair_fg(sys);
    }

    update_cached_geometry(x, sys);
    project_cartesian_gradient(x, sys, g);

    f
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
fn fparc(
    icart: usize,
    first_jcart: u32,
    sys: &mut PackContext,
    pbc: &PbcConstants,
) -> (F, F) {
    let mut result = 0.0;
    let mut local_fdist: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);
    let move_flag = sys.move_flag;

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        let props_j = sys.atom_props[jcart];
        // Skip same molecule
        if hot.props.ibmol == props_j.ibmol && hot.props.ibtype == props_j.ibtype {
            jcart_id = next;
            continue;
        }
        // Skip two fixed atoms (only if the system has any at all).
        if hot.fixed_i && (props_j.flags & ATOM_FLAG_FIXED != 0) {
            jcart_id = next;
            continue;
        }

        let xj = sys.xcart[jcart];
        let (dx, dy, dz) = pbc_wrap_delta(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2], pbc);
        let datom = dx * dx + dy * dy + dz * dz;
        let rsum = hot.props.radius + props_j.radius;
        let tol = rsum * rsum;

        if datom < tol {
            let penalty = datom - tol;
            let scale = hot.props.fscale * props_j.fscale;
            result += scale * penalty * penalty;
            if hot.has_short && (hot.use_short_i || (props_j.flags & ATOM_FLAG_SHORT != 0)) {
                let short_rsum = hot.shrad_i + sys.short_radius[jcart];
                let short_tol = short_rsum * short_rsum;
                if datom < short_tol {
                    let penalty = datom - short_tol;
                    let mut sr_scale = (hot.shscl_i * sys.short_radius_scale[jcart]).sqrt();
                    sr_scale *= (tol * tol) / (short_tol * short_tol);
                    result += scale * sr_scale * penalty * penalty;
                }
            }
        }

        let rsum_ini = hot.props.radius_ini + props_j.radius_ini;
        let tol_ini = rsum_ini * rsum_ini;
        let violation = tol_ini - datom;
        if violation > local_fdist {
            local_fdist = violation;
        }
        if move_flag {
            if violation > sys.fdist_atom[icart] {
                sys.fdist_atom[icart] = violation;
            }
            if violation > sys.fdist_atom[jcart] {
                sys.fdist_atom[jcart] = violation;
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
    }

    update_cached_geometry(x, sys);
    project_cartesian_gradient(x, sys, g);
}

/// Expand all active molecules into `sys.xcart`, applying either constraint values
/// (`computef`) or constraint gradients (`computeg`), and rebuilding the cell list.
#[inline]
fn expand_molecules(x: &[F], sys: &mut PackContext, mode: ExpandMode) -> F {
    let mut f = 0.0;
    let mut ilubar = 0usize;
    let mut ilugan = sys.ntotmol * 3;
    let mut icart = 0usize;

    // Packmol computef/computeg loop only over free types (1..ntype).
    for itype in 0..sys.ntype {
        if !sys.comptype[itype] {
            icart += sys.nmols[itype] * sys.natoms[itype];
            continue;
        }

        for _imol in 0..sys.nmols[itype] {
            let xcm = [x[ilubar], x[ilubar + 1], x[ilubar + 2]];
            let beta = x[ilugan];
            let gama = x[ilugan + 1];
            let teta = x[ilugan + 2];
            let (v1, v2, v3) = eulerrmat(beta, gama, teta);
            let idatom_base = sys.idfirst[itype];

            for iatom in 0..sys.natoms[itype] {
                let idatom = idatom_base + iatom;
                let pos = compcart(&xcm, &sys.coor[idatom], &v1, &v2, &v3);
                sys.xcart[icart] = pos;

                match mode {
                    ExpandMode::F => {
                        f += accumulate_constraint_value(icart, &pos, sys);
                    }
                    ExpandMode::G => {
                        accumulate_constraint_gradient(icart, &pos, sys);
                    }
                    ExpandMode::FG => {
                        f += accumulate_constraint_value(icart, &pos, sys);
                        accumulate_constraint_gradient(icart, &pos, sys);
                    }
                }

                if !sys.init1 {
                    insert_atom_in_cell(icart, &pos, sys);
                }

                icart += 1;
            }

            ilugan += 3;
            ilubar += 3;
        }
    }

    f
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

#[inline(always)]
fn insert_atom_in_cell(icart: usize, pos: &[F; 3], sys: &mut PackContext) {
    let cell = setcell(
        pos,
        &sys.pbc_min,
        &sys.pbc_length,
        &sys.cell_length,
        &sys.ncells,
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

/// Parallel counterpart to [`accumulate_pair_fg`]. Partitions `active_cells`
/// into `N_threads` contiguous chunks, runs each on its own partial
/// gradient buffer, then serially merges the buffers back into
/// `sys.work.gxcar`. Gated on `!move_flag` (the per-atom `fdist_atom`
/// bookkeeping is intentionally skipped).
#[cfg(feature = "rayon")]
fn accumulate_pair_fg_parallel(sys: &mut PackContext) -> (F, F) {
    let ntotat = sys.ntotat;
    let n_threads = rayon::current_num_threads().max(1);
    sys.work.reset_partial_gxcar(n_threads, ntotat);

    // Move partial buffers out so `sys` can be borrowed immutably across threads.
    let mut partials = std::mem::take(&mut sys.work.partial_gxcar);
    let pbc = pbc_constants(sys);
    let sys_ro: &PackContext = sys;
    let active_cells = sys_ro.active_cells.as_slice();
    let chunk_size = active_cells.len().div_ceil(n_threads);

    let (f_total, fdist_max) = partials
        .par_iter_mut()
        .enumerate()
        .map(|(tidx, grad)| {
            let start = tidx * chunk_size;
            let end = ((tidx + 1) * chunk_size).min(active_cells.len());
            let mut f_local: F = 0.0;
            let mut fdist_local: F = 0.0;
            for &icell in &active_cells[start..end] {
                let neighbors = sys_ro.neighbor_cells_g[icell];
                let mut icart_id = sys_ro.latomfirst[icell];
                while icart_id != NONE_IDX {
                    let icart = icart_id as usize;
                    let (df, dfd) =
                        fgparc_stats(icart, sys_ro.latomnext[icart], sys_ro, grad, &pbc);
                    f_local += df;
                    if dfd > fdist_local {
                        fdist_local = dfd;
                    }
                    for &ncell in &neighbors {
                        let (df, dfd) =
                            fgparc_stats(icart, sys_ro.latomfirst[ncell], sys_ro, grad, &pbc);
                        f_local += df;
                        if dfd > fdist_local {
                            fdist_local = dfd;
                        }
                    }
                    icart_id = sys_ro.latomnext[icart];
                }
            }
            (f_local, fdist_local)
        })
        .reduce(
            || (0.0 as F, 0.0 as F),
            |a, b| (a.0 + b.0, a.1.max(b.1)),
        );

    // Merge per-thread gradients back into the main buffer. Small
    // systems stay serial to avoid rayon dispatch overhead dominating
    // the merge itself — the threshold is calibrated via the
    // `partial_gradient_merge` bench.
    if ntotat >= crate::context::work_buffers::MERGE_PARALLEL_THRESHOLD_NTOTAT {
        crate::context::WorkBuffers::merge_partials_parallel(&partials, &mut sys.work.gxcar);
    } else {
        crate::context::WorkBuffers::merge_partials_serial(&partials, &mut sys.work.gxcar);
    }

    sys.work.partial_gxcar = partials;
    (f_total, fdist_max)
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
        let props_j = sys.atom_props[jcart];
        if hot.props.ibmol == props_j.ibmol && hot.props.ibtype == props_j.ibtype {
            jcart_id = next;
            continue;
        }
        if hot.fixed_i && (props_j.flags & ATOM_FLAG_FIXED != 0) {
            jcart_id = next;
            continue;
        }

        let rsum = hot.props.radius + props_j.radius;
        let tol = rsum * rsum;
        let xj = sys.xcart[jcart];
        let (dx, dy, dz) = pbc_wrap_delta(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2], pbc);
        let datom = dx * dx + dy * dy + dz * dz;

        if datom < tol {
            let scale = hot.props.fscale * props_j.fscale;
            let dtemp = scale * 4.0 * (datom - tol);
            let xdiff0 = dtemp * dx;
            let xdiff1 = dtemp * dy;
            let xdiff2 = dtemp * dz;
            let gi = &mut sys.work.gxcar[icart];
            gi[0] += xdiff0;
            gi[1] += xdiff1;
            gi[2] += xdiff2;
            let gj = &mut sys.work.gxcar[jcart];
            gj[0] -= xdiff0;
            gj[1] -= xdiff1;
            gj[2] -= xdiff2;
            if hot.has_short && (hot.use_short_i || (props_j.flags & ATOM_FLAG_SHORT != 0)) {
                let short_rsum = hot.shrad_i + sys.short_radius[jcart];
                let short_tol = short_rsum * short_rsum;
                if datom < short_tol {
                    let mut sr_scale = (hot.shscl_i * sys.short_radius_scale[jcart]).sqrt();
                    sr_scale *= (tol * tol) / (short_tol * short_tol);
                    let dtemp2 = scale * 4.0 * sr_scale * (datom - short_tol);
                    let xdiff0 = dtemp2 * dx;
                    let xdiff1 = dtemp2 * dy;
                    let xdiff2 = dtemp2 * dz;
                    let gi = &mut sys.work.gxcar[icart];
                    gi[0] += xdiff0;
                    gi[1] += xdiff1;
                    gi[2] += xdiff2;
                    let gj = &mut sys.work.gxcar[jcart];
                    gj[0] -= xdiff0;
                    gj[1] -= xdiff1;
                    gj[2] -= xdiff2;
                }
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
fn fgparc(
    icart: usize,
    first_jcart: u32,
    sys: &mut PackContext,
    pbc: &PbcConstants,
) -> (F, F) {
    let mut result = 0.0;
    let mut local_fdist: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);
    let move_flag = sys.move_flag;

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        let props_j = sys.atom_props[jcart];
        if hot.props.ibmol == props_j.ibmol && hot.props.ibtype == props_j.ibtype {
            jcart_id = next;
            continue;
        }
        if hot.fixed_i && (props_j.flags & ATOM_FLAG_FIXED != 0) {
            jcart_id = next;
            continue;
        }

        let xj = sys.xcart[jcart];
        let (dx, dy, dz) = pbc_wrap_delta(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2], pbc);
        let datom = dx * dx + dy * dy + dz * dz;
        let rsum = hot.props.radius + props_j.radius;
        let tol = rsum * rsum;

        if datom < tol {
            let penalty = datom - tol;
            let scale = hot.props.fscale * props_j.fscale;
            result += scale * penalty * penalty;

            let dtemp = scale * 4.0 * penalty;
            let xdiff0 = dtemp * dx;
            let xdiff1 = dtemp * dy;
            let xdiff2 = dtemp * dz;
            let gi = &mut sys.work.gxcar[icart];
            gi[0] += xdiff0;
            gi[1] += xdiff1;
            gi[2] += xdiff2;
            let gj = &mut sys.work.gxcar[jcart];
            gj[0] -= xdiff0;
            gj[1] -= xdiff1;
            gj[2] -= xdiff2;

            if hot.has_short && (hot.use_short_i || (props_j.flags & ATOM_FLAG_SHORT != 0)) {
                let short_rsum = hot.shrad_i + sys.short_radius[jcart];
                let short_tol = short_rsum * short_rsum;
                if datom < short_tol {
                    let short_penalty = datom - short_tol;
                    let mut sr_scale = (hot.shscl_i * sys.short_radius_scale[jcart]).sqrt();
                    sr_scale *= (tol * tol) / (short_tol * short_tol);
                    let sr_pair_scale = scale * sr_scale;
                    result += sr_pair_scale * short_penalty * short_penalty;

                    let dtemp2 = sr_pair_scale * 4.0 * short_penalty;
                    let xdiff0 = dtemp2 * dx;
                    let xdiff1 = dtemp2 * dy;
                    let xdiff2 = dtemp2 * dz;
                    let gi = &mut sys.work.gxcar[icart];
                    gi[0] += xdiff0;
                    gi[1] += xdiff1;
                    gi[2] += xdiff2;
                    let gj = &mut sys.work.gxcar[jcart];
                    gj[0] -= xdiff0;
                    gj[1] -= xdiff1;
                    gj[2] -= xdiff2;
                }
            }
        }

        let rsum_ini = hot.props.radius_ini + props_j.radius_ini;
        let tol_ini = rsum_ini * rsum_ini;
        let violation = tol_ini - datom;
        if violation > local_fdist {
            local_fdist = violation;
        }
        if move_flag {
            if violation > sys.fdist_atom[icart] {
                sys.fdist_atom[icart] = violation;
            }
            if violation > sys.fdist_atom[jcart] {
                sys.fdist_atom[jcart] = violation;
            }
        }

        jcart_id = next;
    }

    (result, local_fdist)
}

/// Parallel-friendly variant of [`fgparc`] — writes the pair gradient into
/// an external `grad` buffer instead of `sys.work.gxcar`, so multiple rayon
/// threads can accumulate into disjoint partial buffers and the caller
/// merges them. `fdist_atom` bookkeeping is intentionally dropped: the
/// parallel fast path is gated on `!move_flag`.
#[cfg(feature = "rayon")]
#[inline(always)]
fn fgparc_stats(
    icart: usize,
    first_jcart: u32,
    sys: &PackContext,
    grad: &mut [[F; 3]],
    pbc: &PbcConstants,
) -> (F, F) {
    let mut result: F = 0.0;
    let mut fdist_max: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        let props_j = sys.atom_props[jcart];
        if hot.props.ibmol == props_j.ibmol && hot.props.ibtype == props_j.ibtype {
            jcart_id = next;
            continue;
        }
        if hot.fixed_i && (props_j.flags & ATOM_FLAG_FIXED != 0) {
            jcart_id = next;
            continue;
        }

        let xj = sys.xcart[jcart];
        let (dx, dy, dz) = pbc_wrap_delta(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2], pbc);
        let datom = dx * dx + dy * dy + dz * dz;
        let rsum = hot.props.radius + props_j.radius;
        let tol = rsum * rsum;

        if datom < tol {
            let penalty = datom - tol;
            let scale = hot.props.fscale * props_j.fscale;
            result += scale * penalty * penalty;

            let dtemp = scale * 4.0 * penalty;
            let xdiff0 = dtemp * dx;
            let xdiff1 = dtemp * dy;
            let xdiff2 = dtemp * dz;
            grad[icart][0] += xdiff0;
            grad[icart][1] += xdiff1;
            grad[icart][2] += xdiff2;
            grad[jcart][0] -= xdiff0;
            grad[jcart][1] -= xdiff1;
            grad[jcart][2] -= xdiff2;

            if hot.has_short && (hot.use_short_i || (props_j.flags & ATOM_FLAG_SHORT != 0)) {
                let short_rsum = hot.shrad_i + sys.short_radius[jcart];
                let short_tol = short_rsum * short_rsum;
                if datom < short_tol {
                    let short_penalty = datom - short_tol;
                    let mut sr_scale = (hot.shscl_i * sys.short_radius_scale[jcart]).sqrt();
                    sr_scale *= (tol * tol) / (short_tol * short_tol);
                    let sr_pair_scale = scale * sr_scale;
                    result += sr_pair_scale * short_penalty * short_penalty;

                    let dtemp2 = sr_pair_scale * 4.0 * short_penalty;
                    let xdiff0 = dtemp2 * dx;
                    let xdiff1 = dtemp2 * dy;
                    let xdiff2 = dtemp2 * dz;
                    grad[icart][0] += xdiff0;
                    grad[icart][1] += xdiff1;
                    grad[icart][2] += xdiff2;
                    grad[jcart][0] -= xdiff0;
                    grad[jcart][1] -= xdiff1;
                    grad[jcart][2] -= xdiff2;
                }
            }
        }

        let rsum_ini = hot.props.radius_ini + props_j.radius_ini;
        let tol_ini = rsum_ini * rsum_ini;
        let violation = tol_ini - datom;
        if violation > fdist_max {
            fdist_max = violation;
        }

        jcart_id = next;
    }

    (result, fdist_max)
}

#[cfg(feature = "rayon")]
#[inline(always)]
fn fparc_stats(
    icart: usize,
    first_jcart: u32,
    sys: &PackContext,
    pbc: &PbcConstants,
) -> (F, F) {
    let mut result: F = 0.0;
    let mut fdist_max: F = 0.0;
    let mut jcart_id = first_jcart;
    let xi = sys.xcart[icart];
    let hot = AtomHotState::load(icart, sys);

    while jcart_id != NONE_IDX {
        let jcart = jcart_id as usize;
        let next = sys.latomnext[jcart];
        let props_j = sys.atom_props[jcart];
        if hot.props.ibmol == props_j.ibmol && hot.props.ibtype == props_j.ibtype {
            jcart_id = next;
            continue;
        }
        if hot.fixed_i && (props_j.flags & ATOM_FLAG_FIXED != 0) {
            jcart_id = next;
            continue;
        }

        let xj = sys.xcart[jcart];
        let (dx, dy, dz) = pbc_wrap_delta(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2], pbc);
        let datom = dx * dx + dy * dy + dz * dz;
        let rsum = hot.props.radius + props_j.radius;
        let tol = rsum * rsum;

        if datom < tol {
            let penalty = datom - tol;
            let scale = hot.props.fscale * props_j.fscale;
            result += scale * penalty * penalty;
            if hot.has_short && (hot.use_short_i || (props_j.flags & ATOM_FLAG_SHORT != 0)) {
                let short_rsum = hot.shrad_i + sys.short_radius[jcart];
                let short_tol = short_rsum * short_rsum;
                if datom < short_tol {
                    let penalty = datom - short_tol;
                    let mut sr_scale = (hot.shscl_i * sys.short_radius_scale[jcart]).sqrt();
                    sr_scale *= (tol * tol) / (short_tol * short_tol);
                    result += scale * sr_scale * penalty * penalty;
                }
            }
        }

        let rsum_ini = hot.props.radius_ini + props_j.radius_ini;
        let tol_ini = rsum_ini * rsum_ini;
        let violation = tol_ini - datom;
        if violation > fdist_max {
            fdist_max = violation;
        }

        jcart_id = next;
    }

    (result, fdist_max)
}

fn project_cartesian_gradient(x: &[F], sys: &mut PackContext, g: &mut [F]) {
    g.iter_mut().for_each(|v| *v = 0.0);

    let mut k1 = 0usize;
    let mut k2 = sys.ntotmol * 3;
    let mut icart = 0usize;

    for itype in 0..sys.ntype {
        if !sys.comptype[itype] {
            icart += sys.nmols[itype] * sys.natoms[itype];
            continue;
        }

        for _imol in 0..sys.nmols[itype] {
            let beta = x[k2];
            let gama = x[k2 + 1];
            let teta = x[k2 + 2];

            let (dv1beta, dv1gama, dv1teta, dv2beta, dv2gama, dv2teta, dv3beta, dv3gama, dv3teta) =
                eulerrmat_derivatives(beta, gama, teta);

            let idatom_base = sys.idfirst[itype];
            for iatom in 0..sys.natoms[itype] {
                let idatom = idatom_base + iatom;
                let cr = sys.coor[idatom];

                for k in 0..3 {
                    g[k1 + k] += sys.work.gxcar[icart][k];
                }

                for k in 0..3 {
                    g[k2] += (cr[0] * dv1beta[k] + cr[1] * dv2beta[k] + cr[2] * dv3beta[k])
                        * sys.work.gxcar[icart][k];
                    g[k2 + 1] += (cr[0] * dv1gama[k] + cr[1] * dv2gama[k] + cr[2] * dv3gama[k])
                        * sys.work.gxcar[icart][k];
                    g[k2 + 2] += (cr[0] * dv1teta[k] + cr[1] * dv2teta[k] + cr[2] * dv3teta[k])
                        * sys.work.gxcar[icart][k];
                }

                icart += 1;
            }

            k2 += 3;
            k1 += 3;
        }
    }
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
