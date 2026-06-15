#![allow(clippy::needless_range_loop)]
//! Finite-difference gradient consistency tests.

use std::sync::Arc;

use molpack::objective::{compute_f, compute_fg, compute_g};
use molpack::{
    AboveGaussianRestraint, AbovePlaneRestraint, AtomRestraint, BelowGaussianRestraint,
    BelowPlaneRestraint, F, InsideBoxRestraint, InsideCubeRestraint, InsideCylinderRestraint,
    InsideEllipsoidRestraint, InsideSphereRestraint, OutsideBoxRestraint, OutsideCubeRestraint,
    OutsideCylinderRestraint, OutsideEllipsoidRestraint, OutsideSphereRestraint, PackContext,
};

// ── helpers ────────────────────────────────────────────────────────────────

/// Central finite-difference gradient for variable `i`.
fn finite_diff(x: &[F], sys: &mut PackContext, i: usize, h: F) -> F {
    let mut xp = x.to_vec();
    let mut xm = x.to_vec();
    xp[i] += h;
    xm[i] -= h;
    let fp = compute_f(&xp, sys);
    let fm = compute_f(&xm, sys);
    (fp - fm) / (2.0 * h)
}

/// Build a minimal PackContext for `nmol` single-atom molecules.
///
/// Also assigns distinct `ibmol[icart]` values so pair-penalty kernels do
/// not skip atom pairs as "same molecule". The prior version left every
/// atom at `ibmol=0`, which silently made `gradient_pair_penalty` test
/// a no-op.
fn single_atom_system(nmol: usize) -> PackContext {
    let ntotat = nmol;
    let mut sys = PackContext::new(ntotat, nmol, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![nmol];
    sys.natoms = vec![1];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    sys.coor = vec![[0.0, 0.0, 0.0]];
    sys.radius = vec![1.0; ntotat];
    sys.radius_ini = vec![1.0; ntotat];
    sys.fscale = vec![1.0; ntotat];
    for i in 0..ntotat {
        sys.ibmol[i] = i;
    }
    sys.sync_atom_props();
    sys
}

fn setup_cells(sys: &mut PackContext, cell_n: usize, cell_len: F) {
    sys.ncells = [cell_n, cell_n, cell_n];
    sys.cell_length = [cell_len; 3];
    sys.pbc_min = [0.0; 3];
    sys.pbc_length = [cell_len * cell_n as F; 3];
    sys.resize_cell_arrays();
}

/// Generic `f` ↔ `fg` finite-difference parity check for a single restraint.
///
/// Installs `restraint` on a one-atom system, places the atom at `pos`, and
/// asserts the analytic gradient matches a central finite difference on all
/// three translational DOF. `pos` MUST be chosen so every penalty branch is
/// active — the helper asserts `f(pos) > 0` so a vacuous "penalty inactive,
/// gradient trivially zero" pass cannot masquerade as a real check.
fn check_restraint_gradient(
    restraint: std::sync::Arc<dyn AtomRestraint>,
    pos: [F; 3],
    h: F,
    tol: F,
    label: &str,
) {
    let mut sys = single_atom_system(1);
    sys.restraints = vec![restraint];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    let mut x = vec![0.0; 6];
    x[0] = pos[0];
    x[1] = pos[1];
    x[2] = pos[2];

    let f0 = compute_f(&x, &mut sys);
    assert!(
        f0 > 1e-9,
        "{label}: penalty inactive at {pos:?} (f={f0}) — test would be vacuous; \
         move the atom so every penalty branch is engaged"
    );

    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < tol,
            "{label} gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── f ↔ fg parity for the remaining restraint kinds ─────────────────────────
//
// The dedicated tests above cover kinds 3/4/5/9/10/12/14. These cover the
// other seven concrete restraints so that every `AtomRestraint` impl has its
// `fg` pinned against its `f` by finite difference. Each atom is positioned
// where the penalty is genuinely active (and, for the linear-penalty box /
// cube kinds, off every median plane so the FD step does not straddle the
// gradient kink).

#[test]
fn gradient_inside_cube_constraint() {
    // Cube [0,4]³; atom past the +x/+y/+z faces → all three upper terms active.
    check_restraint_gradient(
        std::sync::Arc::new(InsideCubeRestraint::new([0.0, 0.0, 0.0], 4.0)),
        [5.0, 5.5, 6.0],
        1e-7,
        1e-5,
        "inside_cube",
    );
}

#[test]
fn gradient_outside_cube_constraint() {
    // Cube [0,4]³; atom inside, off-center toward the min corner.
    check_restraint_gradient(
        std::sync::Arc::new(OutsideCubeRestraint::new([0.0, 0.0, 0.0], 4.0)),
        [1.0, 1.3, 0.7],
        1e-7,
        1e-5,
        "outside_cube",
    );
}

#[test]
fn gradient_outside_box_constraint() {
    // Box [0,4]³; atom inside, off-center toward the min corner.
    check_restraint_gradient(
        std::sync::Arc::new(OutsideBoxRestraint::new([0.0, 0.0, 0.0], [4.0, 4.0, 4.0])),
        [1.1, 0.9, 1.4],
        1e-7,
        1e-5,
        "outside_box",
    );
}

#[test]
fn gradient_outside_sphere_constraint() {
    // Sphere r=3 at origin; atom inside → penalty active.
    check_restraint_gradient(
        std::sync::Arc::new(OutsideSphereRestraint::new([0.0, 0.0, 0.0], 3.0)),
        [1.0, 0.5, -0.4],
        1e-6,
        1e-3,
        "outside_sphere",
    );
}

#[test]
fn gradient_below_plane_constraint() {
    // Plane z = 5; atom above → penalty active.
    check_restraint_gradient(
        std::sync::Arc::new(BelowPlaneRestraint::new([0.0, 0.0, 1.0], 5.0)),
        [0.0, 0.0, 7.0],
        1e-7,
        1e-5,
        "below_plane",
    );
}

/// `OutsideCylinderRestraint` (kind 13). Its `f` is the *product* of the three
/// clamped face penalties (an AND: penalty only when the atom is inside every
/// boundary), and its `fg` is the product rule of that same product — this is
/// the deliberate asymmetry with kind-12 `InsideCylinder` (which *sums*). Place
/// the atom strictly inside all three boundaries so the product is non-zero and
/// every product-rule term contributes.
#[test]
fn gradient_outside_cylinder_constraint() {
    // Cylinder along +x, base at origin, length 4, radius 2. Atom at axial
    // w=2 (inside 0..4) and radial d=0.5 (inside r²=4): all three terms active.
    check_restraint_gradient(
        std::sync::Arc::new(OutsideCylinderRestraint::new(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            2.0,
            4.0,
        )),
        [2.0, 0.5, 0.5],
        1e-6,
        1e-3,
        "outside_cylinder",
    );
}

#[test]
fn gradient_below_gaussian_constraint() {
    // Gaussian hill σ=2, z0=0, height 3, centred at (0,0). At (1,1) the surface
    // sits at ≈2.34; atom at z=4 is above → below-gaussian penalty active.
    check_restraint_gradient(
        std::sync::Arc::new(BelowGaussianRestraint::new(0.0, 0.0, 2.0, 2.0, 0.0, 3.0)),
        [1.0, 1.0, 4.0],
        1e-6,
        1e-3,
        "below_gaussian",
    );
}

// ── pair penalty gradient ──────────────────────────────────────────────────

#[test]
fn gradient_pair_penalty() {
    let mut sys = single_atom_system(2);
    sys.restraints.clear();
    sys.iratom_offsets = vec![0, 0, 0];
    sys.iratom_data.clear();
    setup_cells(&mut sys, 1, 10.0);

    // x = [com0(3), com1(3), euler0(3), euler1(3)]
    let mut x = vec![0.0; 12];
    x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    x[3] = 2.5;
    x[4] = 1.0;
    x[5] = 1.0;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-7;
    for i in 0..6 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-3,
            "pair gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── box constraint gradient ────────────────────────────────────────────────

#[test]
fn gradient_box_constraint() {
    let mut sys = single_atom_system(1);
    sys.restraints = vec![Arc::new(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [false; 3],
    ))];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    let mut x = vec![0.0; 6];
    x[0] = 1.2;
    x[1] = -0.1;
    x[2] = 0.3;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-7;
    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-5,
            "box constraint gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── sphere constraint gradient ─────────────────────────────────────────────

#[test]
fn gradient_sphere_constraint() {
    let mut sys = single_atom_system(1);
    sys.restraints = vec![Arc::new(InsideSphereRestraint::new([0.0, 0.0, 0.0], 3.0))];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    let mut x = vec![0.0; 6];
    x[0] = 4.0;
    x[1] = 1.0;
    x[2] = 0.0;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-7;
    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-4,
            "sphere constraint gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── plane constraint gradient ──────────────────────────────────────────────

#[test]
fn gradient_above_plane_constraint() {
    let mut sys = single_atom_system(1);
    sys.restraints = vec![Arc::new(AbovePlaneRestraint::new([0.0, 0.0, 1.0], 5.0))];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    let mut x = vec![0.0; 6];
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 3.0; // below plane z=5

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-7;
    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-5,
            "above_plane gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── cylinder constraint gradient ───────────────────────────────────────────

/// Finite-difference parity for `InsideCylinderRestraint`. Cylinder is the
/// most algebraically complex single-atom restraint (axis projection +
/// radial distance + finite length), so its hand-rolled `fg` is the most
/// likely to drift relative to `f` if anyone touches `restraint.rs`. Place
/// the atom *outside* the cylinder on every axis so all three penalty
/// terms (`-w`, `w-len`, `d-r²`) are simultaneously active.
#[test]
fn gradient_inside_cylinder_constraint() {
    let mut sys = single_atom_system(1);
    // Cylinder along +x, base at origin, length 4, radius 2.
    sys.restraints = vec![Arc::new(InsideCylinderRestraint::new(
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        2.0,
        4.0,
    ))];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    // Outside on every axis: x past the end cap, off the radial axis.
    let mut x = vec![0.0; 6];
    x[0] = 6.0; // past length=4
    x[1] = 3.5; // outside radius=2
    x[2] = 0.5;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-6;
    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-3,
            "cylinder constraint gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── ellipsoid constraint gradient ──────────────────────────────────────────

/// Finite-difference parity for `InsideEllipsoidRestraint`. The penalty
/// uses anisotropic axes (a/b/c distinct), so the gradient mixes
/// per-axis division by `axis²` — most likely place for an off-by-axis
/// transcription bug.
#[test]
fn gradient_inside_ellipsoid_constraint() {
    let mut sys = single_atom_system(1);
    sys.restraints = vec![Arc::new(InsideEllipsoidRestraint::new(
        [0.0, 0.0, 0.0],
        [3.0, 2.0, 1.5],
        1.0,
    ))];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    // Outside the ellipsoid → penalty active. (3.5, 2.4, 1.7) lies just
    // beyond the surface (a1 + a2 + a3 ≈ 3.18 > 1).
    let mut x = vec![0.0; 6];
    x[0] = 3.5;
    x[1] = 2.4;
    x[2] = 1.7;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-6;
    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-3,
            "ellipsoid constraint gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── outside-ellipsoid constraint gradient (kind 9) ─────────────────────────

/// Finite-difference parity for `OutsideEllipsoidRestraint` (Packmol
/// kind 9). Until the `f` / `fg` `scale2` symmetry fix, `f` returned a
/// raw `v²` while `fg`'s gradient corresponded to `∂(scale2·v²)/∂x` —
/// at the default `scale2 = 0.01` this made the optimizer see a
/// gradient 100× flatter than `f` actually was, so the optimizer
/// declared "converged" while the function value was still high. This
/// test pins the symmetric form and would have caught the original
/// transcription bug.
#[test]
fn gradient_outside_ellipsoid_constraint() {
    let mut sys = single_atom_system(1);
    sys.restraints = vec![Arc::new(OutsideEllipsoidRestraint::new(
        [0.0, 0.0, 0.0],
        [3.0, 2.0, 1.5],
        1.0,
    ))];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    // Inside the ellipsoid → penalty active.
    let mut x = vec![0.0; 6];
    x[0] = 0.5;
    x[1] = 0.3;
    x[2] = -0.4;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-6;
    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-3,
            "outside_ellipsoid gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── gaussian constraint gradient ───────────────────────────────────────────

/// Finite-difference parity for `AboveGaussianRestraint`. The gradient
/// has a non-trivial cross term (`(z-z0)` × Gaussian profile) that is
/// easy to mis-derive; without this guard a refactor of the Fortran
/// kind-14 transcription could go unnoticed.
#[test]
fn gradient_above_gaussian_constraint() {
    let mut sys = single_atom_system(1);
    // Gaussian "hill" centred at (0,0) with σ=2, base z0=0, height 3.
    sys.restraints = vec![Arc::new(AboveGaussianRestraint::new(
        0.0, 0.0, 2.0, 2.0, 0.0, 3.0,
    ))];
    sys.iratom_offsets = vec![0, 1];
    sys.iratom_data = vec![0];
    sys.init1 = true;

    // Below the gaussian surface: at (1,1) the surface sits at
    // 3*exp(-(1+1)/8) ≈ 2.34, so z=0.5 is below → penalty active.
    let mut x = vec![0.0; 6];
    x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 0.5;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-6;
    for i in 0..3 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-3,
            "gaussian constraint gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── rotation gradient (multi-atom molecules) ───────────────────────────────

#[test]
fn gradient_with_rotations() {
    let mut sys = PackContext::new(4, 2, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![2];
    sys.natoms = vec![2];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    sys.coor = vec![[0.0, 0.0, 0.0], [1.0, 0.2, -0.1]];

    sys.radius = vec![1.0; 4];
    sys.radius_ini = vec![1.0; 4];
    sys.fscale = vec![1.0; 4];
    // 2 molecules × 2 atoms — atoms 0,1 belong to mol 0, atoms 2,3 to mol 1.
    sys.ibmol = vec![0, 0, 1, 1];
    sys.ibtype = vec![0; 4];
    sys.sync_atom_props();

    sys.restraints.clear();
    sys.iratom_offsets = vec![0, 0, 0, 0, 0];
    sys.iratom_data.clear();

    setup_cells(&mut sys, 2, 5.0);

    // x = [com0(3), com1(3), euler0(3), euler1(3)]
    let mut x = vec![0.0; 12];
    x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    x[3] = 2.1;
    x[4] = 1.4;
    x[5] = 1.3;
    x[6] = 0.3;
    x[7] = 0.5;
    x[8] = 0.7;
    x[9] = -0.4;
    x[10] = 0.2;
    x[11] = -0.6;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-7;
    for i in 0..x.len() {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 5e-3,
            "rotation gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// ── constraint + pair penalty combined ─────────────────────────────────────

#[test]
fn gradient_combined_constraint_and_pairs() {
    let mut sys = PackContext::new(3, 3, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![3];
    sys.natoms = vec![1];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    sys.coor = vec![[0.0, 0.0, 0.0]];

    sys.radius = vec![1.0; 3];
    sys.radius_ini = vec![1.0; 3];
    sys.fscale = vec![1.0; 3];
    sys.ibmol = vec![0, 1, 2];
    sys.sync_atom_props();

    // Box restraint on all atoms
    sys.restraints = vec![Arc::new(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [5.0, 5.0, 5.0],
        [false; 3],
    ))];
    sys.iratom_offsets = vec![0, 1, 1, 1]; // only first atom has constraint
    sys.iratom_data = vec![0];

    setup_cells(&mut sys, 1, 10.0);

    // x = [com0(3), com1(3), com2(3), euler0(3), euler1(3), euler2(3)]
    let mut x = vec![0.0; 18];
    x[0] = 6.0; // outside box
    x[1] = 2.0;
    x[2] = 2.0;
    x[3] = 3.0;
    x[4] = 2.0;
    x[5] = 2.0;
    x[6] = 3.5;
    x[7] = 2.5;
    x[8] = 2.0;

    let _ = compute_f(&x, &mut sys);
    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    let h = 1e-7;
    for i in 0..9 {
        let gfd = finite_diff(&x, &mut sys, i, h);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-3,
            "combined gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

#[test]
fn fused_function_and_gradient_matches_separate_evaluation() {
    let mut sys = PackContext::new(4, 2, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![2];
    sys.natoms = vec![2];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    sys.coor = vec![[0.0, 0.0, 0.0], [1.0, 0.2, -0.1]];

    sys.radius = vec![1.0; 4];
    sys.radius_ini = vec![1.0; 4];
    sys.fscale = vec![1.0; 4];
    sys.ibmol = vec![0, 0, 1, 1];
    sys.sync_atom_props();

    sys.restraints = vec![Arc::new(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [5.0, 5.0, 5.0],
        [false; 3],
    ))];
    sys.iratom_offsets = vec![0, 1, 1, 2, 2];
    sys.iratom_data = vec![0, 0];
    setup_cells(&mut sys, 2, 5.0);

    let x = vec![1.2, 1.0, 1.1, 2.4, 1.3, 1.2, 0.3, 0.5, 0.7, -0.4, 0.2, -0.6];

    let f_sep = compute_f(&x, &mut sys);
    let mut g_sep = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g_sep);

    let mut g_fused = vec![0.0; x.len()];
    let f_fused = compute_fg(&x, &mut sys, &mut g_fused);

    assert!(
        (f_sep - f_fused).abs() < 1e-10,
        "f mismatch: {f_sep} vs {f_fused}"
    );
    for (i, (&a, &b)) in g_sep.iter().zip(&g_fused).enumerate() {
        let err = (a - b).abs();
        assert!(err < 1e-10, "g mismatch at {i}: {a} vs {b} (err={err})");
    }
}
