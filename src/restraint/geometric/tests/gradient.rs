#![allow(clippy::needless_range_loop)]
//! Finite-difference gradient consistency tests.

use std::sync::Arc;

use crate::F;
use crate::PackContext;
use crate::objective::{compute_f, compute_fg, compute_g};
use crate::restraint::AtomRestraint;
use crate::restraint::geometric::{
    AbovePlaneRestraint, BelowPlaneRestraint, InsideBoxRestraint, InsideCubeRestraint,
    InsideCylinderRestraint, InsideEllipsoidRestraint, InsideSphereRestraint, OutsideBoxRestraint,
    OutsideCubeRestraint, OutsideCylinderRestraint, OutsideEllipsoidRestraint,
    OutsideSphereRestraint,
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
    // One reference conformer per copy (`coor` shares `xcart`'s index space).
    sys.coor = vec![[0.0, 0.0, 0.0]; nmol];
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
    let side = cell_len * cell_n as F;
    sys.simbox = molrs::spatial::simbox::SimBox::cube(side, molrs::types::F3::zeros(3), [false; 3])
        .expect("cell");
    sys.grid = molrs::spatial::neighbors::CellGrid::with_dims([cell_n as u32; 3], [false; 3]);
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

#[test]
fn gradient_with_rotations() {
    let mut sys = PackContext::new(4, 2, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![2];
    sys.natoms = vec![2];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    sys.coor = [[0.0, 0.0, 0.0], [1.0, 0.2, -0.1]].repeat(2);

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
    sys.coor = vec![[0.0, 0.0, 0.0]; 3];

    sys.radius = vec![1.0; 3];
    sys.radius_ini = vec![1.0; 3];
    sys.fscale = vec![1.0; 3];
    sys.ibmol = vec![0, 1, 2];
    sys.sync_atom_props();

    // Box restraint on all atoms
    sys.restraints = vec![Arc::new(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [5.0, 5.0, 5.0],
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
    sys.coor = [[0.0, 0.0, 0.0], [1.0, 0.2, -0.1]].repeat(2);

    sys.radius = vec![1.0; 4];
    sys.radius_ini = vec![1.0; 4];
    sys.fscale = vec![1.0; 4];
    sys.ibmol = vec![0, 0, 1, 1];
    sys.sync_atom_props();

    sys.restraints = vec![Arc::new(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [5.0, 5.0, 5.0],
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

// ── collective restraints through the full objective ───────────────────────
//
// The unit tests next to each collective restraint check its own `f`/`fg`
// against a finite difference on a bare coordinate array. Nothing exercised the
// step that puts it into the objective: resolving each species' slice of
// `xcart`, skipping inactive types, accumulating the coupled gradient into the
// scatter buffer, and projecting that onto the molecules' COM/Euler DOF. That
// path was silently untested — including through a rewrite of it.

#[test]
fn collective_restraint_gradient_matches_finite_difference_through_the_objective() {
    use crate::restraint::GaussianPlane;

    // Five monatomic molecules strung along z, biased toward a Gaussian
    // profile about the plane z = 0. The gradient of a distribution-matching
    // penalty is coupled across every copy, so this also checks that the
    // scatter lands on the right molecule.
    let nmol = 5;
    let mut sys = single_atom_system(nmol);
    // `init1` short-circuits both the pair terms and the collective ones, so it
    // has to be off here — with it on the check would pass while measuring
    // nothing, which is how this path stayed untested.
    sys.init1 = false;
    sys.collective = vec![(
        0usize,
        Arc::new(GaussianPlane::new([0.0, 0.0, 1.0], 0.0, 1.0, 0.0, 3.0)) as Arc<_>,
    )];

    let mut x = vec![0.0; 6 * nmol];
    // Spaced well beyond the 2 Å contact distance so the pair term stays at
    // zero and the finite difference sees the collective term alone.
    for m in 0..nmol {
        x[3 * m] = 0.5 * m as F;
        x[3 * m + 1] = -0.3 * m as F;
        x[3 * m + 2] = 6.0 * m as F - 12.0;
    }

    let f0 = compute_f(&x, &mut sys);
    assert!(
        f0 > 1e-9,
        "collective penalty inactive (f={f0}) — the check would be vacuous"
    );

    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    for i in 0..3 * nmol {
        let gfd = finite_diff(&x, &mut sys, i, 1e-6);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-5,
            "collective gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

// A centroid term reaches the molecules only through the reduction, so the
// per-copy width the objective hands it (`natoms_per_copy`) has to be right.
// A monatomic system cannot tell a correct reduction from one that treats
// every atom as its own copy — both give the same answer — so this case uses
// multi-atom copies deliberately.
#[test]
fn self_separation_gradient_matches_finite_difference_through_the_objective() {
    use crate::restraint::SelfSeparation;

    let nmol = 3;
    let natoms = 3;
    let mut sys = PackContext::new(nmol * natoms, nmol, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![nmol];
    sys.natoms = vec![natoms];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    // A small, compact 3-atom template, one reference conformer per copy.
    sys.coor = [[0.0, 0.0, 0.0], [0.6, 0.0, 0.0], [0.0, 0.6, 0.0]].repeat(nmol);
    sys.radius = vec![1.0; nmol * natoms];
    sys.radius_ini = vec![1.0; nmol * natoms];
    sys.fscale = vec![1.0; nmol * natoms];
    sys.ibmol = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
    sys.ibtype = vec![0; nmol * natoms];
    sys.sync_atom_props();
    sys.restraints.clear();
    sys.iratom_offsets = vec![0; nmol * natoms + 1];
    sys.iratom_data.clear();
    sys.init1 = false;
    setup_cells(&mut sys, 4, 6.0);

    sys.collective = vec![(0usize, Arc::new(SelfSeparation::new(10.0, 1.0)) as Arc<_>)];

    // Centroids 4 Å apart: well inside d_min = 10 so the term is active, and
    // far enough that no two atoms of different copies touch (contact is 2 Å),
    // so the finite difference sees the separation term alone.
    let mut x = vec![0.0; 6 * nmol];
    for m in 0..nmol {
        x[3 * m] = 4.0 * m as F;
        x[3 * m + 1] = 0.4 * m as F;
        x[3 * m + 2] = -0.3 * m as F;
        // Non-trivial orientations: a centroid term must be blind to them, and
        // the projection onto the Euler DOF has to agree that it is.
        x[3 * nmol + 3 * m] = 0.2 * m as F;
        x[3 * nmol + 3 * m + 1] = -0.35 * m as F;
        x[3 * nmol + 3 * m + 2] = 0.15 * m as F;
    }

    let f0 = compute_f(&x, &mut sys);
    assert!(
        f0 > 1e-9,
        "separation penalty inactive (f={f0}) — the check would be vacuous"
    );

    let mut g = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g);

    for i in 0..x.len() {
        let gfd = finite_diff(&x, &mut sys, i, 1e-6);
        let err = (g[i] - gfd).abs();
        assert!(
            err < 1e-5,
            "self-separation gradient mismatch at var {i}: analytic={} fd={gfd} err={err}",
            g[i]
        );
    }
}

#[test]
fn an_inactive_species_contributes_no_collective_gradient() {
    use crate::restraint::GaussianPlane;

    let nmol = 4;
    let mut sys = single_atom_system(nmol);
    sys.init1 = false;
    sys.collective = vec![(
        0usize,
        Arc::new(GaussianPlane::new([0.0, 0.0, 1.0], 0.0, 1.0, 0.0, 3.0)) as Arc<_>,
    )];
    let mut x = vec![0.0; 6 * nmol];
    for m in 0..nmol {
        x[3 * m + 2] = 6.0 * m as F - 9.0;
    }

    let mut g_on = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g_on);
    assert!(
        g_on.iter().any(|v| v.abs() > 1e-9),
        "expected a non-zero gradient while the species is active"
    );

    sys.comptype = vec![false];
    let mut g_off = vec![0.0; x.len()];
    compute_g(&x, &mut sys, &mut g_off);
    assert!(
        g_off.iter().all(|v| v.abs() < 1e-12),
        "an inactive species must contribute nothing"
    );
}

// ── the separation bound reaches the convergence verdict ───────────────────
//
// These are not gradient checks, but they need this file's objective-level
// harness (`single_atom_system` / `setup_cells`) and duplicating it elsewhere
// would be worse. A collective penalty lands in `f`; only a restraint that
// declares itself a *bound* also lands in `frest`, which is what the packer
// tests for convergence. Without that, an unsatisfiable request would be
// silently reported as a successful pack.

/// Three monatomic copies strung along x at `spacing`, with `SelfSeparation`.
fn separation_system(spacing: F, d_min: F) -> (PackContext, Vec<F>) {
    use crate::restraint::SelfSeparation;

    let nmol = 3;
    let mut sys = single_atom_system(nmol);
    sys.init1 = false;
    sys.restraints.clear();
    sys.iratom_offsets = vec![0; nmol + 1];
    sys.iratom_data.clear();
    setup_cells(&mut sys, 4, 8.0);
    sys.collective = vec![(0usize, Arc::new(SelfSeparation::new(d_min, 1.0)) as Arc<_>)];

    let mut x = vec![0.0; 6 * nmol];
    for m in 0..nmol {
        x[3 * m] = spacing * m as F;
    }
    (sys, x)
}

#[test]
fn a_satisfied_separation_bound_leaves_the_verdict_clean() {
    // 6 Å apart, asked for 5 — nothing owed. Atoms are also clear of the 2 Å
    // contact distance, so no other term contributes either.
    let (mut sys, x) = separation_system(6.0, 5.0);
    let f = compute_f(&x, &mut sys);
    assert_eq!(f, 0.0, "no violation, so no penalty");
    assert_eq!(sys.frest, 0.0);
}

#[test]
fn an_unsatisfiable_separation_bound_shows_up_in_frest() {
    // 3 Å apart, asked for 8 — violated, and beyond the 2 Å contact distance so
    // the pair term is silent and `frest` can only be the separation bound.
    let (mut sys, x) = separation_system(3.0, 8.0);
    let f = compute_f(&x, &mut sys);
    assert!(f > 0.0, "the bound is violated, so f must be positive");
    assert!(
        sys.frest > 0.0,
        "a violated bound must reach frest, or the run would report success"
    );
    assert_eq!(sys.fdist, 0.0, "copies are not overlapping");

    // `compute_fg` must agree with `compute_f` on both numbers.
    let mut g = vec![0.0; x.len()];
    let f_fg = compute_fg(&x, &mut sys, &mut g);
    assert!((f - f_fg).abs() < 1e-12);
    let frest_fg = sys.frest;
    assert!(
        (frest_fg - f).abs() < 1e-12,
        "one group: frest == its penalty"
    );

    // The gradient-only path shares the `fg` accumulator but must not touch the
    // verdict — the same split the per-atom value/gradient pair keeps.
    sys.frest = 0.0;
    compute_g(&x, &mut sys, &mut g);
    assert_eq!(sys.frest, 0.0, "compute_g must not write the verdict");
}

#[test]
fn a_distribution_restraint_stays_out_of_the_verdict() {
    // The asymmetry that makes the above safe: a Wasserstein penalty is always
    // positive for a finite sample, so if it counted toward `frest` no pack
    // carrying one could ever converge.
    use crate::restraint::GaussianPlane;

    let nmol = 4;
    let mut sys = single_atom_system(nmol);
    sys.init1 = false;
    sys.restraints.clear();
    sys.iratom_offsets = vec![0; nmol + 1];
    sys.iratom_data.clear();
    setup_cells(&mut sys, 4, 8.0);
    sys.collective = vec![(
        0usize,
        Arc::new(GaussianPlane::new([0.0, 0.0, 1.0], 0.0, 1.0, 0.0, 3.0)) as Arc<_>,
    )];

    let mut x = vec![0.0; 6 * nmol];
    for m in 0..nmol {
        x[3 * m + 2] = 6.0 * m as F - 9.0;
    }

    let f = compute_f(&x, &mut sys);
    assert!(f > 0.0, "the distribution penalty is active");
    assert_eq!(sys.frest, 0.0, "but it must not gate convergence");
}
