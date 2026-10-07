//! Fixtures shared by the growth tests. Each submodule owns one part of
//! `src/grow/`; nothing here runs a growth to completion.

mod driver;
mod field;
mod internal;
mod prior;
mod refusals;

// Fixtures for the chain-growth solver tests
// (`.claude/specs/chain-growth-solver.md`), running on the `CbmcGrow`
// engine.
//
// ── Section: named rejections + result surface ────────────────────────────
//
// Covers the engine seam: the named `GrowError` rejections (no silent
// degradation, spec principle 3) and `State::degraded`. NO growth
// algorithm is exercised here. The rigid-placement layout contract lives in
// `RigidView`'s own tests (`rigid_view_layout` /
// `rigid_view_set_com_out_of_range_panics`).

use crate::grow::field::{BlockKind, OverlapField, Probe};

use crate::grow::internal::InternalTree;
use crate::test_fixtures::{chain_bonds, frame_from_parts, zigzag_coords};

use crate::grow::{GrowConfig, GrowError, TorsionPrior};

use crate::{CbmcGrow, GencanPack, IntraResidual, PackEngine, PackError, RegionRestraint, Target};
use molrs::op::F;

use molrs::core::BondDistanceWeights;

use molrs::core::Block;

use molrs::core::Frame;

use ndarray::Array1;

use rand::rngs::SmallRng;

use rand::{RngExt, SeedableRng};

use std::sync::Arc;

/// Zigzag bead chain as a `molrs::core::Frame` (see [`zigzag_coords`] /
/// [`frame_from_parts`]).
fn chain_frame(n: usize, bond_len: F, with_bonds: bool) -> Frame {
    let bonds = if with_bonds {
        chain_bonds(n)
    } else {
        Vec::new()
    };
    frame_from_parts(&zigzag_coords(n, bond_len), &bonds)
}

/// A generous periodic box that trivially fits two 5-bead chains.
const BOX_MAX: [F; 3] = [20.0, 20.0, 20.0];

// ── Section: Task 2 — InternalTree: internal-coordinate decomposition ──────
//
// Review tests for `src/grow/internal.rs` (spec Design §4a): the round-trip
// (template torsion values must rebuild the template coordinates exactly),
// the random-vars invariant (any free-variable values must preserve every
// bond length and bonded angle — the only detector for a ring-closure bond
// misclassified as a free variable), and the 1-`exclusion_depth` table.

/// Row-major identity rotation: `place_seed` under it, with `origin` at the
/// template position of `seed_atoms()[0]`, reproduces the template's own seed
/// geometry.
const IDENTITY: [[F; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn vsub(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vdot(a: [F; 3], b: [F; 3]) -> F {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vdist(a: [F; 3], b: [F; 3]) -> F {
    let d = vsub(a, b);
    vdot(d, d).sqrt()
}

/// Bonded angle `a–vertex–c` in radians.
fn vangle(a: [F; 3], vertex: [F; 3], c: [F; 3]) -> F {
    let u = vsub(a, vertex);
    let v = vsub(c, vertex);
    (vdot(u, v) / (vdot(u, u).sqrt() * vdot(v, v).sqrt()))
        .clamp(-1.0, 1.0)
        .acos()
}

fn max_abs_dev(a: &[[F; 3]], b: &[[F; 3]]) -> F {
    a.iter()
        .zip(b)
        .flat_map(|(p, q)| (0..3).map(move |k| (p[k] - q[k]).abs()))
        .fold(0.0, F::max)
}

/// Every bonded angle triple `(i, j, k)` — two bonds sharing vertex `j`.
fn bonded_angle_triples(n: usize, bonds: &[(u32, u32)]) -> Vec<(usize, usize, usize)> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(i, j) in bonds {
        adj[i as usize].push(j as usize);
        adj[j as usize].push(i as usize);
    }
    let mut out = Vec::new();
    for (j, nbrs) in adj.iter().enumerate() {
        for a in 0..nbrs.len() {
            for b in a + 1..nbrs.len() {
                out.push((nbrs[a], j, nbrs[b]));
            }
        }
    }
    out
}

/// Rebuild the full coordinate set: seed at the template's own seed position
/// under the identity rotation, then every step under `vars`.
fn rebuild_coords(tree: &InternalTree, template: &[[F; 3]], vars: &[F]) -> Vec<[F; 3]> {
    let mut coords = vec![[F::NAN; 3]; tree.n_atoms()];
    tree.place_seed(template[tree.seed_atoms()[0]], &IDENTITY, &mut coords);
    for k in 0..tree.n_steps() {
        tree.place_step(k, vars, &mut coords);
    }
    coords
}

/// The template's own value for every free variable, collected step by step
/// through the public accessors (`step_var` / `template_var`).
fn template_vars(tree: &InternalTree) -> Vec<F> {
    let mut vars = vec![0.0 as F; tree.n_vars()];
    for k in 0..tree.n_steps() {
        if let Some(v) = tree.step_var(k) {
            vars[v] = tree.template_var(k);
        }
    }
    vars
}

/// Decompose, rebuild with the template's own variable values, and assert the
/// exact round-trip plus the step-partition consistency (seed + steps place
/// every atom exactly once). Returns the tree for extra assertions.
fn assert_roundtrip(coords: &[[F; 3]], bonds: &[(u32, u32)], what: &str) -> InternalTree {
    let tree = InternalTree::from_frame(
        &frame_from_parts(coords, bonds),
        &BondDistanceWeights::from_exclusion_depth(3),
    )
    .unwrap_or_else(|e| panic!("{what}: template must decompose, got {e}"));
    assert_eq!(tree.n_atoms(), coords.len(), "{what}: n_atoms");

    let mut all: Vec<usize> = tree.seed_atoms().to_vec();
    for k in 0..tree.n_steps() {
        all.extend(tree.step_atoms(k));
    }
    all.sort_unstable();
    assert_eq!(
        all,
        (0..tree.n_atoms()).collect::<Vec<_>>(),
        "{what}: seed + step atoms must cover every atom exactly once"
    );

    let rebuilt = rebuild_coords(&tree, coords, &template_vars(&tree));
    let dev = max_abs_dev(coords, &rebuilt);
    assert!(
        dev < 1e-9,
        "{what}: template-vars round-trip ‖Δ‖∞ = {dev:e}, must be < 1e-9"
    );
    tree
}

/// 10-atom zigzag backbone with two one-atom side branches at tetrahedral-ish
/// positions off backbone atoms 3 and 6. Branch bonds are terminal, so the
/// non-terminal backbone bonds stay the only free variables.
fn branched_parts() -> (Vec<[F; 3]>, Vec<(u32, u32)>) {
    let mut coords = zigzag_coords(10, 1.53);
    let mut bonds = chain_bonds(10);
    let c3 = coords[3];
    coords.push([c3[0], c3[1] + 1.44, c3[2] + 0.51]); // atom 10
    bonds.push((3, 10));
    let c6 = coords[6];
    coords.push([c6[0], c6[1] - 1.44, c6[2] - 0.51]); // atom 11
    bonds.push((6, 11));
    (coords, bonds)
}

/// Chair-ish (non-planar) 6-ring with ~1.54 Å bonds plus a 4-atom zigzag tail
/// off ring atom 0. Ring bonds are cyclic and must never become free
/// variables; the three non-terminal tail bonds must.
fn ring_tail_parts() -> (Vec<[F; 3]>, Vec<(u32, u32)>) {
    let pi = std::f64::consts::PI as F;
    let (r, h) = (1.46 as F, 0.25 as F);
    let mut coords: Vec<[F; 3]> = (0..6)
        .map(|i| {
            let phi = i as F * pi / 3.0;
            [
                r * phi.cos(),
                r * phi.sin(),
                if i % 2 == 0 { h } else { -h },
            ]
        })
        .collect();
    let mut bonds: Vec<(u32, u32)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
    let mut tip = coords[0];
    for t in 0..4u32 {
        tip = [
            tip[0] + 1.25,
            tip[1],
            tip[2] + if t % 2 == 0 { 0.9 } else { -0.9 },
        ];
        coords.push(tip); // atoms 6, 7, 8, 9
        bonds.push(if t == 0 { (0, 6) } else { (5 + t, 6 + t) });
    }
    (coords, bonds)
}

// ── Section: TorsionPrior — sampling + C∞ calibration ─────────────────────
//
// The prior's sampling surface (spec Design §4a′):
//
//   TorsionPrior::sample(&self, template_value: F, rng: &mut impl rand::Rng) -> F
//   TorsionPrior::three_state_from_c_inf(c_inf: F, theta: F) -> TorsionPrior
//
// `sample` semantics: Uniform → uniform on (−π, π]; Template{kappa} →
// von-Mises-like around `template_value`; States → pick a state angle by
// normalized weight, exact angle, no jitter (v1).

/// Mean characteristic ratio C_n = ⟨R²⟩ / ((n−1)·b²) over `n_samples`
/// conformers of an isolated chain (end-to-end between the graph ends,
/// no PBC), every free torsion drawn independently from `prior`.
fn sampled_c_n(prior: &TorsionPrior, n_beads: usize, bond: F, n_samples: usize, seed: u64) -> F {
    let template = zigzag_coords(n_beads, bond);
    let tree = InternalTree::from_frame(
        &frame_from_parts(&template, &chain_bonds(n_beads)),
        &BondDistanceWeights::from_exclusion_depth(3),
    )
    .expect("chain template decomposes");
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut vars = vec![0.0 as F; tree.n_vars()];
    let mut sum_r2 = 0.0 as F;
    for _ in 0..n_samples {
        for k in 0..tree.n_steps() {
            if let Some(v) = tree.step_var(k) {
                vars[v] = prior.sample(tree.template_var(k), &mut rng);
            }
        }
        let rebuilt = rebuild_coords(&tree, &template, &vars);
        let r = vsub(rebuilt[n_beads - 1], rebuilt[0]);
        sum_r2 += vdot(r, r);
    }
    sum_r2 / n_samples as F / ((n_beads - 1) as F * bond * bond)
}

// ── Section: GrowStage — constructive growth end-to-end ───────────────────
//
// Exercises `src/grow/driver.rs` (spec Design §4): the `GrowConfig`
// builder surface, the `GrowError::{NoBox, TriclinicCell, FixedTarget}`
// rejections, and the constructive all-grow pack itself.
//
// Grow and GENCAN compose as an explicit chain over a fixed matrix:
// `grow_then_gencan_chaining_over_fixed_matrix` below.

// ── Section: restraint hard rejection + callback wiring ───────────────────
//
// Restraint contract (spec Design §3, ac-007): every candidate atom position is
// checked against the target's restraints via the existing
// `AtomRestraint::f`; `f > 0` is a hard rejection, same treatment as a
// hard-core violation. `frest == 0.0` thereby becomes a CONSTRUCTIVE
// guarantee, exactly like `fdist == 0.0` — strict zero, not `< precision`.
//
// Callback contract (spec Design §4g): one `StepReport` per growth round with
// `loop_idx` = round number (1-based, strictly increasing), `radscale` =
// current hard-core scale (1.0 while undegraded), and fdist/frest = 0.0
// while the hard-rejection regime holds. `Callback::should_stop() == true`
// aborts growth: `pack` still returns Ok, with `converged == false`.

// ── Section: push-off — the explicit free-target chain ────────────────────
//
// When growth ends unconverged (degraded > 0), the engine says so and stops.
// The rigid push-off is the user-explicit chain (placement-seeding spec):
// the SAME free targets go to `GencanPack::with_restart(&grown)`, whose
// phases continue on the coor/x growth wrote (Auhl slow push-off /
// Theodorou–Suter staged relaxation, spec §5.4/§5.7). The seeded run must
// (i) NOT run `initial()` — that re-randomizes every COM/Euler and
// teleports the grown chains before descent even starts — and (ii) run the
// phases with movebad disabled, so molecules move by rigid-body descent
// only.

// ── Section: density-resolved box, CG angle prior ─────────────────────────
//
// Density (spec §7, ac-005): `with_density(rho)` on the shared engine
// settings resolves in stage ① to a CUBIC periodic box `[0, L]³` (all axes
// periodic) with `L³ = total_mass / (N_A · rho)` cm³, converted to Å³ by the
// unit registry (masses in g/mol, rho in g/cm³), the total mass summing over
// ALL targets × their counts. Masses default to element lookup;
// `Target::with_mass(amu)` overrides the per-copy total (the only route for
// element-"X" targets). Named errors:
// `PackError::DensityConflictsWithBox` (density + explicit box/cell) and
// `PackError::UnknownMass { target }` (density given, a target's mass
// unresolvable, no override). Density is solver-agnostic: it belongs to the
// shared engine settings, not to `GrowConfig`.
//
// Grow and GENCAN compose as an explicit chain over a fixed matrix:
// `grow_then_gencan_chaining_over_fixed_matrix` below.
//
// CG angle prior (spec §4a′/§5.5, ac-011): `AnglePrior` makes the bond
// angle a sampling degree of freedom. `AnglePrior::Template` (the default)
// copies template angles verbatim — the AA behavior; `AnglePrior::Wlc`
// (via `wlc_from_c_inf`) is the CG path, calibrated so a discrete worm-like
// chain reproduces c∞ = (1+⟨cosθ′⟩)/(1−⟨cosθ′⟩), ⟨cosθ′⟩ = (c∞−1)/(c∞+1),
// where θ′ is the bond-deflection angle.

fn cube_tris(lo: [F; 3], hi: [F; 3]) -> Vec<[[F; 3]; 3]> {
    let p = |x, y, z| [x, y, z];
    let [x0, y0, z0] = lo;
    let [x1, y1, z1] = hi;
    vec![
        [p(x0, y0, z0), p(x0, y0, z1), p(x0, y1, z1)],
        [p(x0, y0, z0), p(x0, y1, z1), p(x0, y1, z0)],
        [p(x1, y0, z0), p(x1, y1, z0), p(x1, y1, z1)],
        [p(x1, y0, z0), p(x1, y1, z1), p(x1, y0, z1)],
        [p(x0, y0, z0), p(x1, y0, z0), p(x1, y0, z1)],
        [p(x0, y0, z0), p(x1, y0, z1), p(x0, y0, z1)],
        [p(x0, y1, z0), p(x0, y1, z1), p(x1, y1, z1)],
        [p(x0, y1, z0), p(x1, y1, z1), p(x1, y1, z0)],
        [p(x0, y0, z0), p(x0, y1, z0), p(x1, y1, z0)],
        [p(x0, y0, z0), p(x1, y1, z0), p(x1, y0, z0)],
        [p(x0, y0, z1), p(x1, y0, z1), p(x1, y1, z1)],
        [p(x0, y0, z1), p(x1, y1, z1), p(x0, y1, z1)],
    ]
}

// ── Section: the stage identity carried on StepReport ───────────────────────
