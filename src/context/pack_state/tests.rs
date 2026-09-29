//! Contract tests for [`PackState`], [`Placed`] and the crate's single
//! unscaled-verdict primitive [`evaluate_unscaled`].
//!
//! They live inside the crate because both types are `pub(crate)` for the
//! duration of the `stage-pipeline` chain and are therefore invisible to an
//! integration test in `tests/` (a separate crate). The module is mounted from
//! `src/context/pack_state.rs` with `#[cfg(test)] mod tests;` and is collected
//! by the ordinary `--lib` gate:
//!
//! ```text
//! cargo test -p molcrafts-molpack --lib --tests -- pack_state
//! ```
//!
//! Every test in this file names a fixture built by hand through
//! `PackContext::new`, in the same style as `pack_context.rs::geometry_cache_tests` and
//! `restraint::geometric::tests::gradient`, so nothing here depends on the `.inp` front end, on
//! `initial`, or on a solver driver.

use std::sync::Arc;

use molrs::types::F;

use super::{PackState, Placed, evaluate_unscaled};
use crate::constraints::EvalMode;
use crate::context::{PackContext, RigidView};
use crate::numerics::DEFAULT_SCALE2;
use crate::restraint::AtomRestraint;
use crate::restraint::geometric::InsideBoxRestraint;

// ── the shared fixture ─────────────────────────────────────────────────────

/// Free molecules in the fixture.
const NMOL: usize = 6;
/// Atoms per molecule — a dimer.
const NAT: usize = 2;
/// Total atoms (all free; the fixture carries no fixed structure).
const NTOTAT: usize = NMOL * NAT;
/// Edge of the cubic packing cell, in Angstrom.
const BOX: F = 20.0;
/// Placement-stream seed (see [`Lcg`]).
const SEED: u64 = 7;
/// Dimer bond length, in Angstrom.
const BOND: F = 1.5;
/// Unscaled atomic radius (`radius_ini`), in Angstrom.
const RADIUS: F = 3.0;
/// Radius factor standing in for the GENCAN `discale` schedule, i.e. the
/// state in which `radius != radius_ini` and the swap inside the unscaled
/// verdict actually moves data.
const DISCALE: F = 0.5;
/// Low corner of the centre-of-mass draw. Deliberately outside the cell so
/// some atoms violate the box restraint and `frest` is non-zero: with
/// `frest == 0` the `scale` / `scale2` fields would not reach the verdict at
/// all (they enter only through the restraint terms, `src/objective.rs:488`,
/// `:606`, `:733`) and the symmetry test below would be vacuous.
const COM_LO: F = -1.0;
/// Span of the centre-of-mass draw. Tighter than the cell so the six dimers
/// crowd each other and `fdist` is non-zero as well.
const COM_SPAN: F = 16.0;

/// A 64-bit LCG (Knuth's MMIX multiplier / increment) yielding `[0, 1)`.
///
/// Hand-rolled on purpose: the fixture's numbers must not move when the
/// `rand` dependency bumps its generator's stream, because the hard-coded
/// goldens in [`pack_state_regression_unscaled_verdict_golden`] are captured
/// from exactly this sequence.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Next uniform draw in `[0, 1)`, from the top 53 bits of the state.
    fn unit(&mut self) -> F {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// Six rigid dimers in a 20 Angstrom cubic cell, placed from [`SEED`].
///
/// Returns the context and the flat placement vector `x` (`6 * NMOL` values,
/// COM block then Euler block) that produced its `xcart`: the two are written
/// through [`RigidView::write_xcart`], so the lab-frame coordinates and the
/// placement vector describe the same geometry rather than merely coexisting.
///
/// On return `radius == radius_ini` (the growth-path state) and
/// `scale` / `scale2` hold their constructor defaults; the GENCAN-path state
/// is this fixture plus [`scale_radius`].
fn six_dimers() -> (PackContext, Vec<F>) {
    let mut ctx = PackContext::new(NTOTAT, NMOL, 1);
    ctx.ntype_with_fixed = 1;
    ctx.nmols = vec![NMOL];
    ctx.natoms = vec![NAT];
    ctx.idfirst = vec![0];
    ctx.comptype = vec![true];

    // One reference conformer per copy (`coor` shares `xcart`'s index space):
    // each dimer's two atoms sit on the x axis about their own centre.
    ctx.coor = (0..NTOTAT)
        .map(|icart| {
            let sign = if icart % NAT == 0 { -0.5 } else { 0.5 };
            [sign * BOND, 0.0, 0.0]
        })
        .collect();

    ctx.radius = vec![RADIUS; NTOTAT];
    ctx.radius_ini = vec![RADIUS; NTOTAT];
    ctx.fscale = vec![1.0; NTOTAT];
    for (icart, slot) in ctx.ibmol.iter_mut().enumerate() {
        *slot = icart / NAT;
    }
    ctx.ibtype.fill(0);

    let box_restraint: Arc<dyn AtomRestraint> =
        Arc::new(InsideBoxRestraint::new([0.0; 3], [BOX; 3]));
    ctx.restraints = vec![box_restraint];
    ctx.iratom_offsets = (0..=NTOTAT).collect();
    ctx.iratom_data = vec![0; NTOTAT];

    ctx.sizemin = [0.0; 3];
    ctx.sizemax = [BOX; 3];

    ctx.simbox = molrs::spatial::simbox::SimBox::cube(BOX, molrs::types::F3::zeros(3), [false; 3])
        .expect("cubic packing cell");
    ctx.grid = molrs::spatial::neighbors::CellGrid::with_dims([2; 3], [false; 3]);
    ctx.resize_cell_arrays();
    ctx.sync_atom_props();

    let mut rng = Lcg::new(SEED);
    let mut view = RigidView::fresh(NMOL);
    for imol in 0..NMOL {
        let com = [
            COM_LO + rng.unit() * COM_SPAN,
            COM_LO + rng.unit() * COM_SPAN,
            COM_LO + rng.unit() * COM_SPAN,
        ];
        view.set_com(imol, com);
    }
    for imol in 0..NMOL {
        let euler = [
            rng.unit() * std::f64::consts::TAU as F,
            rng.unit() * std::f64::consts::TAU as F,
            rng.unit() * std::f64::consts::TAU as F,
        ];
        view.set_euler(imol, euler);
    }
    view.write_xcart(&mut ctx);

    let x = view.as_slice().to_vec();
    (ctx, x)
}

/// Put the fixture into the GENCAN-path state: `radius = factor * radius_ini`,
/// mirrored into `atom_props` through the setter.
fn scale_radius(ctx: &mut PackContext, factor: F) {
    let scaled: Vec<F> = ctx.radius_ini.iter().map(|r| r * factor).collect();
    for (icart, r) in scaled.into_iter().enumerate() {
        ctx.set_radius(icart, r);
    }
}

fn radius_bits(ctx: &PackContext) -> Vec<u64> {
    ctx.radius.iter().map(|r| r.to_bits()).collect()
}

fn bits(triple: (F, F, F)) -> (u64, u64, u64) {
    (triple.0.to_bits(), triple.1.to_bits(), triple.2.to_bits())
}

/// **The parity oracle: the pre-03 unscaled verdict, owned by this test file.**
///
/// This is the body of `gencan::phases::evaluate_unscaled` as it stood before
/// stage-pipeline-03, transcribed verbatim from `src/gencan/phases.rs:33-45`
/// (captured 2026-09-03 from commit ef87105): swap `radius_ini` into
/// `radius` through `work.radiuswork`, evaluate `FOnly`, swap back. It reads
/// neither `scale` nor `scale2`, which is precisely the asymmetry ac-003
/// exists to remove.
///
/// **Do not "simplify" this into a call to the production function.** The
/// implementer deletes the original definition in this same spec, so a helper
/// that forwarded to the merged `evaluate_unscaled` would make the two ac-004
/// parity tests compare the implementation with itself and certify nothing.
/// An oracle has to be a copy that stops moving; the copy is the point. If
/// this body ever needs to change, the change is a behaviour change and must
/// be argued for, not applied.
fn legacy_unscaled(ctx: &mut PackContext, x: &[F]) -> (F, F, F) {
    ctx.work.radiuswork.copy_from_slice(&ctx.radius);
    // `i` is both the argument to the setter and the index into the source
    // array, so this is not a `needless_range_loop`.
    for i in 0..ctx.ntotat {
        ctx.set_radius(i, ctx.radius_ini[i]);
    }
    let f_total = ctx.evaluate(x, EvalMode::FOnly, None).f_total;
    let fdist = ctx.fdist;
    let frest = ctx.frest;
    for i in 0..ctx.ntotat {
        ctx.set_radius(i, ctx.work.radiuswork[i]);
    }
    (f_total, fdist, frest)
}

// ── 1. wrapping semantics ──────────────────────────────────────────────────

#[test]
fn new_wraps_the_same_fixedatom_and_comptype_storage() {
    let (ctx, _) = six_dimers();
    assert!(!ctx.fixedatom[3], "fixture starts with every atom free");
    assert!(ctx.comptype[0], "fixture starts with its one type active");

    let mut state = PackState::new(ctx, NMOL);

    // `PackState` wraps rather than extracts: a write through `ctx_mut` is
    // visible through `ctx`, because there is one storage, not a mirror.
    state.ctx_mut().set_fixed_atom(3, true);
    state.ctx_mut().comptype[0] = false;

    assert!(
        state.ctx().fixedatom[3],
        "fixedatom must be the wrapped context's own vector, not a copy"
    );
    assert!(
        !state.ctx().comptype[0],
        "comptype must be the wrapped context's own vector, not a copy"
    );
    assert_eq!(
        state.ctx().fixedatom.len(),
        NTOTAT,
        "wrapping must not resize the anchored-atom set"
    );
}

#[test]
fn new_starts_with_placed_none_and_a_fresh_rigid_view() {
    let (ctx, _) = six_dimers();
    let state = PackState::new(ctx, NMOL);

    assert_eq!(
        state.placed(),
        Placed::None,
        "a freshly wrapped state has placed nothing"
    );
    assert_eq!(
        state.rigid().nmol(),
        NMOL,
        "the rigid slot is sized from the constructor argument"
    );
    assert_eq!(
        state.rigid().as_slice().len(),
        6 * NMOL,
        "a fresh view holds 3 COM + 3 Euler values per molecule"
    );
    assert!(
        state.rigid().as_slice().iter().all(|v| *v == 0.0),
        "RigidView::fresh is zeroed"
    );
}

// ── 2. the shape marker ────────────────────────────────────────────────────

#[test]
fn set_placed_moves_the_shape_marker_to_all() {
    let (ctx, _) = six_dimers();
    let mut state = PackState::new(ctx, NMOL);

    state.set_placed(Placed::All);

    assert_eq!(state.placed(), Placed::All);
}

// ── 3. the split borrow ────────────────────────────────────────────────────

#[test]
fn rigid_split_mut_hands_out_two_disjoint_mutable_borrows() {
    let (ctx, _) = six_dimers();
    let mut state = PackState::new(ctx, NMOL);

    // Both `&mut` must be live at once — this is the whole point of the
    // accessor, and a signature that returns them sequentially would not
    // compile here.
    {
        let (ctx, rigid) = state.rigid_split_mut();
        ctx.scale = 0.25;
        rigid.set_com(1, [1.5, 2.5, 3.5]);
        ctx.scale2 = 0.125;
        rigid.set_euler(1, [0.25, 0.5, 0.75]);
    }

    let (ctx, rigid) = state.into_parts();
    assert_eq!(ctx.scale.to_bits(), (0.25 as F).to_bits());
    assert_eq!(ctx.scale2.to_bits(), (0.125 as F).to_bits());
    assert_eq!(rigid.com(1), [1.5, 2.5, 3.5]);
    assert_eq!(rigid.euler(1), [0.25, 0.5, 0.75]);
}

// ── 4. giving the parts back ───────────────────────────────────────────────

#[test]
fn into_parts_returns_the_wrapped_context_and_view() {
    let (ctx, _) = six_dimers();
    let mut state = PackState::new(ctx, NMOL);
    state.ctx_mut().frest = 4.25;
    state.rigid_split_mut().1.set_com(0, [7.0, 8.0, 9.0]);

    let (ctx, rigid) = state.into_parts();

    assert_eq!(
        ctx.ntotat, NTOTAT,
        "into_parts hands back the context that was wrapped"
    );
    assert_eq!(ctx.frest.to_bits(), (4.25 as F).to_bits());
    assert_eq!(rigid.nmol(), NMOL);
    assert_eq!(rigid.com(0), [7.0, 8.0, 9.0]);
}

// ── 5. the geometry-cache forwarder ────────────────────────────────────────

#[test]
fn invalidate_geometry_cache_forwards_to_the_context() {
    let (mut ctx, x) = six_dimers();
    // One evaluation populates the cached Cartesian expansion; the same probe
    // `pack_context.rs::geometry_cache_tests` uses.
    let _ = ctx.evaluate(&x, EvalMode::FOnly, None);
    assert!(
        ctx.work.cached_geometry.is_some(),
        "fixture must warm the cache, or the test below is vacuous"
    );

    let mut state = PackState::new(ctx, NMOL);
    state.invalidate_geometry_cache();

    assert!(
        state.ctx().work.cached_geometry.is_none(),
        "PackState::invalidate_geometry_cache must reach PackContext's cache"
    );
}

// ── 6. ac-003: save/restore symmetry ───────────────────────────────────────

#[test]
fn evaluate_unscaled_restores_scale_and_radius() {
    let (mut ctx, x) = six_dimers();
    scale_radius(&mut ctx, DISCALE);
    // An artificial state: neither field holds its default, so a helper that
    // only *sets* them (and never restores) is caught, and so is one that
    // never sets them at all.
    ctx.scale = 0.7;
    ctx.scale2 = 0.05;

    let scale_before = ctx.scale.to_bits();
    let scale2_before = ctx.scale2.to_bits();
    let radius_before = radius_bits(&ctx);

    let got = evaluate_unscaled(&mut ctx, &x);

    // (a) `scale` / `scale2` are restored as symmetrically as `radius`.
    assert_eq!(
        ctx.scale.to_bits(),
        scale_before,
        "evaluate_unscaled must give the caller's scale back"
    );
    assert_eq!(
        ctx.scale2.to_bits(),
        scale2_before,
        "evaluate_unscaled must give the caller's scale2 back"
    );
    // (b) the radius swap does not leak.
    assert_eq!(
        radius_bits(&ctx),
        radius_before,
        "evaluate_unscaled must give the caller's radius back"
    );

    // (c) the returned triple is the *unscaled* verdict: the one the same
    // fixture produces at scale = 1.0, scale2 = DEFAULT_SCALE2 and
    // radius = radius_ini, evaluated straight through the objective. Without
    // this the assertions above would pass on a helper that simply never
    // touches `scale` / `scale2` and hands back the artificial-scale verdict.
    let (mut reference, _) = six_dimers();
    reference.scale = 1.0;
    reference.scale2 = DEFAULT_SCALE2;
    assert_eq!(
        radius_bits(&reference),
        reference
            .radius_ini
            .iter()
            .map(|r| r.to_bits())
            .collect::<Vec<u64>>(),
        "the reference must evaluate at the unscaled radii"
    );
    let f_total = reference.evaluate(&x, EvalMode::FOnly, None).f_total;
    let want = (f_total, reference.fdist, reference.frest);

    assert_eq!(
        bits(got),
        bits(want),
        "unscaled verdict mismatch: got {got:?}, want {want:?}"
    );
}

// ── 7. ac-004: value parity with the pre-move implementation ───────────────

#[test]
fn evaluate_unscaled_matches_legacy_on_gencan_fixture() {
    let (mut merged, x) = six_dimers();
    scale_radius(&mut merged, DISCALE);
    let (mut legacy, _) = six_dimers();
    scale_radius(&mut legacy, DISCALE);

    let got = evaluate_unscaled(&mut merged, &x);
    let want = legacy_unscaled(&mut legacy, &x);

    assert_eq!(
        bits(got),
        bits(want),
        "GENCAN fixture (radius scaled by discale): got {got:?}, want {want:?}"
    );
    assert_eq!(
        radius_bits(&merged),
        radius_bits(&legacy),
        "post-call radius must match the pre-move implementation bit for bit"
    );
}

#[test]
fn evaluate_unscaled_matches_legacy_on_growth_fixture() {
    let (mut merged, x) = six_dimers();
    let (mut legacy, _) = six_dimers();
    assert_eq!(
        radius_bits(&merged),
        merged
            .radius_ini
            .iter()
            .map(|r| r.to_bits())
            .collect::<Vec<u64>>(),
        "the growth fixture runs at radius == radius_ini"
    );

    let got = evaluate_unscaled(&mut merged, &x);
    let want = legacy_unscaled(&mut legacy, &x);

    assert_eq!(
        bits(got),
        bits(want),
        "growth fixture (radius == radius_ini): got {got:?}, want {want:?}"
    );
    assert_eq!(
        radius_bits(&merged),
        radius_bits(&legacy),
        "post-call radius must match the pre-move implementation bit for bit"
    );
}

#[test]
fn pack_state_evaluate_unscaled_forwards_to_the_free_function() {
    let (ctx, x) = six_dimers();
    let mut state = PackState::new(ctx, NMOL);
    let via_state = state.evaluate_unscaled(&x);

    let (mut bare, _) = six_dimers();
    let via_free = evaluate_unscaled(&mut bare, &x);

    assert_eq!(
        bits(via_state),
        bits(via_free),
        "the method must be a thin forward, not a second implementation"
    );
}

// ── edge: the degenerate context ───────────────────────────────────────────

#[test]
fn evaluate_unscaled_on_empty_context_returns_zeros() {
    let mut state = PackState::new(PackContext::new(0, 0, 0), 0);
    assert_eq!(state.rigid().nmol(), 0);

    let (f_total, fdist, frest) = state.evaluate_unscaled(&[]);

    assert!(f_total.is_finite() && fdist.is_finite() && frest.is_finite());
    assert_eq!(f_total.to_bits(), (0.0 as F).to_bits());
    assert_eq!(fdist.to_bits(), (0.0 as F).to_bits());
    assert_eq!(frest.to_bits(), (0.0 as F).to_bits());
}
