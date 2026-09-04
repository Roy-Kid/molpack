//! Contract tests for the invariant layer
//! (`.claude/specs/stage-pipeline-06-combinators.md`, `src/invariant.rs`).
//!
//! Three types are owned here (law § 11): [`Layers`](molpack::Layers) — the
//! repair-cost ladder L0–L5 as a bit set, living next to its only consumer
//! [`Invariant::layer`](molpack::Invariant::layer) — the
//! [`Invariant`](molpack::Invariant) trait itself, and the crate's first
//! concrete invariant [`RestraintsSatisfied`](molpack::RestraintsSatisfied).
//! What the *combinators* do with a violation
//! (`OnViolation::Fail` / `Rerun`) belongs to the pipeline and stays in
//! `tests/pipeline.rs`; this file never builds a `Pipeline` and never runs an
//! algorithm.
//!
//! **The fixture is a hand-built state, on purpose.** `RestraintsSatisfied`
//! reads two context fields — `frest` and `frest_atom` — and nothing else, so
//! a `PackState` wrapping `PackContext::new(ntotat, nmol, ntype)` with those
//! two fields written directly is the whole of its input. Booting a packing
//! run to produce them would make this file's answers depend on the GENCAN
//! schedule, which owns none of the behaviour under test.
//!
//! Deterministic by construction: no RNG, no wall clock, no filesystem, no
//! network, no third-party oracle.
//!
//! Single-file gate:
//!
//! ```text
//! cargo test -p molcrafts-molpack --test invariant
//! ```

use molpack::{F, Invariant, Layers, PackContext, PackState, RestraintsSatisfied};

// ── the ladder, in rung order ─────────────────────────────────────────────

/// The six rungs, cheapest-to-repair last, in the order `src/stage.rs`'s
/// module prose lists them. Indexing this array is the rung number, which is
/// what lets the loops below say "rung `i`" without a second table.
const LADDER: [Layers; 6] = [
    Layers::L0_CONNECTIVITY,
    Layers::L1_TOPOLOGICAL_STATE,
    Layers::L2_CHAIN_STATISTICS,
    Layers::L3_DENSITY,
    Layers::L4_LOCAL_OVERLAPS,
    Layers::L5_LOCAL_GEOMETRY,
];

// ── the hand-built fixture ────────────────────────────────────────────────

/// A run state whose context reports `frest` as its largest restraint
/// violation and `frest_atom[i]` for each `(i, value)` pair given.
///
/// `PackContext::new` sizes `frest_atom` to `ntotat` and zeroes it, so every
/// index not named here is exactly `0.0` — the "this atom is not part of the
/// violation" value `RestraintsSatisfied` filters on.
fn state_with_restraint_residual(
    ntotat: usize,
    nmol: usize,
    frest: F,
    per_atom: &[(usize, F)],
) -> PackState {
    let mut state = PackState::new(PackContext::new(ntotat, nmol, 1), nmol);
    let ctx = state.ctx_mut();
    ctx.frest = frest;
    for &(icart, value) in per_atom {
        ctx.frest_atom[icart] = value;
    }
    state
}

// ── 1. Layers: the bit set ────────────────────────────────────────────────

/// `contains` answers membership: a rung contains itself, contains nothing
/// else, and a union contains both of its operands. The empty set contains no
/// rung at all — which is what makes "did this invariant declare a layer?" a
/// question with an answer.
#[test]
fn layers_contains_reports_membership() {
    assert!(
        Layers::L3_DENSITY.contains(Layers::L3_DENSITY),
        "a rung must contain itself"
    );
    assert!(
        !Layers::L3_DENSITY.contains(Layers::L4_LOCAL_OVERLAPS),
        "L3 density must not claim to contain L4 local overlaps — the rungs \
         are what tells a caller how expensive the defect is to repair"
    );

    let both = Layers::L3_DENSITY.union(Layers::L4_LOCAL_OVERLAPS);
    assert!(
        both.contains(Layers::L3_DENSITY) && both.contains(Layers::L4_LOCAL_OVERLAPS),
        "a union contains both operands"
    );
    assert!(
        !both.contains(Layers::L0_CONNECTIVITY),
        "a union must not acquire a rung nobody put in it"
    );
    assert!(
        !Layers::EMPTY.contains(Layers::L3_DENSITY),
        "the empty layer set contains no rung"
    );
}

/// `union` is the set union: commutative, idempotent, and inert against
/// [`Layers::EMPTY`]. An invariant that declares two rungs must get the same
/// answer whichever order the caller wrote them in.
#[test]
fn layers_union_is_commutative_idempotent_and_has_an_identity() {
    let (a, b) = (Layers::L1_TOPOLOGICAL_STATE, Layers::L5_LOCAL_GEOMETRY);

    assert_eq!(a.union(b), b.union(a), "union must be commutative");
    assert_eq!(a.union(a), a, "union with itself changes nothing");
    assert_eq!(a.union(Layers::EMPTY), a, "EMPTY is the union identity");
    assert_eq!(
        Layers::EMPTY.union(Layers::EMPTY),
        Layers::EMPTY,
        "the empty set unions to itself"
    );
}

/// The six constants are six *distinct* bits: no rung contains another, and
/// adding a second rung always changes the set. A ladder whose rungs shared a
/// bit could not report which layer a violation sits on, which is the one
/// thing `Invariant::layer` exists for.
#[test]
fn the_six_layer_constants_are_disjoint_bits() {
    for (i, &a) in LADDER.iter().enumerate() {
        assert_ne!(a, Layers::EMPTY, "rung L{i} must be a set bit, not EMPTY");
        assert!(a.contains(a), "rung L{i} must contain itself");
        for (j, &b) in LADDER.iter().enumerate() {
            if i == j {
                continue;
            }
            assert!(
                !a.contains(b),
                "rung L{i} reports that it contains rung L{j} — the six \
                 constants must be disjoint bits"
            );
            assert_ne!(
                a.union(b),
                a,
                "adding rung L{j} to rung L{i} changed nothing — the two \
                 constants share a bit"
            );
        }
    }
}

/// `name` renders three cases: one rung gives that rung's own name (which
/// carries its rung tag, e.g. `"L3 density"`), several give `"mixed"`, none
/// gives `"none"`. The string is what the `InvariantViolated` error carries,
/// so `error.rs` never learns the `Layers` type.
#[test]
fn layers_name_renders_single_mixed_and_none() {
    assert_eq!(
        Layers::EMPTY.name(),
        "none",
        "an empty layer set renders as `none`"
    );
    assert_eq!(
        Layers::L3_DENSITY.union(Layers::L4_LOCAL_OVERLAPS).name(),
        "mixed",
        "two rungs render as `mixed` — the string names the set, it does not \
         pick a winner"
    );

    let mut seen: Vec<&'static str> = Vec::with_capacity(LADDER.len());
    for (i, &layer) in LADDER.iter().enumerate() {
        let name = layer.name();
        assert!(!name.is_empty(), "rung L{i} rendered an empty name");
        assert_ne!(name, "none", "rung L{i} must not render as the empty set");
        assert_ne!(name, "mixed", "a single rung must not render as `mixed`");
        assert!(
            name.contains(&format!("L{i}")),
            "rung L{i} rendered as `{name}`, which does not carry its rung \
             tag — the rung number is how a reader locates the defect on the \
             repair-cost ladder"
        );
        seen.push(name);
    }
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(
        seen.len(),
        LADDER.len(),
        "two rungs rendered the same name: {seen:?}"
    );
}

// ── 2. RestraintsSatisfied: what it declares ──────────────────────────────

/// The built-in invariant declares the density rung and a non-empty name.
///
/// L3 and not L4: a restraint says *where* a molecule belongs — a spatial
/// distribution defect, repairable only by moving whole molecules — while L4
/// is the local overlap a short descent removes.
#[test]
fn restraints_satisfied_declares_the_density_layer() {
    let invariant = RestraintsSatisfied::new(1e-2);

    assert_eq!(
        invariant.layer(),
        Layers::L3_DENSITY,
        "a restraint violation is a spatial-distribution defect (L3), not a \
         local overlap (L4)"
    );
    assert!(
        !invariant.name().is_empty(),
        "the invariant's name is what the InvariantViolated error points at"
    );
}

// ── 3. RestraintsSatisfied: the check ─────────────────────────────────────

/// A state whose largest restraint violation is at or below the tolerance has
/// nothing to report — including the boundary, which is *satisfied*: the
/// violation fires on `frest > tolerance`, strictly.
#[test]
fn restraints_satisfied_finds_nothing_within_tolerance() {
    let state = state_with_restraint_residual(6, 2, 0.25, &[(2, 0.25)]);

    assert!(
        RestraintsSatisfied::new(1.0).check(&state).is_empty(),
        "frest = 0.25 is well inside a tolerance of 1.0"
    );
    assert!(
        RestraintsSatisfied::new(0.25).check(&state).is_empty(),
        "frest exactly at the tolerance is satisfied — the check is `>`, not \
         `>=`, so a caller can ask for `at most this much`"
    );
}

/// Above the tolerance the invariant reports exactly one violation naming the
/// atoms that carry a restraint residual: one `Violation` for one broken
/// invariant, not one per atom.
#[test]
fn restraints_satisfied_names_the_atoms_above_tolerance() {
    let state = state_with_restraint_residual(6, 2, 0.5, &[(2, 0.5), (4, 0.125)]);

    let violations = RestraintsSatisfied::new(0.1).check(&state);

    assert_eq!(
        violations.len(),
        1,
        "one broken invariant is one violation, got {violations:?}"
    );
    assert_eq!(
        violations[0].atoms,
        vec![2, 4],
        "the violation names every atom with a non-zero restraint residual, \
         in index order"
    );
    assert!(
        !violations[0].what.is_empty(),
        "a violation must say what went wrong; an empty `what` reaches the \
         user as a blank error"
    );
}

// ── 4. regression scenario — hard-coded goldens ───────────────────────────

/// Hard-coded golden: the restraint residual a real molpack run leaves
/// behind, replayed through `RestraintsSatisfied` (ac-009).
///
/// Provenance of `GOLDEN_FREST`: molpack's own box-free GENCAN fixture — 60
/// waters held by `InsideBoxRestraint([0,0,0], [14,14,14])`, seed 42,
/// tolerance 2.0, 20 outer loops — captured 2026-09-03 from the build at
/// commit 77cba83 by printing `State::frest` in its shortest
/// round-tripping form. It is the same literal
/// `tests/pipeline.rs::pipeline_regression_single_stage_gencan_golden` pins,
/// copied by hand rather than shared, so neither file can silently move the
/// other's answer. No third-party program was involved.
///
/// **Why the state is hand-built and not the run's.** A `PackState` cannot be
/// rebuilt from a `State` — the result carries a frame and a placement
/// solution, not a context — so the golden is replayed onto a context of the
/// fixture's shape (60 molecules × 3 atoms) with the residual written into
/// `frest` and three atoms carrying it. What is pinned is the invariant's
/// answer to that number: one violation below a 1e-6 tolerance, none below
/// 1e-3, with the residual readable back to 1e-12.
#[test]
fn invariant_regression_restraints_satisfied_golden() {
    /// Largest restraint violation at termination of the box-free fixture.
    const GOLDEN_FREST: F = 0.000_301_283_868_137_773_27;
    /// The three atoms carrying the residual in this replay (chosen indices
    /// inside the fixture's 180-atom span: first water's O, eighth water's
    /// first H, last water's last H).
    const GOLDEN_ATOMS: [usize; 3] = [0, 22, 179];
    /// 60 waters × 3 atoms — the fixture's shape.
    const NTOTAT: usize = 180;
    const NMOL: usize = 60;
    /// Exact-quantity tolerance for a value that is copied, not recomputed.
    const TOL: F = 1e-12;

    let residual: Vec<(usize, F)> = GOLDEN_ATOMS
        .iter()
        .map(|&icart| (icart, GOLDEN_FREST))
        .collect();
    let state = state_with_restraint_residual(NTOTAT, NMOL, GOLDEN_FREST, &residual);

    assert!(
        (state.ctx().frest - GOLDEN_FREST).abs() < TOL,
        "the state reports frest = {} but the golden is {GOLDEN_FREST}",
        state.ctx().frest
    );

    let violations = RestraintsSatisfied::new(1e-6).check(&state);
    assert_eq!(
        violations.len(),
        1,
        "frest = {GOLDEN_FREST} is above a 1e-6 tolerance, so exactly one \
         violation is expected, got {violations:?}"
    );
    assert_eq!(
        violations[0].atoms,
        GOLDEN_ATOMS.to_vec(),
        "the violation must name the three atoms carrying the residual"
    );
    assert!(
        violations[0].what.contains(&format!("{GOLDEN_FREST}")),
        "the violation message must render the residual it was built from \
         (a rounded message hides drift), got: {}",
        violations[0].what
    );

    assert!(
        RestraintsSatisfied::new(1e-3).check(&state).is_empty(),
        "frest = {GOLDEN_FREST} is below a 1e-3 tolerance, so the same state \
         is clean at the looser ruler"
    );
}
