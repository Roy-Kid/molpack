//! Run invariants: [`Invariant`](crate::Invariant) and its documentation.

use molrs::op::F;

use crate::system::PackState;

/// A set of rungs on the repair-cost ladder (see the module docs).
///
/// A set and not an enum: an invariant may guard a defect that spans rungs.
/// The six constants are single, disjoint bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Layers(u8);

impl Layers {
    /// No rung at all — renders as `"none"`.
    pub const EMPTY: Layers = Layers(0);
    /// L0: which atoms are bonded to which.
    pub const L0_CONNECTIVITY: Layers = Layers(1 << 0);
    /// L1: knots, entanglement, catenation.
    pub const L1_TOPOLOGICAL_STATE: Layers = Layers(1 << 1);
    /// L2: end-to-end distance, radius of gyration, orientation.
    pub const L2_CHAIN_STATISTICS: Layers = Layers(1 << 2);
    /// L3: density and its homogeneity — where a molecule belongs.
    pub const L3_DENSITY: Layers = Layers(1 << 3);
    /// L4: local overlaps between neighbours.
    pub const L4_LOCAL_OVERLAPS: Layers = Layers(1 << 4);
    /// L5: bond lengths and angles.
    pub const L5_LOCAL_GEOMETRY: Layers = Layers(1 << 5);

    /// The one table, read by [`name`](Self::name) and nothing else.
    const LADDER: [(Layers, &'static str); 6] = [
        (Layers::L0_CONNECTIVITY, "L0 connectivity"),
        (Layers::L1_TOPOLOGICAL_STATE, "L1 topological state"),
        (Layers::L2_CHAIN_STATISTICS, "L2 chain statistics"),
        (Layers::L3_DENSITY, "L3 density"),
        (Layers::L4_LOCAL_OVERLAPS, "L4 local overlaps"),
        (Layers::L5_LOCAL_GEOMETRY, "L5 local geometry"),
    ];

    /// Subset test: every rung of `other` is also in `self`.
    pub fn contains(self, other: Layers) -> bool {
        self.0 & other.0 == other.0
    }

    /// The union of two rung sets.
    pub fn union(self, other: Layers) -> Layers {
        Layers(self.0 | other.0)
    }

    /// The set rendered for a human: one rung gives that rung's name
    /// (`"L3 density"`), several give `"mixed"`, none gives `"none"`.
    ///
    /// The rendering is what
    /// [`PackError::InvariantViolated`](crate::PackError::InvariantViolated)
    /// carries, which is how the error layer stays free of a dependency on
    /// this module.
    pub fn name(self) -> &'static str {
        if self == Layers::EMPTY {
            return "none";
        }
        Layers::LADDER
            .iter()
            .find(|(rung, _)| *rung == self)
            .map_or("mixed", |(_, name)| *name)
    }
}

/// A property of a packed state a caller can demand, and a guard enforce.
///
/// `Send`, like [`Stage`](crate::Stage): a run may move across threads.
///
/// What a caller may demand of a packed state, and how expensive a defect is
/// to repair.
///
/// An [`Invariant`] answers one question about a [`PackState`] — "is this
/// property still true?" — and, when it is not, says which atoms are involved
/// ([`Violation`]) and which rung of the repair-cost ladder the defect sits on
/// ([`Layers`]). The consumer is
/// [`Guarded`](crate::OnViolation): it reruns the same
/// stage or fails by name, and the layer is what its error reports.
///
/// # The repair-cost ladder
///
/// Structural defects in a packed configuration are not equally expensive to
/// repair. The crate orders them as a six-rung ladder, cheapest to fix at the
/// bottom:
///
/// | Rung | Defect | Why it sits there |
/// |---|---|---|
/// | **L0** | Connectivity — which atoms are bonded to which | Nothing downstream can repair a wrong bond graph; it is decided when the template is read and never again. |
/// | **L1** | Topological state — knots, entanglement, catenation | Undoing a knot needs a chain to pass through itself; no local move and no minimizer reaches it. |
/// | **L2** | Large-scale chain statistics — end-to-end distance, radius of gyration, orientation | Fixing these means re-growing a chain: reptation-scale motion, far beyond a packing run. |
/// | **L3** | Density and its homogeneity | Repairable only by moving whole molecules between regions — global, slow, but mechanical. |
/// | **L4** | Local overlaps between neighbours | The classic push-off: a short descent on the shared objective removes them. |
/// | **L5** | Bond lengths and angles | The cheapest of all — the user's force field fixes these in the first steps of minimization. |
///
/// The ladder is [`Layers`], a bit set, and it lives here because [`Invariant::layer`] is its only
/// reader. A stage's own declarations ([`Requires`](crate::Requires) /
/// [`Guarantees`](crate::Guarantees)) are still about placement shape and
/// acquire no rung.
///
/// # One ruler
///
/// An invariant reads the **shared objective's** verdict off the state — the
/// numbers the system owns after a stage returns — never a second metric of
/// its own. [`RestraintsSatisfied`] therefore compares `frest`, the largest
/// restraint violation the objective computed, against its tolerance; it does
/// not re-derive restraint residuals, because two rulers for one quantity is
/// exactly what the seam exists to prevent.
///
/// **Known limit.** `frest_atom` — the per-atom attribution
/// [`RestraintsSatisfied`] names its atoms from — is written only while
/// `movebad` is picking the worst molecules (`src/movebad.rs`); the
/// end-of-stage evaluation fills `frest` and leaves it alone. A live run can
/// therefore report a violation whose `atoms` list is empty. That is the
/// honest answer: the alternative is a second restraint metric computed here.
///
/// # Writing your own
///
/// Same shape as [`AtomRestraint`](crate::AtomRestraint) and molrs's
/// [`Region`](molrs::core::Region): implement the `pub trait` on
/// your own `pub struct` and hand it to
/// [`Pipeline::with_guarded`](crate::Pipeline::with_guarded) in a
/// `Vec<Box<dyn Invariant>>`. No wrapper enum, no registry.
pub trait Invariant: Send {
    /// Short identifier for reports and for
    /// [`PackError::InvariantViolated`](crate::PackError::InvariantViolated).
    fn name(&self) -> &'static str;

    /// The rung(s) a violation of this invariant sits on — how expensive the
    /// defect is to repair, which is what tells a caller whether rerunning
    /// the same stage can even help.
    fn layer(&self) -> Layers;

    /// Every way `state` breaks this invariant, or an empty vector when it
    /// holds. One broken invariant is normally one [`Violation`]; the vector
    /// is for an invariant that can name independent defects.
    fn check(&self, state: &PackState) -> Vec<Violation>;
}

/// One way a state broke an invariant: the atoms involved and what happened.
#[derive(Debug, Clone)]
pub struct Violation {
    /// The atoms the defect involves, as indices into the run's own atoms.
    /// May be empty when the invariant knows *that* it broke but not *where*.
    pub atoms: Vec<usize>,
    /// What went wrong, rendered for a human — the text that reaches the user
    /// through a log line or an error.
    pub what: String,
}

/// The largest restraint violation the shared objective left on the state is
/// within `tolerance`.
///
/// L3, not L4: a restraint says *where* a molecule belongs — a
/// spatial-distribution defect repairable only by moving whole molecules —
/// while a local overlap is what a short descent removes.
#[derive(Debug, Clone)]
pub struct RestraintsSatisfied {
    tolerance: F,
}

impl RestraintsSatisfied {
    /// Satisfied while `frest <= tolerance`; strictly above it is a
    /// violation, so a caller can ask for "at most this much".
    pub fn new(tolerance: F) -> Self {
        Self { tolerance }
    }
}

impl Invariant for RestraintsSatisfied {
    fn name(&self) -> &'static str {
        "restraints-satisfied"
    }

    fn layer(&self) -> Layers {
        Layers::L3_DENSITY
    }

    /// Reads `frest` — the shared objective's own number — and attributes it
    /// to the atoms carrying a residual. See the module docs for why that
    /// list can be empty on a live run.
    fn check(&self, state: &PackState) -> Vec<Violation> {
        let (frest, tolerance) = (state.sys().frest, self.tolerance);
        if frest <= tolerance {
            return Vec::new();
        }
        let atoms = state
            .sys()
            .frest_atom
            .iter()
            .enumerate()
            .filter(|&(_, &residual)| residual > 0.0)
            .map(|(icart, _)| icart)
            .collect();
        vec![Violation {
            atoms,
            what: format!("frest {frest} exceeds tolerance {tolerance}"),
        }]
    }
}

#[cfg(test)]
mod tests {
    //! Contract tests for the invariant layer.
    //!
    //! Three types are owned here (law § 11): [`Layers`](crate::Layers) — the
    //! repair-cost ladder L0–L5 as a bit set, living next to its only consumer
    //! [`Invariant::layer`](crate::Invariant::layer) — the
    //! [`Invariant`](crate::Invariant) trait itself, and the crate's first
    //! concrete invariant [`RestraintsSatisfied`](crate::RestraintsSatisfied).
    //! What the *combinators* do with a violation
    //! (`OnViolation::Fail` / `Rerun`) belongs to the pipeline and stays in
    //! `pipeline::tests`; this file never builds a `Pipeline` and never runs an
    //! algorithm.
    //!
    //! **The fixture is a hand-built state, on purpose.** `RestraintsSatisfied`
    //! reads two system fields — `frest` and `frest_atom` — and nothing else, so
    //! a `PackState` wrapping `PackSystem::new(ntotat, nmol, ntype)` with those
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
    //! cargo test -p molcrafts-molpack --lib
    //! ```

    use crate::{Invariant, Layers, PackState, PackSystem, RestraintsSatisfied};
    use molrs::op::F;

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

    /// A run state whose system reports `frest` as its largest restraint
    /// violation and `frest_atom[i]` for each `(i, value)` pair given.
    ///
    /// `PackSystem::new` sizes `frest_atom` to `ntotat` and zeroes it, so every
    /// index not named here is exactly `0.0` — the "this atom is not part of the
    /// violation" value `RestraintsSatisfied` filters on.
    fn state_with_restraint_residual(
        ntotat: usize,
        nmol: usize,
        frest: F,
        per_atom: &[(usize, F)],
    ) -> PackState {
        let mut state = PackState::new(PackSystem::new(ntotat, nmol, 1), nmol);
        let sys = state.sys_mut();
        sys.frest = frest;
        for &(icart, value) in per_atom {
            sys.frest_atom[icart] = value;
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
}
