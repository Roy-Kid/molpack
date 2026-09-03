//! What a caller may demand of a packed state, and how expensive a defect is
//! to repair.
//!
//! An [`Invariant`] answers one question about a [`PackState`] — "is this
//! property still true?" — and, when it is not, says which atoms are involved
//! ([`Violation`]) and which rung of the repair-cost ladder the defect sits on
//! ([`Layers`]). The consumer is
//! [`Guarded`](crate::pipeline::combinators::OnViolation): it reruns the same
//! stage or fails by name, and the layer is what its error reports.
//!
//! # The repair-cost ladder
//!
//! Structural defects in a packed configuration are not equally expensive to
//! repair. The crate orders them as a six-rung ladder, cheapest to fix at the
//! bottom:
//!
//! | Rung | Defect | Why it sits there |
//! |---|---|---|
//! | **L0** | Connectivity — which atoms are bonded to which | Nothing downstream can repair a wrong bond graph; it is decided when the template is read and never again. |
//! | **L1** | Topological state — knots, entanglement, catenation | Undoing a knot needs a chain to pass through itself; no local move and no minimizer reaches it. |
//! | **L2** | Large-scale chain statistics — end-to-end distance, radius of gyration, orientation | Fixing these means re-growing a chain: reptation-scale motion, far beyond a packing run. |
//! | **L3** | Density and its homogeneity | Repairable only by moving whole molecules between regions — global, slow, but mechanical. |
//! | **L4** | Local overlaps between neighbours | The classic push-off: a short descent on the shared objective removes them. |
//! | **L5** | Bond lengths and angles | The cheapest of all — the user's force field fixes these in the first steps of minimization. |
//!
//! The ladder is prose in [`stage`](crate::stage) no longer: it is [`Layers`],
//! a bit set, and it lives here because [`Invariant::layer`] is its only
//! reader. A stage's own declarations ([`Requires`](crate::Requires) /
//! [`Guarantees`](crate::Guarantees)) are still about placement shape and
//! acquire no rung.
//!
//! # One ruler
//!
//! An invariant reads the **shared objective's** verdict off the state — the
//! numbers the context owns after a stage returns — never a second metric of
//! its own. [`RestraintsSatisfied`] therefore compares `frest`, the largest
//! restraint violation the objective computed, against its tolerance; it does
//! not re-derive restraint residuals, because two rulers for one quantity is
//! exactly what the seam exists to prevent.
//!
//! **Known limit.** `frest_atom` — the per-atom attribution
//! [`RestraintsSatisfied`] names its atoms from — is written only while
//! `movebad` is picking the worst molecules (`src/movebad.rs`); the
//! end-of-stage evaluation fills `frest` and leaves it alone. A live run can
//! therefore report a violation whose `atoms` list is empty. That is the
//! honest answer: the alternative is a second restraint metric computed here.
//!
//! # Writing your own
//!
//! Same shape as [`AtomRestraint`](crate::AtomRestraint) and
//! [`Region`](crate::Region): implement the `pub trait` on your own `pub
//! struct` and hand it to
//! [`Pipeline::with_guarded`](crate::Pipeline::with_guarded) in a
//! `Vec<Box<dyn Invariant>>`. No wrapper enum, no registry.

use molrs::types::F;

use crate::context::PackState;

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
        let (frest, tolerance) = (state.ctx().frest, self.tolerance);
        if frest <= tolerance {
            return Vec::new();
        }
        let atoms = state
            .ctx()
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
