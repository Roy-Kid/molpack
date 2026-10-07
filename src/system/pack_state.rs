//! [`PackState`]: the state one packing run carries between stages.
//!
//! The type **wraps** a [`PackSystem`] instead of extracting fields out of
//! it. Not one field moves, so the shared objective's read set is unchanged,
//! `sys.xcart` stays the single lab-frame home of the coordinates, and
//! `sys.fixedatom` / `sys.comptype` are read through the wrapper rather than
//! mirrored — a mirror can drift, which is why
//! `PackSystem::debug_assert_atom_props_sync` exists at all. On top of the
//! system the state adds the two things a chain of stages needs and the
//! system has no home for: the placement shape marker [`Placed`] and the
//! rigid placement slot [`RigidView`].
//!
//! # What the state deliberately does not carry
//!
//! * **No `topology` field.** Nothing here produces or consumes one; a
//!   system-level `Topology::tile` arrives with the deferred refinement work,
//!   which brings its producer and its consumer at the same time. An empty
//!   field now would be complexity nothing pays for.
//! * **No per-atom placed bitset.** `OverlapField::is_placed`
//!   (`src/grow/field.rs`) is this crate's live placement bitset, maintained
//!   incrementally by insert / retract. A second one would be a second truth
//!   for one predicate inside a single run. Chained stages need only the
//!   state-level shape marker [`Placed`], advanced by a stage's guarantees.
//! * **No `Option` around `rigid`.** [`PackState::new`] installs
//!   `RigidView::fresh(nmol)` immediately, so "is there a view in the slot?"
//!   is not a representable question.
//!
//! # Why `Debug` is hand-written
//!
//! [`PackSystem`] holds `Arc<dyn Restraint>` values, a cell grid and a
//! frame, and does not implement `Debug`, so `#[derive(Debug)]` on a struct
//! that owns one cannot compile. Rather than grow the system with a derive
//! it does not otherwise need, this module writes the impl by hand: it prints
//! `placed` and the view's molecule count and elides the system.

use std::fmt;

use molrs::op::F;

use crate::Objective;
use crate::eval::EvalMode;
use crate::system::DEFAULT_SCALE2;
use crate::system::{PackSystem, RigidView};

/// The shape of the placements a [`PackState`] currently holds.
///
/// A state-level marker, not a per-atom set: it answers "has every free
/// molecule been placed yet?", which is what a chain of stages checks before
/// running the next one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placed {
    /// Nothing has been placed yet — a freshly wrapped system.
    None,
    /// Every free molecule the system lays out has a placement.
    All,
}

/// A [`PackSystem`] plus the two pieces of run state a stage chain needs.
pub struct PackState {
    /// The wrapped system — the authority for geometry, radii and restraints.
    sys: PackSystem,
    /// The shape marker; see [`Placed`].
    placed: Placed,
    /// The rigid degrees of freedom, always present (never an `Option`).
    rigid: RigidView,
}

impl fmt::Debug for PackState {
    /// Prints the marker and the view's molecule count; the system is
    /// elided because it does not implement `Debug` (see the module docs).
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PackState")
            .field("placed", &self.placed)
            .field("nmol", &self.rigid.nmol())
            .finish_non_exhaustive()
    }
}

impl PackState {
    /// Wrap `sys` into a fresh run state for `nmol` free molecules: a zeroed
    /// [`RigidView`] in the slot and [`Placed::None`] as the marker.
    pub fn new(sys: PackSystem, nmol: usize) -> Self {
        Self {
            sys,
            placed: Placed::None,
            rigid: RigidView::fresh(nmol),
        }
    }

    /// The wrapped packing system.
    pub fn sys(&self) -> &PackSystem {
        &self.sys
    }

    /// The wrapped packing system, mutably.
    pub fn sys_mut(&mut self) -> &mut PackSystem {
        &mut self.sys
    }

    /// The current placement shape marker.
    pub fn placed(&self) -> Placed {
        self.placed
    }

    /// Advance the placement shape marker.
    pub fn set_placed(&mut self, placed: Placed) {
        self.placed = placed;
    }

    /// The rigid placement vector.
    pub fn rigid(&self) -> &RigidView {
        &self.rigid
    }

    /// System and view as two disjoint mutable borrows, for the many
    /// operations that write placements while reading or updating the
    /// system (`write_xcart`, `capture_from_xcart`, a stage's inner loop).
    pub fn rigid_split_mut(&mut self) -> (&mut PackSystem, &mut RigidView) {
        (&mut self.sys, &mut self.rigid)
    }

    /// Give the two owned parts back, consuming the state.
    pub fn into_parts(self) -> (PackSystem, RigidView) {
        (self.sys, self.rigid)
    }

    /// Drop the system's cached Cartesian expansion.
    ///
    /// A forward to [`PackSystem::invalidate_geometry_cache`], not a second
    /// implementation: the cache and its invalidation semantics belong to the
    /// system. It sits here because the caller is a stage boundary, which
    /// holds the state and not the bare system.
    pub fn invalidate_geometry_cache(&mut self) {
        self.sys.invalidate_geometry_cache();
    }

    /// Evaluate the shared objective once at unscaled radii.
    ///
    /// A thin forward to the free `evaluate_unscaled` below (crate-private,
    /// so this is code font and not a link) — one body, two
    /// spellings, so a caller holding a `&mut PackSystem` and a caller
    /// holding a `&mut PackState` cannot disagree about the verdict.
    pub fn evaluate_unscaled(&mut self, x: &[F]) -> (F, F, F) {
        evaluate_unscaled(&mut self.sys, x)
    }
}

/// Evaluate the shared objective once at **unscaled** radii, restoring every
/// field this function writes before it returns.
///
/// Returns `(f_total, fdist, frest)` from that evaluation — the triple the
/// GENCAN main loop feeds to `flast` / `fimp` / callback `StepReport`, and the
/// one both growth drivers turn into their `StageOutcome`. On return
/// `sys.fdist` / `sys.frest` / `sys.fdist_atom` / `sys.frest_atom` still
/// describe this unscaled evaluation — that radius-dependent inner state is
/// what a caller asks for — while `sys.radius`, `sys.scale` and `sys.scale2`
/// hold exactly the values the caller had on entry.
///
/// # Order of operations
///
/// ```text
/// save scale / scale2 → set the unscaled pair
///   save radius into work.radiuswork, swap in radius_ini
///     EvalMode::FOnly
///   restore radius from work.radiuswork
/// restore scale / scale2
/// ```
///
/// # Why both groups are restored
///
/// `scale` and `scale2` are `pub` fields of [`PackSystem`], which the crate
/// root re-exports, so a caller may legitimately hold values of its own in
/// them. The contract above is that this function gives the caller's values
/// back; restoring `radius` but not the scale pair would make that promise
/// true of only half the fields written here. Symmetry costs two stack slots
/// and removes an asymmetry from the contract.
///
/// # Why the save/restore is bitwise inert
///
/// Outside this function the crate writes `scale` / `scale2` in one place,
/// `pack::initial` (the port of Packmol `initial.f90:50-51`), and it writes
/// the `PackSystem::new` defaults (`1.0` and `DEFAULT_SCALE2` == `0.01`).
/// Setting the pair here therefore writes what every caller already holds,
/// and restoring it writes those same values back. The objective kernels
/// bind `scale2` into a local; they read it, they do not write it.
///
/// The radius swap is inert on those paths too: `radius` is scaled only by
/// the GENCAN schedule (`src/pack/gencan/phases.rs`) and by `movebad`, which
/// restores it, so growth evaluates at `radius == radius_ini` and the swap
/// moves the same numbers out and back. `work.radiuswork` is sized `ntotat`
/// by `WorkBuffers`, so it always has room for the whole array.
///
/// # What the swap moves, and what it does not
///
/// `fdist` accumulates the pairwise violation computed **unconditionally**
/// from `radius_ini` (the pair kernel's `rsum_ini`, outside the `overlap`
/// branch), and `frest` comes from the restraint terms, which read `scale` /
/// `scale2` but no radius. The radius swap therefore moves `f_total` alone —
/// which is why the regression golden can pin this triple and be pinning it
/// for the right reason.
///
/// It lives here because the pipeline layer may not import `gencan/`, and it
/// is `pub(crate)`: nothing outside the crate needs it.
pub(crate) fn evaluate_unscaled(sys: &mut PackSystem, x: &[F]) -> (F, F, F) {
    let scale = sys.scale;
    let scale2 = sys.scale2;
    sys.scale = 1.0;
    sys.scale2 = DEFAULT_SCALE2;

    sys.work.radiuswork.copy_from_slice(&sys.radius);
    for i in 0..sys.ntotat {
        sys.set_radius(i, sys.radius_ini[i]);
    }
    let f_total = sys.evaluate(x, EvalMode::FOnly, None).f_total;
    let fdist = sys.fdist;
    let frest = sys.frest;
    for i in 0..sys.ntotat {
        sys.set_radius(i, sys.work.radiuswork[i]);
    }

    sys.scale = scale;
    sys.scale2 = scale2;
    (f_total, fdist, frest)
}

#[cfg(test)]
mod tests;
