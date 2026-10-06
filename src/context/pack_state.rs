//! [`PackState`]: the state one packing run carries between stages.
//!
//! The type **wraps** a [`PackContext`] instead of extracting fields out of
//! it. Not one field moves, so the shared objective's read set is unchanged,
//! `ctx.xcart` stays the single lab-frame home of the coordinates, and
//! `ctx.fixedatom` / `ctx.comptype` are read through the wrapper rather than
//! mirrored — that mirror has drifted before, which is why
//! `PackContext::debug_assert_atom_props_sync` exists at all. On top of the
//! context the state adds the two things a chain of stages needs and the
//! context has no home for: the placement shape marker [`Placed`] and the
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
//! [`PackContext`] holds `Arc<dyn Restraint>` values, a cell grid and a
//! frame, and does not implement `Debug`, so `#[derive(Debug)]` on a struct
//! that owns one cannot compile. Rather than grow the context with a derive
//! it does not otherwise need, this module writes the impl by hand: it prints
//! `placed` and the view's molecule count and elides the context.

use std::fmt;

use molrs::op::types::F;

use crate::Objective;
use crate::context::DEFAULT_SCALE2;
use crate::context::{PackContext, RigidView};
use crate::eval::EvalMode;

/// The shape of the placements a [`PackState`] currently holds.
///
/// A state-level marker, not a per-atom set: it answers "has every free
/// molecule been placed yet?", which is what a chain of stages checks before
/// running the next one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placed {
    /// Nothing has been placed yet — a freshly wrapped context.
    None,
    /// Every free molecule the context lays out has a placement.
    All,
}

/// A [`PackContext`] plus the two pieces of run state a stage chain needs.
pub struct PackState {
    /// The wrapped context — the authority for geometry, radii and restraints.
    ctx: PackContext,
    /// The shape marker; see [`Placed`].
    placed: Placed,
    /// The rigid degrees of freedom, always present (never an `Option`).
    rigid: RigidView,
}

impl fmt::Debug for PackState {
    /// Prints the marker and the view's molecule count; the context is
    /// elided because it does not implement `Debug` (see the module docs).
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PackState")
            .field("placed", &self.placed)
            .field("nmol", &self.rigid.nmol())
            .finish_non_exhaustive()
    }
}

impl PackState {
    /// Wrap `ctx` into a fresh run state for `nmol` free molecules: a zeroed
    /// [`RigidView`] in the slot and [`Placed::None`] as the marker.
    pub fn new(ctx: PackContext, nmol: usize) -> Self {
        Self {
            ctx,
            placed: Placed::None,
            rigid: RigidView::fresh(nmol),
        }
    }

    /// The wrapped context.
    pub fn ctx(&self) -> &PackContext {
        &self.ctx
    }

    /// The wrapped context, mutably.
    pub fn ctx_mut(&mut self) -> &mut PackContext {
        &mut self.ctx
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

    /// Context and view as two disjoint mutable borrows, for the many
    /// operations that write placements while reading or updating the
    /// context (`write_xcart`, `capture_from_xcart`, a stage's inner loop).
    pub fn rigid_split_mut(&mut self) -> (&mut PackContext, &mut RigidView) {
        (&mut self.ctx, &mut self.rigid)
    }

    /// Give the two owned parts back, consuming the state.
    pub fn into_parts(self) -> (PackContext, RigidView) {
        (self.ctx, self.rigid)
    }

    /// Drop the context's cached Cartesian expansion.
    ///
    /// A forward to [`PackContext::invalidate_geometry_cache`], not a second
    /// implementation: the cache and its invalidation semantics belong to the
    /// context. It sits here because the caller is a stage boundary, which
    /// holds the state and not the bare context.
    pub fn invalidate_geometry_cache(&mut self) {
        self.ctx.invalidate_geometry_cache();
    }

    /// Evaluate the shared objective once at unscaled radii.
    ///
    /// A thin forward to the free `evaluate_unscaled` below (crate-private,
    /// so this is code font and not a link) — one body, two
    /// spellings, so a caller holding a `&mut PackContext` and a caller
    /// holding a `&mut PackState` cannot disagree about the verdict.
    pub fn evaluate_unscaled(&mut self, x: &[F]) -> (F, F, F) {
        evaluate_unscaled(&mut self.ctx, x)
    }
}

/// Evaluate the shared objective once at **unscaled** radii, restoring every
/// field this function writes before it returns.
///
/// Returns `(f_total, fdist, frest)` from that evaluation — the triple the
/// GENCAN main loop feeds to `flast` / `fimp` / handler `StepInfo`, and the
/// one both growth drivers turn into their `StageOutcome`. On return
/// `ctx.fdist` / `ctx.frest` / `ctx.fdist_atom` / `ctx.frest_atom` still
/// describe this unscaled evaluation — that radius-dependent inner state is
/// what a caller asks for — while `ctx.radius`, `ctx.scale` and `ctx.scale2`
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
/// `scale` and `scale2` are `pub` fields of [`PackContext`], which the crate
/// root re-exports, so a caller may legitimately hold values of its own in
/// them. The contract above is that this function gives the caller's values
/// back; restoring `radius` but not the scale pair would make that promise
/// true of only half the fields written here. Symmetry costs two stack slots
/// and removes an asymmetry from the contract.
///
/// # Why the merge is bitwise inert at today's callers
///
/// The crate writes `scale` / `scale2` in exactly three places, and all three
/// write the constructor defaults (`1.0` and `DEFAULT_SCALE2` == `0.01`, the
/// struct literal at `src/context/pack_context.rs:370-371`):
///
/// * `src/initial.rs:338-339`, the port of Packmol `initial.f90:50-51`;
/// * the final-verdict site of the continuum growth driver
///   (`src/grow/driver.rs`), now absorbed into this function;
/// * the same site in the lattice growth driver
///   (`src/grow/lattice/mod.rs`), likewise absorbed.
///
/// Setting the pair here therefore writes what those callers already held,
/// and restoring it writes those same values back: the save/restore pair is
/// inert at every call site that exists today, which is what lets the two
/// growth idioms collapse into this one without moving a bit.
/// (`src/objective.rs:414,622` bind `scale2` into a local; they read it, they
/// do not write it.)
///
/// The radius swap is inert on those paths too: `radius` is scaled only by
/// the GENCAN schedule (`src/gencan/phases.rs`) and by `movebad`, which
/// restores it, so growth evaluates at `radius == radius_ini` and the swap
/// moves the same numbers out and back. `work.radiuswork` is sized `ntotat`
/// by `WorkBuffers`, so it always has room for the whole array.
///
/// # What the swap moves, and what it does not
///
/// `fdist` accumulates the pairwise violation computed **unconditionally**
/// from `radius_ini` (`src/objective.rs:292-300`, outside the `overlap`
/// branch), and `frest` comes from the restraint terms, which read `scale` /
/// `scale2` but no radius. The radius swap therefore moves `f_total` alone —
/// which is why the regression golden can pin this triple and be pinning it
/// for the right reason.
///
/// # A published path is withdrawn here (booked)
///
/// This function used to be `pub fn evaluate_unscaled` in
/// `src/gencan/phases.rs`, reachable as
/// `molpack::gencan::phases::evaluate_unscaled`. It lives here because the
/// pipeline layer may not import `gencan/`, and it is `pub(crate)` because
/// this chain does not pay for a published symbol before its seam lands —
/// so the old path leaves the public surface. The docs site still places the
/// function in `phases.rs` (`docs/architecture.md:37`,
/// `docs/extending.md:460`); those two lines belong to the chain's
/// documentation task, not to this module.
pub(crate) fn evaluate_unscaled(ctx: &mut PackContext, x: &[F]) -> (F, F, F) {
    let scale = ctx.scale;
    let scale2 = ctx.scale2;
    ctx.scale = 1.0;
    ctx.scale2 = DEFAULT_SCALE2;

    ctx.work.radiuswork.copy_from_slice(&ctx.radius);
    for i in 0..ctx.ntotat {
        ctx.set_radius(i, ctx.radius_ini[i]);
    }
    let f_total = ctx.evaluate(x, EvalMode::FOnly, None).f_total;
    let fdist = ctx.fdist;
    let frest = ctx.frest;
    for i in 0..ctx.ntotat {
        ctx.set_radius(i, ctx.work.radiuswork[i]);
    }

    ctx.scale = scale;
    ctx.scale2 = scale2;
    (f_total, fdist, frest)
}

#[cfg(test)]
mod tests;
