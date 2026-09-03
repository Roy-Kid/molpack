//! The rigid degrees of freedom of one packing run, and nothing else.
//!
//! [`RigidView`] owns the flat placement vector every solver optimizes: for
//! each free molecule a centre of mass and three Euler angles. The layout is
//! the Packmol convention the whole crate already speaks —
//!
//! ```text
//! x = [ com(0) com(1) … com(n-1) | euler(0) euler(1) … euler(n-1) ]
//!       3 values each              3 values each
//! ```
//!
//! COM block first, Euler block second, three values per molecule each, so
//! the vector is `6 * nmol` long. A molecule's index into this vector is its
//! index in the context layout (`idfirst` / `natoms` / `nmols`, type-major
//! then copy-major), which is also `xcart`'s and `coor`'s molecule order:
//! one index space, no translation table.
//!
//! Two coordinate frames meet here, both stored on a [`PackContext`].
//! `ctx.coor` holds each copy's *reference conformer* — its atoms' positions
//! measured from that copy's own centre of mass, with no orientation applied,
//! so it describes the molecule's shape and nothing else. `ctx.xcart` holds
//! *lab-frame* coordinates: where those atoms actually sit in the packing
//! cell. The three Euler angles are the three-parameter description of the
//! rotation matrix `R` that carries the first frame into the second (see
//! [`crate::euler`]), so a molecule's placement is exactly `(com, euler)` and
//! its shape is exactly its block of `coor`.
//!
//! Besides the accessors the view owns the two conversions between those
//! frames:
//!
//! - [`RigidView::write_xcart`] expands `(com, euler)` plus the per-copy
//!   reference conformer into lab-frame coordinates
//!   (`xcart = com + R(euler) · coor`) — the outbound rebuild;
//! - [`RigidView::capture_from_xcart`] reads lab-frame coordinates back into
//!   `(com, euler = 0)` plus a centred conformer in `coor` — the growth
//!   writeback.
//!
//! [`RigidView::install_seed`] is the third and last place the view touches a
//! [`PackContext`], and the only constructing one: it takes a placement
//! solution from a previous run as two bare slices, so this module never
//! names the `entry` layer that produced them (dependencies run
//! `entry → context`, never back).
//!
//! # What this view deliberately does not carry
//!
//! There is **no flag recording that these placements were fed in from
//! outside**. That fact is state, not geometry, and its home is
//! `PackState.placed`. Storing it here as well would create a second
//! representation that nothing maintains: the GENCAN phases write solutions
//! through [`RigidView::set_com`] and [`RigidView::as_mut_slice`] and would
//! never update such a flag, so the two copies would diverge on the first
//! iteration. The view knows the `6 * nmol` numbers it holds, and only those.
//!
//! # Why the accessors bounds-check explicitly
//!
//! Offset arithmetic alone does not catch an out-of-range molecule. With
//! `i == nmol`, `3 * i` lands on the first slot of the *Euler* block — a
//! silent cross-block write into molecule 0 — while the Euler accessors trap
//! only because they run off the end of the buffer. The asymmetry is why
//! every index accessor asserts `i < nmol` up front: an out-of-range molecule
//! is a programming error and must be unrepresentable rather than quietly
//! corrupting another molecule's placement.

use crate::context::PackContext;
use crate::euler::{compcart, eulerrmat};
use molrs::types::F;

/// The rigid placement vector of one run: three COM plus three Euler values
/// per free molecule, in the flat layout documented at module level.
///
/// The view owns its buffer, so its length and its molecule count can never
/// disagree. Construct it with [`RigidView::fresh`] (a zeroed run) or with
/// [`RigidView::install_seed`] (continuing from another run's solution).
#[derive(Debug, Clone)]
pub struct RigidView {
    /// `6 * nmol` variables: COM block first, then the Euler block.
    x: Vec<F>,
    /// Number of molecules addressed — `x.len() / 6`, kept as the layout's
    /// stride so the Euler block's base is a single field read.
    nmol: usize,
}

impl RigidView {
    /// A zeroed view for `nmol` molecules: every COM at the origin, every
    /// Euler angle zero.
    pub fn fresh(nmol: usize) -> Self {
        Self {
            x: vec![0.0 as F; 6 * nmol],
            nmol,
        }
    }

    /// Number of molecules this view addresses.
    pub fn nmol(&self) -> usize {
        self.nmol
    }

    /// Centre of mass of molecule `i`.
    ///
    /// # Panics
    ///
    /// Panics when `i >= nmol` (see the module docs on bounds checking).
    pub fn com(&self, i: usize) -> [F; 3] {
        assert!(i < self.nmol, "molecule {i} is outside 0..{}", self.nmol);
        let o = 3 * i;
        [self.x[o], self.x[o + 1], self.x[o + 2]]
    }

    /// Set the centre of mass of molecule `i`.
    ///
    /// # Panics
    ///
    /// Panics when `i >= nmol` — without this check the write would land in
    /// the Euler block of another molecule instead of failing.
    pub fn set_com(&mut self, i: usize, com: [F; 3]) {
        assert!(i < self.nmol, "molecule {i} is outside 0..{}", self.nmol);
        let o = 3 * i;
        self.x[o..o + 3].copy_from_slice(&com);
    }

    /// Euler angles `(beta, gamma, theta)` of molecule `i`.
    ///
    /// # Panics
    ///
    /// Panics when `i >= nmol`.
    pub fn euler(&self, i: usize) -> [F; 3] {
        assert!(i < self.nmol, "molecule {i} is outside 0..{}", self.nmol);
        let o = 3 * self.nmol + 3 * i;
        [self.x[o], self.x[o + 1], self.x[o + 2]]
    }

    /// Set the Euler angles of molecule `i`.
    ///
    /// # Panics
    ///
    /// Panics when `i >= nmol`.
    pub fn set_euler(&mut self, i: usize, euler: [F; 3]) {
        assert!(i < self.nmol, "molecule {i} is outside 0..{}", self.nmol);
        let o = 3 * self.nmol + 3 * i;
        self.x[o..o + 3].copy_from_slice(&euler);
    }

    /// The backing flat vector, for handing to the shared objective.
    pub fn as_slice(&self) -> &[F] {
        &self.x
    }

    /// Mutable access to the backing flat vector. The GENCAN path drives its
    /// phase machinery over the raw layout; growth keeps to the typed
    /// accessors above.
    pub fn as_mut_slice(&mut self) -> &mut [F] {
        &mut self.x
    }

    /// Expand the rigid DOF into lab-frame coordinates:
    /// `xcart = com + R(euler) · coor` for every atom of every copy.
    ///
    /// Molecules are enumerated in the context's own layout order, and
    /// `coor` shares `xcart`'s index space, so each copy is rebuilt from its
    /// own reference conformer (copies diverge once an in-loop optimizer has
    /// relaxed them independently).
    ///
    /// This is the one rebuild in the crate — the outbound direction at the
    /// end of a run, the inbound one that materializes a seeded state before
    /// push-off, and the initialization passes all call it; no caller
    /// re-derives `xcart` from the placements on its own.
    ///
    /// # Panics
    ///
    /// Panics when the context lays out more free molecules than the view
    /// holds (`sum(ctx.nmols[..ctx.ntype]) > nmol`): view and context must
    /// have been sized from the same run. Debug builds additionally assert
    /// that every atom rebuilt here is a free one — a fixed structure's
    /// coordinates are given, not placed, and must never be overwritten from
    /// `(com, euler)`.
    pub fn write_xcart(&self, ctx: &mut PackContext) {
        let mut imol = 0usize;
        let mut icart = 0usize;

        for itype in 0..ctx.ntype {
            for _ in 0..ctx.nmols[itype] {
                let xcm = self.com(imol);
                let [beta, gama, teta] = self.euler(imol);
                let (v1, v2, v3) = eulerrmat(beta, gama, teta);

                for _ in 0..ctx.natoms[itype] {
                    let pos = compcart(&xcm, &ctx.coor[icart], &v1, &v2, &v3);
                    ctx.xcart[icart] = pos;
                    // Packmol's initial.f90 sets fixedatom=false on every
                    // free atom here, but in Rust that bit is already false
                    // from construction and `sync_atom_props` has been
                    // called — writing it again would desync `atom_props`.
                    debug_assert!(!ctx.fixedatom[icart]);
                    icart += 1;
                }

                imol += 1;
            }
        }
    }

    /// Build a view from a previous run's placement solution and install that
    /// run's conformers into `ctx`.
    ///
    /// The solution arrives as plain data — the flat placement vector and the
    /// per-copy reference conformers — so this layer never depends on the
    /// snapshot type the entry layer hands out; the caller unpacks it. Both
    /// halves are copied verbatim (bitwise), which is what makes chaining one
    /// run into the next reproduce the first run's geometry exactly. Atoms
    /// past `coor.len()` (the fixed structures) are left untouched.
    ///
    /// The molecule count comes from `ctx.ntotmol`, the number of free
    /// molecules the context was built for — the same count that sizes the
    /// run's placement vector.
    ///
    /// # Panics
    ///
    /// Panics when `x.len() != 6 * ctx.ntotmol`, or when `coor` is longer
    /// than the context's coordinate array.
    pub fn install_seed(x: &[F], coor: &[[F; 3]], ctx: &mut PackContext) -> Self {
        let nmol = ctx.ntotmol;
        assert_eq!(
            x.len(),
            6 * nmol,
            "placement vector holds 6 variables per molecule (3 COM + 3 Euler)"
        );
        ctx.coor[..coor.len()].copy_from_slice(coor);
        Self {
            x: x.to_vec(),
            nmol,
        }
    }

    /// Capture lab-frame coordinates back into the rigid DOF — the writeback
    /// contract the growth drivers produce.
    ///
    /// Per copy: the COM is the centroid of that copy's `xcart` atoms, the
    /// centred conformer (`xcart − com`) goes into that copy's own `coor`
    /// block, and the Euler angles are zeroed. Growth builds a shape rather
    /// than rotating a template, so the whole conformation lives in `coor`
    /// and the rotation is the identity; a stale rotation from an earlier
    /// stage must therefore be cleared, not preserved.
    ///
    /// Both halves belong to one call: writing the centred conformer is the
    /// exact inverse of [`RigidView::write_xcart`], and splitting them would
    /// leave the caller a step it could forget, with `coor` and the COM
    /// describing different geometry.
    ///
    /// `xcart` itself is read, never written: it is the authoritative
    /// lab-frame state at this point, and every driver syncs its own buffers
    /// into it before capturing.
    ///
    /// # Panics
    ///
    /// Panics when the context lays out more free molecules than the view
    /// holds (`sum(ctx.nmols[..ctx.ntype]) > nmol`), the same sizing contract
    /// as [`RigidView::write_xcart`].
    pub fn capture_from_xcart(&mut self, ctx: &mut PackContext) {
        let mut imol = 0usize;

        for itype in 0..ctx.ntype {
            let na = ctx.natoms[itype];
            for icopy in 0..ctx.nmols[itype] {
                let base = ctx.idfirst[itype] + icopy * na;

                let mut com = [0.0 as F; 3];
                for a in 0..na {
                    let p = ctx.xcart[base + a];
                    for k in 0..3 {
                        com[k] += p[k];
                    }
                }
                for v in com.iter_mut() {
                    *v /= na as F;
                }

                for a in 0..na {
                    let p = ctx.xcart[base + a];
                    ctx.coor[base + a] = [p[0] - com[0], p[1] - com[1], p[2] - com[2]];
                }

                self.set_com(imol, com);
                self.set_euler(imol, [0.0, 0.0, 0.0]);
                imol += 1;
            }
        }
    }
}
