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
//! rotation matrix `R` that carries the first frame into the second (the
//! `euler` module), so a molecule's placement is exactly `(com, euler)` and
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
use molrs::op::F;
use molrs::op::centroid;
use molrs::op::vec3::sub;

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
            let unit_weights = vec![1.0; na];
            for icopy in 0..ctx.nmols[itype] {
                let base = ctx.idfirst[itype] + icopy * na;
                let atoms = base..base + na;

                // At unit weights the sum runs in atom order and is divided
                // by `na` once: bit for bit the plain mean.
                let com = centroid(&ctx.xcart[atoms.clone()], &unit_weights)
                    .expect("a molecule type has atoms");

                for i in atoms {
                    ctx.coor[i] = sub(ctx.xcart[i], com);
                }

                self.set_com(imol, com);
                self.set_euler(imol, [0.0, 0.0, 0.0]);
                imol += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    //! Contract tests for `src/context/rigid_view.rs`.
    //!
    //! `RigidView` is the single home of the rigid degrees of freedom: the flat
    //! `6 * nmol` placement vector (COM block first, Euler block second, three
    //! values per molecule each) plus the three operations that cross between it
    //! and a `PackContext` —
    //!
    //! - `write_xcart` — `xcart = com + R(euler) · coor` (the outbound rebuild,
    //!   inherited verbatim from `initial::init_xcart_from_x`),
    //! - `install_seed` — inject a seed placement + conformer as plain data,
    //! - `capture_from_xcart` — read lab-frame coordinates back into
    //!   `(com, euler = 0)` plus a centered conformer in `ctx.coor` (the growth
    //!   writeback contract).
    //!
    //! Fixtures build a `PackContext` directly (`PackContext::new` + the public
    //! layout fields), the same way `pack_context.rs::geometry_cache_tests` does — no engine,
    //! no solver, no growth driver. Categories: basics, edge cases, immutability
    //! (a fresh view is all zeros; a clone is an independent snapshot). No
    //! physics is asserted here: the view is pure bookkeeping over a
    //! coordinate layout.

    use crate::context::RigidView as ContextRigidView;
    use crate::euler::{compcart, eulerrmat};
    use crate::{PackContext, RigidView};
    use molrs::op::F;

    /// The crate-root re-export and the `context` path name one type, not two.
    /// Compile-time only.
    fn _both_paths_are_one_type(v: RigidView) -> ContextRigidView {
        v
    }

    // ── fixtures ───────────────────────────────────────────────────────────────

    /// Two three-atom copies of one type: `ntotat = 6`, `ntotmol = 2`, layout
    /// `idfirst = [0]`. `coor` holds one reference conformer **per copy**,
    /// sharing `xcart`'s index space (the `PackContext` convention).
    fn two_copies_of_three() -> PackContext {
        let mut ctx = PackContext::new(6, 2, 1);
        ctx.nmols = vec![2];
        ctx.natoms = vec![3];
        ctx.idfirst = vec![0];
        ctx.coor = vec![[0.0; 3]; 6];
        ctx
    }

    /// The regression conformer: two copies whose reference blocks are already
    /// centred (each block sums to exactly zero), so the COM survives the
    /// `write_xcart` → `capture_from_xcart` round trip.
    const CENTERED_COOR: [[F; 3]; 6] = [
        [1.0, 0.0, -0.5],
        [-0.5, 0.25, 0.75],
        [-0.5, -0.25, -0.25],
        [0.5, -1.5, 0.25],
        [-1.25, 0.75, -1.0],
        [0.75, 0.75, 0.75],
    ];

    /// Centroid of `coords` computed with the growth writeback's own arithmetic
    /// (`src/grow/driver.rs`): accumulate component-wise in atom order, then
    /// divide once by the atom count. Same order, same bits.
    fn centroid(coords: &[[F; 3]]) -> [F; 3] {
        let mut com = [0.0 as F; 3];
        for p in coords {
            for k in 0..3 {
                com[k] += p[k];
            }
        }
        for v in com.iter_mut() {
            *v /= coords.len() as F;
        }
        com
    }

    // ── Category: basics — the flat-vector layout contract ─────────────────────

    /// `fresh(nmol)` replaces `vec![0.0; 6 * nmol]`: the view owns its buffer,
    /// it is `6 * nmol` long, and every degree of freedom starts at zero.
    #[test]
    fn rigid_view_fresh_is_zeroed_and_sized() {
        let view = RigidView::fresh(4);
        assert_eq!(view.nmol(), 4, "fresh(4) addresses 4 molecules");
        assert_eq!(
            view.as_slice().len(),
            6 * 4,
            "the flat vector holds 6 variables per molecule (3 COM + 3 Euler)"
        );
        assert!(
            view.as_slice().iter().all(|&v| v == 0.0),
            "a fresh view is all zeros: {:?}",
            view.as_slice()
        );
        for i in 0..4 {
            assert_eq!(view.com(i), [0.0; 3], "fresh com({i})");
            assert_eq!(view.euler(i), [0.0; 3], "fresh euler({i})");
        }
    }

    /// The view writes the Packmol flat-vector convention: COM of molecule `i` at
    /// `x[3*i .. 3*i+3]`, Euler of molecule `i` at
    /// `x[3*nmol + 3*i .. 3*nmol + 3*i + 3]`.
    ///
    /// Ported from `grow::tests::placements_view_layout`, which pinned the same
    /// contract on the deleted `PlacementsMut`.
    #[test]
    fn rigid_view_layout() {
        let nmol = 3;
        let mut view = RigidView::fresh(nmol);
        assert_eq!(view.nmol(), nmol);
        view.set_com(1, [1.0, 2.0, 3.0]);
        view.set_euler(2, [0.1, 0.2, 0.3]);
        // Read back through the typed accessors.
        assert_eq!(view.com(1), [1.0, 2.0, 3.0]);
        assert_eq!(view.euler(2), [0.1, 0.2, 0.3]);
        // Untouched copies read zero.
        assert_eq!(view.com(0), [0.0; 3]);
        assert_eq!(view.euler(0), [0.0; 3]);

        // Raw slots: COM block first, Euler block at 3*nmol.
        let x = view.as_slice();
        assert_eq!(&x[3..6], &[1.0, 2.0, 3.0], "com(1) slot");
        assert_eq!(
            &x[15..18],
            &[0.1, 0.2, 0.3],
            "euler(2) slot = x[3*3 + 3*2 ..]"
        );
        // Every other slot untouched.
        for (i, &v) in x.iter().enumerate() {
            if !(3..6).contains(&i) && !(15..18).contains(&i) {
                assert_eq!(v, 0.0, "slot {i} must be untouched");
            }
        }
    }

    /// `as_mut_slice` is the raw door the GENCAN phases drive; writes through it
    /// are visible to the typed accessors and vice versa — one buffer, not two.
    #[test]
    fn rigid_view_as_mut_slice_shares_one_buffer() {
        let mut view = RigidView::fresh(2);
        view.as_mut_slice()[4] = 7.5; // com(1).y
        view.as_mut_slice()[6] = -1.25; // euler(0).beta
        assert_eq!(view.com(1), [0.0, 7.5, 0.0]);
        assert_eq!(view.euler(0), [-1.25, 0.0, 0.0]);

        view.set_com(0, [1.0, 2.0, 3.0]);
        assert_eq!(&view.as_slice()[0..3], &[1.0, 2.0, 3.0]);
        assert_eq!(view.as_mut_slice().len(), 12);
    }

    // ── Category: basics — write_xcart (com + R(euler) · coor) ─────────────────

    /// `write_xcart` expands the rigid DOF into lab-frame coordinates with the
    /// crate's own Euler convention: for every copy, `xcart = com + R · coor`,
    /// enumerated over the `(idfirst, natoms, nmols)` layout.
    ///
    /// The expectation is recomputed here from `euler::eulerrmat` + `compcart` —
    /// the same arithmetic in the same order, so equality is exact.
    #[test]
    fn rigid_view_write_xcart_composes_com_and_rotation() {
        let mut ctx = two_copies_of_three();
        ctx.coor = CENTERED_COOR.to_vec();

        let mut view = RigidView::fresh(2);
        let coms = [[10.0, -3.0, 2.5], [-4.0, 7.25, 0.5]];
        let eulers = [[0.3, -0.7, 1.1], [2.0, 0.5, -1.25]];
        for i in 0..2 {
            view.set_com(i, coms[i]);
            view.set_euler(i, eulers[i]);
        }

        view.write_xcart(&mut ctx);

        for imol in 0..2 {
            let (v1, v2, v3) = eulerrmat(eulers[imol][0], eulers[imol][1], eulers[imol][2]);
            for a in 0..3 {
                let icart = 3 * imol + a;
                let want = compcart(&coms[imol], &CENTERED_COOR[icart], &v1, &v2, &v3);
                for (k, (got, want)) in ctx.xcart[icart].iter().zip(&want).enumerate() {
                    assert_eq!(
                        got, want,
                        "atom {icart} component {k}: xcart must be com + R(euler)·coor"
                    );
                }
            }
        }
    }

    /// Zero Euler angles are the identity rotation, so `write_xcart` degenerates
    /// to a pure translation of the stored conformer. This is the exact shape the
    /// growth path relies on: `capture_from_xcart` leaves `euler = 0`, and the
    /// outbound rebuild must then reproduce the captured coordinates.
    #[test]
    fn rigid_view_write_xcart_identity_rotation_is_translation() {
        let mut ctx = two_copies_of_three();
        ctx.coor = CENTERED_COOR.to_vec();

        let mut view = RigidView::fresh(2);
        view.set_com(0, [1.0, 2.0, 4.0]);
        view.set_com(1, [-8.0, 0.5, 16.0]);

        view.write_xcart(&mut ctx);

        for imol in 0..2 {
            let com = view.com(imol);
            for a in 0..3 {
                let icart = 3 * imol + a;
                let want = [
                    com[0] + CENTERED_COOR[icart][0],
                    com[1] + CENTERED_COOR[icart][1],
                    com[2] + CENTERED_COOR[icart][2],
                ];
                for (k, (got, want)) in ctx.xcart[icart].iter().zip(&want).enumerate() {
                    assert_eq!(
                        got, want,
                        "atom {icart} component {k}: euler = 0 must be a pure translation"
                    );
                }
            }
        }
    }

    // ── Category: basics — install_seed (plain-data injection) ─────────────────

    /// `install_seed` is the constructing form: the seed arrives as two bare
    /// slices (never as an `entry`-layer type), the conformer is copied into the
    /// leading `coor` block bitwise, and the returned view carries the seed
    /// placements verbatim. Atoms past the seed's length (fixed targets) keep
    /// whatever `coor` held.
    #[test]
    fn rigid_view_install_seed_copies_conformer_and_placements() {
        let mut ctx = PackContext::new(8, 2, 1);
        ctx.nmols = vec![2];
        ctx.natoms = vec![3];
        ctx.idfirst = vec![0];
        // Six free-atom slots plus two trailing slots that must not be touched.
        ctx.coor = vec![[9.0, 9.0, 9.0]; 8];

        let seed_coor: [[F; 3]; 6] = CENTERED_COOR;
        let seed_x: [F; 12] = [
            1.5, -2.5, 3.5, // com of molecule 0
            -0.25, 0.75, 8.0, // com of molecule 1
            0.1, 0.2, 0.3, // euler of molecule 0
            -0.4, 0.5, -0.6, // euler of molecule 1
        ];

        let view = RigidView::install_seed(&seed_x, &seed_coor, &mut ctx);

        assert_eq!(view.nmol(), 2, "the seed covers both free copies");
        assert_eq!(
            view.as_slice(),
            &seed_x[..],
            "the seed placements are injected verbatim (zero-conversion chaining)"
        );
        for (i, want) in seed_coor.iter().enumerate() {
            for (k, (got, want)) in ctx.coor[i].iter().zip(want).enumerate() {
                assert_eq!(
                    got.to_bits(),
                    want.to_bits(),
                    "coor[{i}] component {k} must be a bitwise copy of the seed"
                );
            }
        }
        for i in 6..8 {
            assert_eq!(
                ctx.coor[i],
                [9.0, 9.0, 9.0],
                "coor[{i}] is past the seed and must be left alone"
            );
        }
        // The typed accessors read the injected layout, not a re-derivation.
        assert_eq!(view.com(1), [-0.25, 0.75, 8.0]);
        assert_eq!(view.euler(0), [0.1, 0.2, 0.3]);
    }

    // ── Category: basics — capture_from_xcart (the growth writeback) ───────────

    /// The writeback contract, lifted from the growth drivers: per copy the COM
    /// is the centroid of that copy's lab-frame atoms, `ctx.coor` receives the
    /// centered conformer (`xcart − com`), and the Euler angles are zeroed —
    /// growth produces no rigid rotation, the whole shape lives in `coor`.
    ///
    /// Coordinates are deliberately non-representable in binary, so the assertion
    /// also pins the summation order (accumulate in atom order, divide once).
    #[test]
    fn rigid_view_capture_from_xcart_centers_each_copy() {
        let mut ctx = two_copies_of_three();
        let xcart: [[F; 3]; 6] = [
            [0.1, 0.2, 0.3],
            [1.3, -0.7, 2.9],
            [-0.4, 3.1, 0.55],
            [7.7, 1.1, -2.2],
            [8.3, 0.9, -1.05],
            [9.15, 2.35, -3.4],
        ];
        ctx.xcart = xcart.to_vec();
        // Pre-fill `coor` with a sentinel so the writeback is visible.
        ctx.coor = vec![[42.0; 3]; 6];

        let mut view = RigidView::fresh(2);
        // A stale rotation must be cleared, not preserved.
        view.set_euler(0, [1.0, 2.0, 3.0]);
        view.set_euler(1, [-1.0, -2.0, -3.0]);

        view.capture_from_xcart(&mut ctx);

        for imol in 0..2 {
            let block = &xcart[3 * imol..3 * imol + 3];
            let com = centroid(block);
            let got_com = view.com(imol);
            for (k, (got, want)) in got_com.iter().zip(&com).enumerate() {
                assert_eq!(
                    got, want,
                    "molecule {imol} COM component {k} must be the centroid of its atoms"
                );
            }
            assert_eq!(
                view.euler(imol),
                [0.0; 3],
                "molecule {imol}: growth writes no rigid rotation, so euler is zeroed"
            );
            for a in 0..3 {
                let icart = 3 * imol + a;
                for k in 0..3 {
                    assert_eq!(
                        ctx.coor[icart][k],
                        xcart[icart][k] - com[k],
                        "coor[{icart}] component {k} must be the centered conformer"
                    );
                }
            }
        }
        assert_eq!(
            ctx.xcart,
            xcart.to_vec(),
            "capture reads xcart; it must not rewrite it"
        );
    }

    // ── Category: edge cases ───────────────────────────────────────────────────

    /// `fresh(0)`: an empty view is legal (a run with no free molecules), it has
    /// no slots, and both context crossings are no-ops.
    #[test]
    fn rigid_view_fresh_zero_molecules() {
        let mut view = RigidView::fresh(0);
        assert_eq!(view.nmol(), 0);
        assert!(view.as_slice().is_empty(), "fresh(0) has no slots");

        let mut ctx = PackContext::new(0, 0, 0);
        ctx.nmols = Vec::new();
        ctx.natoms = Vec::new();
        ctx.idfirst = Vec::new();
        ctx.coor = Vec::new();

        view.write_xcart(&mut ctx);
        assert!(ctx.xcart.is_empty(), "nothing to expand");
        view.capture_from_xcart(&mut ctx);
        assert!(ctx.coor.is_empty(), "nothing to capture");
        assert_eq!(view.nmol(), 0);
    }

    /// A molecule index outside `0..nmol` is a programming error, not a runtime
    /// condition, and the view must say so.
    ///
    /// Ported from `grow::tests::placements_view_rejects_bad_len`: the deleted
    /// `PlacementsMut` could be handed a mis-sized backing slice, while
    /// `RigidView::fresh` owns its buffer (illegal state unrepresentable), so the
    /// only reachable length error left is an out-of-range molecule.
    ///
    /// This needs an explicit `i < nmol` check on the accessors. The offset
    /// arithmetic alone does NOT catch it: `set_com(nmol, ..)` lands on
    /// `x[3*nmol .. 3*nmol+3]`, which is molecule 0's **Euler** slot — a
    /// silent cross-block write, the worst possible outcome (`set_euler` out of
    /// range does trap, since it runs off the end of the buffer; that asymmetry
    /// is exactly why the check has to be explicit).
    #[test]
    #[should_panic]
    fn rigid_view_set_com_out_of_range_panics() {
        let mut view = RigidView::fresh(2);
        view.set_com(2, [1.0, 2.0, 3.0]);
    }

    /// The mirror of the COM case on the Euler block: an out-of-range molecule is
    /// refused, never wrapped around into another molecule's slot.
    #[test]
    #[should_panic]
    fn rigid_view_set_euler_out_of_range_panics() {
        let mut view = RigidView::fresh(2);
        view.set_euler(2, [0.1, 0.2, 0.3]);
    }

    /// A one-atom copy has its COM exactly on the atom and a conformer of exactly
    /// zero — no drift from the centroid division.
    #[test]
    fn rigid_view_capture_from_xcart_single_atom_copy() {
        let mut ctx = PackContext::new(2, 2, 1);
        ctx.nmols = vec![2];
        ctx.natoms = vec![1];
        ctx.idfirst = vec![0];
        ctx.coor = vec![[7.0; 3]; 2];
        ctx.xcart = vec![[1.25, -3.5, 0.125], [-9.0, 0.0, 4.75]];

        let mut view = RigidView::fresh(2);
        view.capture_from_xcart(&mut ctx);

        assert_eq!(view.com(0), [1.25, -3.5, 0.125]);
        assert_eq!(view.com(1), [-9.0, 0.0, 4.75]);
        assert_eq!(ctx.coor[0], [0.0; 3], "a single atom sits on its own COM");
        assert_eq!(ctx.coor[1], [0.0; 3]);
        assert_eq!(view.euler(0), [0.0; 3]);
        assert_eq!(view.euler(1), [0.0; 3]);
    }

    /// With several types of different sizes, the view enumerates the context's
    /// `(idfirst, natoms, nmols)` layout, so a molecule index in `x` is exactly
    /// the molecule's index in the `xcart` layout — type-major, copy-major.
    #[test]
    fn rigid_view_capture_from_xcart_multitype_layout() {
        // type 0: two copies of 2 atoms (icart 0..4); type 1: one copy of 3
        // atoms (icart 4..7).
        let mut ctx = PackContext::new(7, 3, 2);
        ctx.nmols = vec![2, 1];
        ctx.natoms = vec![2, 3];
        ctx.idfirst = vec![0, 4];
        ctx.coor = vec![[0.0; 3]; 7];
        ctx.xcart = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0], // molecule 0 → centroid [1, 0, 0]
            [0.0, 10.0, 0.0],
            [0.0, 14.0, 0.0], // molecule 1 → centroid [0, 12, 0]
            [0.0, 0.0, 30.0],
            [3.0, 0.0, 30.0],
            [0.0, 3.0, 30.0], // molecule 2 → centroid [1, 1, 30]
        ];

        let mut view = RigidView::fresh(3);
        view.capture_from_xcart(&mut ctx);

        assert_eq!(view.com(0), [1.0, 0.0, 0.0], "molecule 0 = type 0, copy 0");
        assert_eq!(view.com(1), [0.0, 12.0, 0.0], "molecule 1 = type 0, copy 1");
        assert_eq!(view.com(2), [1.0, 1.0, 30.0], "molecule 2 = type 1, copy 0");
        // The centered conformer lands in each copy's OWN `coor` block.
        assert_eq!(ctx.coor[2], [0.0, -2.0, 0.0], "copy 1's first atom");
        assert_eq!(ctx.coor[3], [0.0, 2.0, 0.0], "copy 1's second atom");
        assert_eq!(ctx.coor[4], [-1.0, -1.0, 0.0], "type 1's first atom");
        for i in 0..3 {
            assert_eq!(view.euler(i), [0.0; 3], "molecule {i} euler");
        }
    }

    // ── Category: immutability ─────────────────────────────────────────────────

    /// The view owns its buffer, so a clone is an independent snapshot: writing
    /// through the clone leaves the original untouched. (`Debug` is part of the
    /// pinned surface too.)
    #[test]
    fn rigid_view_clone_is_independent_snapshot() {
        let mut original = RigidView::fresh(2);
        original.set_com(0, [1.0, 2.0, 3.0]);
        original.set_euler(1, [0.4, 0.5, 0.6]);

        let mut clone = original.clone();
        assert_eq!(
            clone.as_slice(),
            original.as_slice(),
            "a clone starts equal to its source"
        );

        clone.set_com(0, [-1.0, -2.0, -3.0]);
        assert_eq!(
            original.com(0),
            [1.0, 2.0, 3.0],
            "writing through the clone must not reach the original"
        );
        assert_eq!(clone.com(0), [-1.0, -2.0, -3.0]);
        assert_eq!(
            original.euler(1),
            [0.4, 0.5, 0.6],
            "the untouched half of the original is unchanged"
        );
        assert!(!format!("{original:?}").is_empty(), "RigidView is Debug");
    }
}
