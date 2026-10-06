//! Pairwise separation: a lower bound on the distance between two copies of the
//! same species.
//!
//! The packer's own pair term keeps *atoms* from overlapping, and it stops
//! caring the moment two molecules are no longer touching. Nothing in it
//! distinguishes two copies of one species from a copy of each of two — so a
//! species is free to pile its copies into one corner as long as they do not
//! interpenetrate. This term is the missing statement: *these molecules also
//! keep their distance from each other*.

use molrs::spatial::neighbors::CellGrid;
use molrs::spatial::simbox::{Mic, SimBox};
use molrs::types::F;

use super::com;
use super::{GroupCtx, Restraint};

/// Keep every pair of copies of one species at least `d_min` apart, measured
/// **centre to centre**.
///
/// The penalty is silent above `d_min` and grows below it:
///
/// ```text
///   E = λ · ∑_{c<c'}  (d_min − D)²        for D < d_min,   D = ‖mic(R_c − R_c')‖
/// ```
///
/// where `R_c` is copy `c`'s geometric centroid. It is **quadratic in the
/// length by which the bound is missed** — the same shape as a geometric
/// restraint's penalty (`InsideBoxRestraint` is `scale·(overshoot)²`), and
/// deliberately not the pair term's quartic-in-length form. This restraint
/// reports into `frest`, so it has to be commensurate with the other things
/// there: `precision = 0.01` then means "within 0.1 Å" for a separation bound
/// exactly as it does for a box wall, and that reading does not drift with
/// `d_min`. The gradient is continuous at contact, where the penalty and its
/// slope both vanish.
///
/// Distances use the minimum-image convention from [`GroupCtx`], so copies on
/// opposite sides of a periodic boundary are pushed apart exactly as neighbours
/// in the interior are.
///
/// Two exactly coincident centroids have no defined direction to separate
/// along and so receive no force — the same measure-zero convention the radial
/// geometries use. The pair term keeps that case from arising in practice.
///
/// # Centre to centre, not atom to atom
///
/// For a compact molecule the centroid is a faithful stand-in for the whole.
/// For a long or branched one it is not: two chains can interdigitate with
/// distant centroids, or share a centroid while their atoms stay apart. This
/// term states what it measures; it does not claim to bound atom–atom contact.
///
/// # Feasibility
///
/// Nothing here checks that `count` copies at `d_min` actually fit in the
/// region. An impossible request simply does not converge, and the run reports
/// it through `frest` like any other unsatisfied restraint.
///
/// # Example
///
/// ```
/// # use molpack::Target;
/// # use molpack::restraint::SelfSeparation;
/// // 50 ions that must stay 12 Å apart from one another.
/// let ions = Target::from_coords(&[[0.0; 3]], &[1.5], 50)
///     .with_collective_restraint(SelfSeparation::new(12.0, 1.0));
/// ```
#[derive(Debug, Clone)]
pub struct SelfSeparation {
    d_min: F,
    strength: F,
}

impl SelfSeparation {
    /// - `d_min` — minimum centre-to-centre distance between two copies.
    /// - `strength` — overall multiplier `λ`; `1.0` weights a shortfall like a
    ///   geometric restraint weights an equal overshoot, below `1.0` makes the
    ///   bound softer.
    ///
    /// # Panics
    /// If `d_min` or `strength` is not positive.
    pub fn new(d_min: F, strength: F) -> Self {
        assert!(d_min > 0.0, "SelfSeparation d_min must be positive");
        assert!(strength > 0.0, "SelfSeparation strength must be positive");
        Self { d_min, strength }
    }

    /// The minimum centre-to-centre distance this restraint asks for.
    pub fn d_min(&self) -> F {
        self.d_min
    }
}

impl SelfSeparation {
    /// Visit every unordered pair of sites closer than `d_min`, exactly once.
    ///
    /// One home for "which pairs are close", so the value-only and fused paths
    /// cannot drift apart on which pairs they consider or which image they
    /// measure.
    ///
    /// # Why a partition and not a double loop
    ///
    /// The bound is local — it is silent above `d_min` — but a double loop is
    /// not: it pays for every pair whether or not the two are anywhere near
    /// each other, so it costs the same on a converged configuration as on a
    /// clumped one. That is `O(N²)` in the number of *copies*, and at melt
    /// scale it dominates everything else the objective does (measured: 6x the
    /// whole evaluation at 1k copies, 68x at 10k).
    ///
    /// Binning the centres into cells at least `d_min` wide reduces it to the
    /// pairs that can actually be within the bound: each cell against itself
    /// and its forward neighbours. This is the same partition-and-stencil the
    /// packer's own pair loop uses, on the same primitive, one level up — over
    /// per-copy centres instead of atoms, at `d_min` instead of contact
    /// distance.
    #[inline]
    fn for_each_close_pair<V>(&self, sites: &[[F; 3]], ctx: &GroupCtx<'_>, mut visit: V)
    where
        V: FnMut(usize, usize, [F; 3], F),
    {
        let n = sites.len();
        if n < 2 {
            return;
        }
        let d2 = self.d_min * self.d_min;

        let grid = partition(ctx.cell, self.d_min, n);
        let ncells = grid.n_cells();

        // A partition only pays when it excludes something. The 3x3x3 stencil
        // reaches every cell of a partition with 27 or fewer, so the sweep
        // would examine every pair anyway and the binning would be pure
        // overhead — measured: at 200 copies in a box only three cells wide it
        // cost 1.5x the plain double loop. Below the margin, take the double
        // loop; it is also the whole story for a handful of copies, where
        // N²/2 is nothing.
        if ncells < GRID_MIN_CELLS {
            for a in 0..n {
                for b in a + 1..n {
                    if let Some((delta, dd)) = closer_than(sites, &ctx.mic, a, b, d2) {
                        visit(a, b, delta, dd);
                    }
                }
            }
            return;
        }

        // Counting sort of the sites into cells: `starts` is the CSR row index,
        // `order` the site ids grouped by cell.
        let cell_of: Vec<u32> = sites
            .iter()
            .map(|s| grid.cell_of(ctx.cell, *s) as u32)
            .collect();
        let mut starts = vec![0u32; ncells + 1];
        for &c in &cell_of {
            starts[c as usize + 1] += 1;
        }
        for i in 0..ncells {
            starts[i + 1] += starts[i];
        }
        let mut cursor = starts.clone();
        let mut order = vec![0u32; n];
        for (i, &c) in cell_of.iter().enumerate() {
            order[cursor[c as usize] as usize] = i as u32;
            cursor[c as usize] += 1;
        }

        let mut stencil = [0usize; 27];
        for icell in 0..ncells {
            let (lo, hi) = (starts[icell] as usize, starts[icell + 1] as usize);
            if lo == hi {
                continue;
            }
            // Within the cell: every unordered pair once.
            for a in lo..hi {
                for b in a + 1..hi {
                    let (i, j) = (order[a] as usize, order[b] as usize);
                    if let Some((delta, dd)) = closer_than(sites, &ctx.mic, i, j, d2) {
                        visit(i, j, delta, dd);
                    }
                }
            }
            // Forward neighbours: each adjacent cell pair is visited once, so
            // every cross-cell site pair is too.
            let nn = grid.stencil_forward(icell, &mut stencil);
            for &nc in &stencil[..nn] {
                let (nlo, nhi) = (starts[nc] as usize, starts[nc + 1] as usize);
                for a in lo..hi {
                    for b in nlo..nhi {
                        let (i, j) = (order[a] as usize, order[b] as usize);
                        if let Some((delta, dd)) = closer_than(sites, &ctx.mic, i, j, d2) {
                            visit(i, j, delta, dd);
                        }
                    }
                }
            }
        }
    }
}

impl Restraint for SelfSeparation {
    fn f(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>) -> F {
        let sites = com::centroids(coords, ctx.natoms_per_copy);
        let mut sum = 0.0;
        self.for_each_close_pair(&sites, &ctx, |_, _, _, dd| {
            let shortfall = self.d_min - dd.sqrt();
            sum += shortfall * shortfall;
        });
        self.strength * sum
    }

    fn fg(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>, grads: &mut [[F; 3]]) -> F {
        let sites = com::centroids(coords, ctx.natoms_per_copy);
        let mut dsites = vec![[0.0 as F; 3]; sites.len()];
        let mut sum = 0.0;

        self.for_each_close_pair(&sites, &ctx, |c, other, delta, dd| {
            let dist = dd.sqrt();
            let shortfall = self.d_min - dist;
            sum += shortfall * shortfall;
            if dist <= D_GUARD {
                // Coincident centres: no direction to separate along.
                return;
            }
            // ∂/∂R_c (d_min − D)² = −2 (d_min − D) Δ/D, and −that for R_c'.
            let w = -2.0 * self.strength * shortfall / dist;
            for k in 0..3 {
                let contrib = w * delta[k];
                dsites[c][k] += contrib;
                dsites[other][k] -= contrib;
            }
        });

        com::scatter(&dsites, ctx.natoms_per_copy, grads);
        self.strength * sum
    }

    /// A separation bound is met or not met, and its penalty is exactly zero
    /// once met — so it belongs in the convergence verdict.
    fn is_bound(&self) -> bool {
        true
    }

    fn name(&self) -> &'static str {
        "SelfSeparation"
    }
}

/// Minimum-image displacement and squared distance for one pair, when they are
/// closer than `d_min` (`D² < d2`).
///
/// A free function rather than a closure over the visitor: nesting the two made
/// the visitor call opaque to the inliner and cost ~37 instructions on every
/// pair examined — enough to make the direct sweep 2.4x slower than the loop it
/// replaced (measured at 200 copies).
#[inline(always)]
fn closer_than(sites: &[[F; 3]], mic: &Mic, a: usize, b: usize, d2: F) -> Option<([F; 3], F)> {
    let delta = mic.apply([
        sites[a][0] - sites[b][0],
        sites[a][1] - sites[b][1],
        sites[a][2] - sites[b][2],
    ]);
    let dd = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
    (dd < d2).then_some((delta, dd))
}

/// Fewer cells than this and the 3x3x3 stencil sees the whole partition, so
/// binning cannot exclude a single pair. 27 is the break-even count (the
/// stencil is 27 cells of the `ncells` available); the margin above it pays for
/// the binning pass itself.
const GRID_MIN_CELLS: usize = 64;

/// Below this separation the direction between two centres is treated as
/// undefined, so no force is applied — the same convention the radial
/// distribution geometries use at their centre.
const D_GUARD: F = 1e-9;

/// A partition of `cell` whose cells are **at least `d_min` wide**, so the
/// 3x3x3 stencil around a cell covers everything within the bound.
///
/// Capped at one cell per site (at least 8): a large box with a small `d_min`
/// and few copies would otherwise allocate a cell table orders of magnitude
/// bigger than the point set it indexes, and rebuild it on every evaluation.
/// Coarser cells only ever add candidate pairs, never drop one, so the cap
/// costs selectivity and never correctness.
fn partition(cell: &SimBox, d_min: F, nsites: usize) -> CellGrid {
    CellGrid::for_cutoff_capped(cell, d_min, nsites.max(8))
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{assert_fd_grad_in, rng_uniform};
    use super::*;
    use molrs::spatial::simbox::SimBox;

    fn cube(side: F, periodic: bool) -> SimBox {
        SimBox::cube(side, molrs::types::F3::zeros(3), [periodic; 3]).expect("test box")
    }

    /// A context whose minimum image and partition come from the *same* box —
    /// the pairing the objective always supplies. Handing a periodic `Mic` a
    /// non-periodic partition would make the grid miss exactly the pairs the
    /// `Mic` was added to catch.
    fn ctx<'a>(cell: &'a SimBox, natoms_per_copy: usize) -> GroupCtx<'a> {
        GroupCtx {
            scale: 1.0,
            scale2: 1.0,
            natoms_per_copy,
            cell,
            mic: cell.mic().simplified(),
        }
    }

    /// The `O(N²)` definition of the penalty, kept as the reference the
    /// partitioned sweep is checked against.
    fn brute_force(r: &SelfSeparation, coords: &[[F; 3]], c: GroupCtx<'_>) -> F {
        let sites = com::centroids(coords, c.natoms_per_copy);
        let d2 = r.d_min * r.d_min;
        let mut sum = 0.0;
        for (i, a) in sites.iter().enumerate() {
            for b in &sites[i + 1..] {
                let delta = c.mic.apply([a[0] - b[0], a[1] - b[1], a[2] - b[2]]);
                let dd = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
                if dd < d2 {
                    let shortfall = r.d_min - dd.sqrt();
                    sum += shortfall * shortfall;
                }
            }
        }
        r.strength * sum
    }

    fn random_sites(n: usize, side: F, seed: u64) -> Vec<[F; 3]> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                [
                    rng_uniform(&mut s, 0.0, side),
                    rng_uniform(&mut s, 0.0, side),
                    rng_uniform(&mut s, 0.0, side),
                ]
            })
            .collect()
    }

    // ── the partitioned sweep must be the double loop ──────────────────────

    #[test]
    fn partitioned_sweep_matches_brute_force_free() {
        // Several densities: d_min well below the spacing (few pairs), around
        // it, and above the box (every pair). The last one is what catches a
        // stencil that silently drops far cells.
        let side = 40.0;
        for &n in &[2usize, 17, 200] {
            for &d_min in &[3.0 as F, 12.0, 80.0] {
                let r = SelfSeparation::new(d_min, 1.0);
                let cell = cube(side, false);
                let sites = random_sites(n, side, 0x1234 + n as u64);
                let got = r.f(&sites, ctx(&cell, 1));
                let want = brute_force(&r, &sites, ctx(&cell, 1));
                assert!(
                    (got - want).abs() < 1e-9 * want.max(1.0),
                    "n={n} d_min={d_min}: partitioned {got} != brute force {want}"
                );
            }
        }
    }

    #[test]
    fn partitioned_sweep_matches_brute_force_periodic() {
        let side = 40.0;
        for &n in &[2usize, 17, 200] {
            for &d_min in &[3.0 as F, 12.0, 19.0] {
                let r = SelfSeparation::new(d_min, 1.0);
                let cell = cube(side, true);
                let sites = random_sites(n, side, 0xbeef + n as u64);
                let got = r.f(&sites, ctx(&cell, 1));
                let want = brute_force(&r, &sites, ctx(&cell, 1));
                assert!(
                    (got - want).abs() < 1e-9 * want.max(1.0),
                    "n={n} d_min={d_min}: partitioned {got} != brute force {want}"
                );
            }
        }
    }

    #[test]
    fn partitioned_sweep_matches_brute_force_for_multi_atom_copies() {
        let side = 30.0;
        let r = SelfSeparation::new(9.0, 1.0);
        let cell = cube(side, true);
        // 40 copies of a 3-atom molecule.
        let coords = random_sites(120, side, 0xc0ffee);
        let got = r.f(&coords, ctx(&cell, 3));
        let want = brute_force(&r, &coords, ctx(&cell, 3));
        assert!((got - want).abs() < 1e-9 * want.max(1.0), "{got} != {want}");
    }

    #[test]
    fn sites_outside_the_box_are_still_paired() {
        // Centres drift outside the cell during optimisation; a non-periodic
        // partition clamps them into the edge cells, which must not lose a
        // pair that is genuinely within the bound.
        let cell = cube(20.0, false);
        let r = SelfSeparation::new(6.0, 1.0);
        let sites = [
            [-35.0, 5.0, 5.0],
            [-32.0, 5.0, 5.0],
            [55.0, 5.0, 5.0],
            [57.0, 5.0, 5.0],
        ];
        let got = r.f(&sites, ctx(&cell, 1));
        let want = brute_force(&r, &sites, ctx(&cell, 1));
        assert!(want > 0.0, "the fixture must actually violate the bound");
        assert!((got - want).abs() < 1e-9 * want, "{got} != {want}");
    }

    #[test]
    fn a_huge_box_with_few_copies_stays_cheap() {
        // `for_cutoff` alone would ask for (5000/2)^3 = 1.5e10 cells here. The
        // cap must keep the table proportional to the point set; the answer is
        // unchanged either way.
        let cell = cube(5_000.0, false);
        let r = SelfSeparation::new(2.0, 1.0);
        assert!(partition(&cell, 2.0, 10).n_cells() <= 10);
        let sites = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [900.0, 0.0, 0.0]];
        let got = r.f(&sites, ctx(&cell, 1));
        let want = brute_force(&r, &sites, ctx(&cell, 1));
        assert!((got - want).abs() < 1e-12, "{got} != {want}");
    }

    // ── behaviour ─────────────────────────────────────────────────────────

    #[test]
    fn silent_when_every_copy_is_far_enough() {
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(5.0, 1.0);
        let coords = [[0.0, 0.0, 0.0], [6.0, 0.0, 0.0], [12.0, 0.0, 0.0]];
        assert_eq!(r.f(&coords, ctx(&cell, 1)), 0.0);
    }

    #[test]
    fn penalises_a_clump_and_grows_as_it_tightens() {
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(5.0, 1.0);
        let loose = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let tight = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let loose_e = r.f(&loose, ctx(&cell, 1));
        assert!(loose_e > 0.0);
        assert!(r.f(&tight, ctx(&cell, 1)) > loose_e);
    }

    #[test]
    fn strength_scales_the_value_linearly() {
        let cell = cube(100.0, false);
        let coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let one = SelfSeparation::new(5.0, 1.0).f(&coords, ctx(&cell, 1));
        let three = SelfSeparation::new(5.0, 3.0).f(&coords, ctx(&cell, 1));
        assert!((three - 3.0 * one).abs() < 1e-12);
    }

    #[test]
    fn measures_molecules_not_atoms() {
        // Two 2-atom copies whose centroids sit 10 apart: silent at d_min = 5,
        // even though individual atoms of the two copies are only 8 apart.
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(5.0, 1.0);
        let coords = [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [9.0, 0.0, 0.0],
            [11.0, 0.0, 0.0],
        ];
        assert_eq!(r.f(&coords, ctx(&cell, 2)), 0.0);
    }

    #[test]
    fn sees_across_a_periodic_boundary() {
        // Centres at x = 1 and x = 19 in a 20 Å box are 2 apart, not 18.
        let r = SelfSeparation::new(5.0, 1.0);
        let coords = [[1.0, 0.0, 0.0], [19.0, 0.0, 0.0]];
        let free = cube(20.0, false);
        let pbc = cube(20.0, true);
        assert_eq!(r.f(&coords, ctx(&free, 1)), 0.0, "free: 18 apart");
        assert!(
            r.f(&coords, ctx(&pbc, 1)) > 0.0,
            "periodic: 2 apart, must be penalised"
        );
    }

    // ── gradient ──────────────────────────────────────────────────────────

    #[test]
    fn gradient_matches_finite_differences_monatomic() {
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(6.0, 1.3);
        let coords = random_sites(8, 10.0, 0x5eed_1234);
        assert_fd_grad_in(&r, &coords, ctx(&cell, 1));
    }

    #[test]
    fn gradient_matches_finite_differences_multi_atom_copies() {
        // 5 copies x 3 atoms — exercises the 1/m scatter.
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(7.0, 1.0);
        let coords = random_sites(15, 12.0, 0xabcd_0001);
        assert_fd_grad_in(&r, &coords, ctx(&cell, 3));
    }

    #[test]
    fn gradient_matches_finite_differences_across_a_periodic_boundary() {
        // Deliberately straddling x = 0 / x = 20 so the minimum-image branch is
        // the one under test. Kept clear of the half-box distance, where the
        // image choice flips and the penalty is not differentiable.
        let cell = cube(20.0, true);
        let r = SelfSeparation::new(6.0, 1.0);
        let coords = [
            [1.0, 5.0, 5.0],
            [19.0, 5.2, 4.8],
            [18.0, 15.0, 5.0],
            [2.5, 15.3, 5.4],
        ];
        assert_fd_grad_in(&r, &coords, ctx(&cell, 1));
    }

    #[test]
    fn gradient_is_zero_when_nothing_is_violated() {
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(2.0, 1.0);
        let coords = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let mut grads = [[0.0 as F; 3]; 2];
        assert_eq!(r.fg(&coords, ctx(&cell, 1), &mut grads), 0.0);
        assert_eq!(grads, [[0.0; 3]; 2]);
    }

    #[test]
    fn gradient_pushes_a_pair_apart() {
        // Descent is -∇E; for the left copy that must point further left.
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(5.0, 1.0);
        let coords = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut grads = [[0.0 as F; 3]; 2];
        r.fg(&coords, ctx(&cell, 1), &mut grads);
        assert!(grads[0][0] > 0.0, "-grad moves copy 0 toward -x");
        assert!(grads[1][0] < 0.0, "-grad moves copy 1 toward +x");
    }

    #[test]
    fn pair_forces_cancel() {
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(5.0, 1.0);
        let coords = [[0.0, 0.0, 0.0], [1.0, 2.0, 0.5]];
        let mut grads = [[0.0 as F; 3]; 2];
        r.fg(&coords, ctx(&cell, 1), &mut grads);
        let [ga, gb] = grads;
        for (a, b) in ga.iter().zip(gb.iter()) {
            assert!((a + b).abs() < 1e-12);
        }
    }

    #[test]
    fn value_and_fused_value_agree() {
        let cell = cube(100.0, false);
        let r = SelfSeparation::new(6.0, 0.7);
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 3.0, 0.0],
            [2.0, 1.0, 1.0],
        ];
        let mut grads = [[0.0 as F; 3]; 4];
        let fused = r.fg(&coords, ctx(&cell, 1), &mut grads);
        assert!((fused - r.f(&coords, ctx(&cell, 1))).abs() < 1e-12);
    }

    // ── contract ──────────────────────────────────────────────────────────

    #[test]
    fn declares_itself_a_bound() {
        // The counterpart assertion — that a distribution target is *not* a
        // bound — lives next to `GaussianPlane`; the asymmetry is the point.
        assert!(SelfSeparation::new(5.0, 1.0).is_bound());
    }

    #[test]
    #[should_panic(expected = "d_min must be positive")]
    fn rejects_non_positive_distance() {
        SelfSeparation::new(0.0, 1.0);
    }

    #[test]
    #[should_panic(expected = "strength must be positive")]
    fn rejects_non_positive_strength() {
        SelfSeparation::new(1.0, 0.0);
    }
}
