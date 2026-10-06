//! Self-avoiding walks on the diamond lattice.
//!
//! Sites live in integer units of a/4: A-sites have all-even coordinates
//! with `x+y+z ≡ 0 (mod 4)` and bond via `+T_k`; B-sites have all-odd
//! coordinates with `≡ 3 (mod 4)` and bond via `−T_k`,
//! `T = {(1,1,1),(1,−1,−1),(−1,1,−1),(−1,−1,1)}`. Successive bonds meet at
//! exactly 109.47°, and the three non-backtracking continuations are the
//! trans/gauche± RIS states — trans ⇔ the new bond repeats the bond two
//! back. Sterics: with the occupancy guard on, a site may be taken only when
//! none of its four lattice neighbours holds a non-bonded atom, which keeps
//! every non-bonded pair at ≥ the 2nd-neighbour distance (`a/√2`).
//!
//! Kept here rather than built on molrs's `builder::SelfAvoidingWalk` on
//! purpose (module-responsibility ruling 10): the walk draws from molpack's
//! own counter-based RNG streams in a fixed order, scores continuations with
//! the RIS weights and the occupancy guard above, and is pinned bit for bit by
//! the lattice-growth goldens. A different walker — even a correct one —
//! consumes the streams differently and moves every grown coordinate.

use std::collections::{HashMap, HashSet};

use molrs::op::types::F;
use rand::rngs::SmallRng;

use crate::random::uniform01;

/// The four A-site bond vectors, in a/4 units.
pub(crate) const T_STEPS: [[i64; 3]; 4] = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]];

/// Diamond lattice over an orthorhombic periodic box: per-axis extents in
/// a/4 units (each a multiple of 4 so the sublattice classes survive the
/// wrap) and the Å-per-unit scale that maps sites back into the box.
pub(crate) struct DiamondLattice {
    pub(crate) m: [i64; 3],
    pub(crate) scale: [F; 3],
    pub(crate) origin: [F; 3],
}

impl DiamondLattice {
    /// Fit the lattice to the box: the lattice constant is taken from the
    /// mean backbone bond length (`a = 4·b/√3`) and then per-axis adjusted
    /// to the nearest commensurate multiple of 4 units.
    ///
    /// Decoration seats backbone atoms on the sites, so **this is where the
    /// molecule's backbone bond length and angles come from**. The per-axis
    /// commensurability moves the step off the template's mean bond by
    /// however much the box rounds — a percent when the cell divides kindly,
    /// a few when it does not — and stretches it anisotropically in a
    /// non-cubic cell, so the angles are exactly tetrahedral only when the
    /// three axes round the same way.
    pub(crate) fn fit(origin: [F; 3], lengths: [F; 3], mean_bond: F) -> Self {
        let unit = mean_bond * 4.0 / (3.0 as F).sqrt() / 4.0; // a/4 in Å
        let mut m = [0i64; 3];
        let mut scale = [0.0 as F; 3];
        for k in 0..3 {
            let raw = (lengths[k] / unit / 4.0).round().max(2.0) as i64;
            m[k] = 4 * raw;
            scale[k] = lengths[k] / m[k] as F;
        }
        Self { m, scale, origin }
    }

    /// Continuum position of an (unwrapped) site.
    pub(crate) fn to_continuum(&self, p: [i64; 3]) -> [F; 3] {
        [
            self.origin[0] + p[0] as F * self.scale[0],
            self.origin[1] + p[1] as F * self.scale[1],
            self.origin[2] + p[2] as F * self.scale[2],
        ]
    }

    pub(crate) fn wrap(&self, p: [i64; 3]) -> [i64; 3] {
        [
            p[0].rem_euclid(self.m[0]),
            p[1].rem_euclid(self.m[1]),
            p[2].rem_euclid(self.m[2]),
        ]
    }

    pub(crate) fn is_a(p: [i64; 3]) -> bool {
        (p[0] + p[1] + p[2]).rem_euclid(4) == 0
    }

    pub(crate) fn neighbors(&self, p: [i64; 3]) -> [[i64; 3]; 4] {
        let s: i64 = if Self::is_a(p) { 1 } else { -1 };
        std::array::from_fn(|k| {
            let t = T_STEPS[k];
            self.wrap([p[0] + s * t[0], p[1] + s * t[1], p[2] + s * t[2]])
        })
    }

    /// Draw a random valid A-site (all-even, sum ≡ 0 mod 4).
    fn random_site(&self, rng: &mut SmallRng) -> [i64; 3] {
        loop {
            let p = [
                2 * (uniform01(rng) * (self.m[0] / 2) as F) as i64,
                2 * (uniform01(rng) * (self.m[1] / 2) as F) as i64,
                2 * (uniform01(rng) * (self.m[2] / 2) as F) as i64,
            ];
            if (p[0] + p[1] + p[2]).rem_euclid(4) == 0 {
                return self.wrap(p);
            }
        }
    }

    /// Every A-site (`x+y+z ≡ 0 (mod 4)`, all even) and B-site
    /// (`≡ 3 (mod 4)`, all odd) in the wrapped cell.
    pub(crate) fn iter_sites(&self) -> impl Iterator<Item = [i64; 3]> + '_ {
        let m = self.m;
        (0..m[0]).flat_map(move |x| {
            (0..m[1]).flat_map(move |y| {
                (0..m[2]).filter_map(move |z| {
                    let even = x % 2 == 0 && y % 2 == 0 && z % 2 == 0;
                    let odd = x % 2 != 0 && y % 2 != 0 && z % 2 != 0;
                    let s = (x + y + z).rem_euclid(4);
                    if (even && s == 0) || (odd && s == 3) {
                        Some([x, y, z])
                    } else {
                        None
                    }
                })
            })
        })
    }

    pub(crate) fn has_allowed_a_site(&self, blocked: &HashSet<[i64; 3]>) -> bool {
        self.iter_sites()
            .any(|p| Self::is_a(p) && !blocked.contains(&p))
    }
}

/// Occupancy over wrapped sites. `blocked` is the complement of
/// Region ∩ lattice (sites outside the attached region); it is not
/// occupancy — a blocked neighbour does not trip the occupancy guard.
pub(crate) struct SawField {
    occ: HashMap<[i64; 3], u32>,
    blocked: HashSet<[i64; 3]>,
}

impl SawField {
    pub(crate) fn new() -> Self {
        Self {
            occ: HashMap::new(),
            blocked: HashSet::new(),
        }
    }

    pub(crate) fn set_blocked(&mut self, blocked: HashSet<[i64; 3]>) {
        self.blocked = blocked;
    }

    fn free_for(
        &self,
        lat: &DiamondLattice,
        p: [i64; 3],
        bonded: Option<[i64; 3]>,
        guard: bool,
    ) -> bool {
        if self.blocked.contains(&p) || self.occ.contains_key(&p) {
            return false;
        }
        if guard {
            for q in lat.neighbors(p) {
                if Some(q) != bonded && self.occ.contains_key(&q) {
                    return false;
                }
            }
        }
        true
    }
}

/// RIS weights for the on-lattice walk.
pub(crate) struct RisWeights {
    pub(crate) p_t: F,
    pub(crate) p_g: F,
}

/// Gauche handedness of the bond triple `(u, v, w)`: 0 for trans
/// (`w == u`), else the sign of the triple product `u · (v × w)`.
fn g_state(u: [i64; 3], v: [i64; 3], w: [i64; 3]) -> i8 {
    if w == u {
        return 0;
    }
    let cx = [
        v[1] * w[2] - v[2] * w[1],
        v[2] * w[0] - v[0] * w[2],
        v[0] * w[1] - v[1] * w[0],
    ];
    let det = u[0] * cx[0] + u[1] * cx[1] + u[2] * cx[2];
    if det > 0 { 1 } else { -1 }
}

/// One finished walk: unwrapped sites plus the wrapped sites it occupies.
pub(crate) struct Walk {
    pub(crate) sites: Vec<[i64; 3]>,
}

fn abs_step(lat: &DiamondLattice, from: [i64; 3], t: [i64; 3]) -> [i64; 3] {
    let s: i64 = if DiamondLattice::is_a(lat.wrap(from)) {
        1
    } else {
        -1
    };
    [from[0] + s * t[0], from[1] + s * t[1], from[2] + s * t[2]]
}

fn t_index(lat: &DiamondLattice, from: [i64; 3], to: [i64; 3]) -> Option<usize> {
    T_STEPS.iter().position(|&t| abs_step(lat, from, t) == to)
}

fn slot_angle(st: i8) -> F {
    match st {
        0 => PI,
        1 => PI / 3.0,
        _ => -PI / 3.0,
    }
}

const PI: F = std::f64::consts::PI as F;

mod walk;
pub(crate) use walk::{forced_zigzag, grow_walk};

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::type_complexity)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    fn linear_tree(n: usize) -> (Vec<Option<usize>>, Vec<Vec<usize>>, Vec<Option<(usize, F)>>) {
        let mut parent = vec![None; n];
        let mut children = vec![Vec::new(); n];
        let mut follows = vec![None; n];
        for j in 1..n {
            parent[j] = Some(j - 1);
            children[j - 1].push(j);
        }
        for j in 3..n {
            follows[j] = Some((j - 3, 0.0));
        }
        (parent, children, follows)
    }

    #[test]
    fn t_steps_are_tetrahedral() {
        for t in T_STEPS {
            assert_eq!(t[0] * t[0] + t[1] * t[1] + t[2] * t[2], 3);
        }
        for (i, ti) in T_STEPS.iter().enumerate() {
            for (j, tj) in T_STEPS.iter().enumerate() {
                if i == j {
                    continue;
                }
                assert_eq!(ti[0] * tj[0] + ti[1] * tj[1] + ti[2] * tj[2], -1);
            }
        }
    }

    #[test]
    fn grow_walk_k1_sites_are_neighbours() {
        let lat = DiamondLattice::fit([0.0; 3], [40.0; 3], 1.53);
        let mut field = SawField::new();
        let mut rng = SmallRng::seed_from_u64(7);
        let (parent, children, follows) = linear_tree(12);
        let weights = RisWeights {
            p_t: 1.0 / 3.0,
            p_g: 1.0 / 3.0,
        };
        let walk = grow_walk(
            &lat, &mut field, 0, &parent, &children, &follows, &weights, true, 20_000, 200,
            &mut rng,
        )
        .expect("linear walk");
        assert_eq!(walk.sites.len(), 12);
        for j in 1..12 {
            let p = parent[j].unwrap();
            assert!(
                t_index(&lat, walk.sites[p], walk.sites[j]).is_some(),
                "atom {j} is not a diamond neighbour of its parent"
            );
        }
    }

    /// A dead end recoils instead of discarding the chain: one attempt
    /// (`max_reseed = 1`) still completes a long guarded walk, which simple
    /// sampling — abandon at the first dead end — almost never does.
    #[test]
    fn grow_walk_recoils_out_of_dead_ends() {
        let lat = DiamondLattice::fit([0.0; 3], [20.0; 3], 1.53);
        let (parent, children, follows) = linear_tree(300);
        let weights = RisWeights { p_t: 0.6, p_g: 0.2 };
        for seed in 0..8 {
            let mut field = SawField::new();
            let mut rng = SmallRng::seed_from_u64(seed);
            let walk = grow_walk(
                &lat, &mut field, 0, &parent, &children, &follows, &weights, true, 20_000, 1,
                &mut rng,
            );
            let walk = walk.unwrap_or_else(|| panic!("seed {seed}: the single attempt gave up"));
            for j in 1..300 {
                let p = parent[j].unwrap();
                assert!(t_index(&lat, walk.sites[p], walk.sites[j]).is_some());
            }
            let mut wrapped: Vec<[i64; 3]> = walk.sites.iter().map(|&s| lat.wrap(s)).collect();
            wrapped.sort_unstable();
            wrapped.dedup();
            assert_eq!(wrapped.len(), 300, "seed {seed}: the walk revisits a site");
            assert_eq!(field.occ.len(), 300, "seed {seed}: recoil leaked occupancy");
        }
    }

    #[test]
    fn forced_zigzag_embeds_a_star() {
        let lat = DiamondLattice::fit([0.0; 3], [40.0; 3], 1.53);
        let mut field = SawField::new();
        let mut rng = SmallRng::seed_from_u64(3);
        // Leaf-rooted 4-arm star: 0 leaf, 1 centre, 2..=4 leaves.
        let parent = vec![None, Some(0), Some(1), Some(1), Some(1)];
        let children = vec![vec![1], vec![2, 3, 4], vec![], vec![], vec![]];
        let follows = vec![None; 5];
        let walk = forced_zigzag(&lat, &mut field, 0, &parent, &children, &follows, &mut rng)
            .expect("star");
        assert_eq!(walk.sites.len(), 5);
        for j in 1..5 {
            let p = parent[j].unwrap();
            assert!(t_index(&lat, walk.sites[p], walk.sites[j]).is_some());
        }
        // Centre (1) has 4 neighbours: parent 0 + 3 children.
        let mut dirs = Vec::new();
        for &c in &[0usize, 2, 3, 4] {
            dirs.push(t_index(&lat, walk.sites[1], walk.sites[c]).expect("spoke"));
        }
        dirs.sort_unstable();
        dirs.dedup();
        assert_eq!(dirs.len(), 4);
    }

    #[test]
    fn walk_skips_blocked_half_space() {
        let lat = DiamondLattice::fit([0.0; 3], [20.0; 3], 1.53);
        let mut blocked = HashSet::new();
        for p in lat.iter_sites() {
            if lat.to_continuum(p)[0] > 10.0 {
                blocked.insert(p);
            }
        }
        assert!(!blocked.is_empty());
        assert!(lat.has_allowed_a_site(&blocked));
        let mut field = SawField::new();
        field.set_blocked(blocked);
        let mut rng = SmallRng::seed_from_u64(7);
        let (parent, children, follows) = linear_tree(8);
        let weights = RisWeights {
            p_t: 1.0 / 3.0,
            p_g: 1.0 / 3.0,
        };
        let walk = grow_walk(
            &lat, &mut field, 0, &parent, &children, &follows, &weights, true, 20_000, 200,
            &mut rng,
        )
        .expect("walk inside x ≤ 10");
        for s in walk.sites {
            let x = lat.to_continuum(lat.wrap(s))[0];
            assert!(x <= 10.0 + 1e-9, "site continuum x={x} crossed the mask");
        }
    }
}
