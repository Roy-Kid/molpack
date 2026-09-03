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

use std::collections::HashMap;

use molrs::types::F;
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
    /// to the nearest commensurate multiple of 4 units (the on-lattice bond
    /// length enters only the topology, never the decorated geometry).
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

    fn wrap(&self, p: [i64; 3]) -> [i64; 3] {
        [
            p[0].rem_euclid(self.m[0]),
            p[1].rem_euclid(self.m[1]),
            p[2].rem_euclid(self.m[2]),
        ]
    }

    fn is_a(p: [i64; 3]) -> bool {
        (p[0] + p[1] + p[2]).rem_euclid(4) == 0
    }

    fn neighbors(&self, p: [i64; 3]) -> [[i64; 3]; 4] {
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
}

/// Occupancy over wrapped sites.
pub(crate) struct SawField {
    occ: HashMap<[i64; 3], u32>,
}

impl SawField {
    pub(crate) fn new() -> Self {
        Self {
            occ: HashMap::new(),
        }
    }

    fn free_for(
        &self,
        lat: &DiamondLattice,
        p: [i64; 3],
        bonded: Option<[i64; 3]>,
        guard: bool,
    ) -> bool {
        if self.occ.contains_key(&p) {
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

/// Grow one self-avoiding walk of `n` sites. Returns `None` when the walk
/// could not complete within the recoil/reseed budget — the caller decides
/// the escape (guard off, then forced zigzag) and its honest accounting.
#[allow(clippy::too_many_arguments)]
pub(crate) fn grow_walk(
    lat: &DiamondLattice,
    field: &mut SawField,
    chain_id: u32,
    n: usize,
    weights: &RisWeights,
    guard: bool,
    max_backtrack: usize,
    max_reseed: usize,
    rng: &mut SmallRng,
) -> Option<Walk> {
    debug_assert!(n >= 2);
    for _ in 0..max_reseed {
        // Seed pair: a free site plus a free neighbour.
        let mut seeded = None;
        for _ in 0..2000 {
            let p0 = lat.random_site(rng);
            if !field.free_for(lat, p0, None, guard) {
                continue;
            }
            let nbs = lat.neighbors(p0);
            let free: Vec<[i64; 3]> = nbs
                .into_iter()
                .filter(|&q| field.free_for(lat, q, Some(p0), guard))
                .collect();
            if !free.is_empty() {
                let pick = (uniform01(rng) * free.len() as F) as usize;
                seeded = Some((p0, free[pick.min(free.len() - 1)]));
                break;
            }
        }
        let (p0, p1) = seeded?;

        // wrapped positions for occupancy; unwrapped for geometry.
        let mut wrapped = vec![p0, p1];
        let mut abs = vec![p0, {
            let s: i64 = if DiamondLattice::is_a(p0) { 1 } else { -1 };
            // Reconstruct the actual step (p1 is wrapped): find the T that
            // maps p0 to p1 under the wrap.
            let mut step = [0i64; 3];
            for t in T_STEPS {
                let q = lat.wrap([p0[0] + s * t[0], p0[1] + s * t[1], p0[2] + s * t[2]]);
                if q == p1 {
                    step = [s * t[0], s * t[1], s * t[2]];
                    break;
                }
            }
            [p0[0] + step[0], p0[1] + step[1], p0[2] + step[2]]
        }];
        let mut steps: Vec<[i64; 3]> = vec![[
            abs[1][0] - abs[0][0],
            abs[1][1] - abs[0][1],
            abs[1][2] - abs[0][2],
        ]];
        let mut states: Vec<i8> = Vec::new();
        field.occ.insert(p0, chain_id);
        field.occ.insert(p1, chain_id);
        let mut tried: Vec<Vec<[i64; 3]>> = vec![Vec::new(), Vec::new()];
        let mut backtracks = 0usize;

        while wrapped.len() < n {
            let cur_w = *wrapped.last().unwrap();
            let cur_a = *abs.last().unwrap();
            let prev_w = wrapped[wrapped.len() - 2];
            let s: i64 = if DiamondLattice::is_a(cur_w) { 1 } else { -1 };

            // Candidates: the 3 non-backtracking continuations.
            let mut cands: Vec<([i64; 3], [i64; 3], i8, F)> = Vec::new(); // (wrapped, step, state, w)
            for t in T_STEPS {
                let step = [s * t[0], s * t[1], s * t[2]];
                let q = lat.wrap([cur_w[0] + step[0], cur_w[1] + step[1], cur_w[2] + step[2]]);
                if q == prev_w || tried.last().unwrap().contains(&q) {
                    continue;
                }
                if !field.free_for(lat, q, Some(cur_w), guard) {
                    continue;
                }
                let (state, w) = if steps.len() >= 2 {
                    let st = g_state(steps[steps.len() - 2], steps[steps.len() - 1], step);
                    // Pentane exclusion: adjacent opposite gauches are the
                    // g+g− collision RIS forbids — never propose them.
                    if let Some(&prev) = states.last()
                        && prev != 0
                        && st != 0
                        && st != prev
                    {
                        continue;
                    }
                    (st, if st == 0 { weights.p_t } else { weights.p_g })
                } else {
                    (0, 1.0) // second bond fixes the angle only
                };
                cands.push((q, step, state, w));
            }

            if cands.is_empty() {
                backtracks += 1;
                if backtracks > max_backtrack || wrapped.len() <= 2 {
                    break; // give up this attempt, reseed
                }
                let dead = wrapped.pop().unwrap();
                abs.pop();
                steps.pop();
                states.pop();
                field.occ.remove(&dead);
                tried.pop();
                tried.last_mut().unwrap().push(dead);
                continue;
            }

            let total: F = cands.iter().map(|c| c.3).sum();
            let mut ticket = uniform01(rng) * total;
            let mut pick = cands.len() - 1;
            for (i, c) in cands.iter().enumerate() {
                ticket -= c.3;
                if ticket <= 0.0 {
                    pick = i;
                    break;
                }
            }
            let (q, step, state, _) = cands[pick];
            wrapped.push(q);
            abs.push([cur_a[0] + step[0], cur_a[1] + step[1], cur_a[2] + step[2]]);
            steps.push(step);
            if steps.len() >= 3 {
                states.push(state);
            }
            field.occ.insert(q, chain_id);
            tried.push(Vec::new());
        }

        if wrapped.len() == n {
            return Some(Walk { sites: abs });
        }
        for p in wrapped {
            field.occ.remove(&p);
        }
    }
    None
}

/// Last-resort escape: an all-trans zigzag from a random seed, ignoring
/// occupancy entirely (still recorded, so later chains see it). Never fails;
/// the caller counts it as a relaxation of the constructive guarantee.
pub(crate) fn forced_zigzag(
    lat: &DiamondLattice,
    field: &mut SawField,
    chain_id: u32,
    n: usize,
    rng: &mut SmallRng,
) -> Walk {
    let p0 = lat.random_site(rng);
    let dir_a = (uniform01(rng) * 4.0) as usize % 4;
    let dir_b = (dir_a + 1 + (uniform01(rng) * 3.0) as usize % 3) % 4;
    let mut abs = vec![p0];
    let mut cur = p0;
    for i in 1..n {
        let t = if i % 2 == 1 {
            T_STEPS[dir_a]
        } else {
            T_STEPS[dir_b]
        };
        let s: i64 = if DiamondLattice::is_a(lat.wrap(cur)) {
            1
        } else {
            -1
        };
        cur = [cur[0] + s * t[0], cur[1] + s * t[1], cur[2] + s * t[2]];
        abs.push(cur);
    }
    for p in &abs {
        field.occ.insert(lat.wrap(*p), chain_id);
    }
    Walk { sites: abs }
}
