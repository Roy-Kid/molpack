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

use std::collections::{HashMap, HashSet};

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

fn wrap_pi(x: F) -> F {
    let mut v = x % (2.0 * PI);
    if v > PI {
        v -= 2.0 * PI;
    } else if v <= -PI {
        v += 2.0 * PI;
    }
    v
}

/// Grow one self-avoiding tree walk. `parent` / `children` / `follows` are
/// the InternalTree heavy projection (`Backbone`). Returns `None` when the
/// walk could not complete within the recoil/reseed budget.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub(crate) fn grow_walk(
    lat: &DiamondLattice,
    field: &mut SawField,
    chain_id: u32,
    parent: &[Option<usize>],
    children: &[Vec<usize>],
    follows: &[Option<(usize, F)>],
    weights: &RisWeights,
    guard: bool,
    max_backtrack: usize,
    max_reseed: usize,
    rng: &mut SmallRng,
) -> Option<Walk> {
    let n = parent.len();
    debug_assert!(n >= 3);
    debug_assert_eq!(children.len(), n);
    debug_assert_eq!(follows.len(), n);

    for _ in 0..max_reseed {
        let mut abs: Vec<Option<[i64; 3]>> = vec![None; n];
        let mut occupied: Vec<[i64; 3]> = Vec::new();

        let clear = |field: &mut SawField, occupied: &[[i64; 3]]| {
            for p in occupied {
                field.occ.remove(p);
            }
        };

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
        let Some((p0, p1w)) = seeded else {
            continue;
        };
        let p1 = {
            let mut found = None;
            for t in T_STEPS {
                let q = abs_step(lat, p0, t);
                if lat.wrap(q) == p1w {
                    found = Some(q);
                    break;
                }
            }
            found.expect("neighbour is a T_STEPS image")
        };
        abs[0] = Some(p0);
        abs[1] = Some(p1);
        field.occ.insert(p0, chain_id);
        field.occ.insert(p1w, chain_id);
        occupied.push(p0);
        occupied.push(p1w);

        // Align[2]: one non-backtracking continuation from atom 1.
        let mut placed2 = false;
        for t in T_STEPS {
            let qabs = abs_step(lat, p1, t);
            let qw = lat.wrap(qabs);
            if qw == p0 {
                continue;
            }
            if !field.free_for(lat, qw, Some(p1w), guard) {
                continue;
            }
            abs[2] = Some(qabs);
            field.occ.insert(qw, chain_id);
            occupied.push(qw);
            placed2 = true;
            break;
        }
        if !placed2 {
            clear(field, &occupied);
            continue;
        }

        let bond_into = |abs: &[Option<[i64; 3]>], j: usize| -> Option<[i64; 3]> {
            let p = parent[j]?;
            let a = abs[j]?;
            let b = abs[p]?;
            Some([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
        };

        let mut backtracks = 0usize;
        let mut failed = false;
        'place: for p in 0..n {
            let unplaced: Vec<usize> = children[p]
                .iter()
                .copied()
                .filter(|&c| abs[c].is_none())
                .collect();
            if unplaced.is_empty() {
                continue;
            }
            let Some(pabs) = abs[p] else {
                failed = true;
                break;
            };
            let pw = lat.wrap(pabs);
            let incoming = parent[p].and_then(|gp| abs[gp]);

            let mut used = [false; 4];
            if let Some(gp) = incoming
                && let Some(idx) = t_index(lat, pabs, gp)
            {
                used[idx] = true;
            }
            for &c in &children[p] {
                if let Some(ca) = abs[c]
                    && let Some(idx) = t_index(lat, pabs, ca)
                {
                    used[idx] = true;
                }
            }

            let gp_step = parent[p].and_then(|_gp| bond_into(&abs, p));
            let ggp_step = parent[p].and_then(|gp| bond_into(&abs, gp));

            let hooked = unplaced
                .iter()
                .copied()
                .filter(|&c| follows[c].is_some())
                .min();

            let mut assign: Vec<(usize, [i64; 3], [i64; 3])> = Vec::new(); // (child, abs, wrapped)
            let mut used_now = used;

            if let Some(h) = hooked {
                struct Cand {
                    ti: usize,
                    qabs: [i64; 3],
                    qw: [i64; 3],
                    wgt: F,
                }
                let mut cands: Vec<Cand> = Vec::new();
                for (ti, &t) in T_STEPS.iter().enumerate() {
                    if used_now[ti] {
                        continue;
                    }
                    let qabs = abs_step(lat, pabs, t);
                    let qw = lat.wrap(qabs);
                    if !field.free_for(lat, qw, Some(pw), guard) {
                        continue;
                    }
                    let step = [qabs[0] - pabs[0], qabs[1] - pabs[1], qabs[2] - pabs[2]];
                    let (_state, wgt) = if let (Some(u), Some(v)) = (ggp_step, gp_step) {
                        let st = g_state(u, v, step);
                        if let Some(prev) = parent[p]
                            .and_then(|gp| parent[gp])
                            .and_then(|ggp| bond_into(&abs, ggp))
                            .map(|uu| g_state(uu, u, v))
                            && prev != 0
                            && st != 0
                            && st != prev
                        {
                            continue;
                        }
                        (st, if st == 0 { weights.p_t } else { weights.p_g })
                    } else {
                        (0, 1.0)
                    };
                    cands.push(Cand { ti, qabs, qw, wgt });
                }
                if cands.is_empty() {
                    backtracks += 1;
                    failed = true;
                    break 'place;
                }
                let total: F = cands.iter().map(|c| c.wgt).sum();
                let mut ticket = uniform01(rng) * total;
                let mut pick = cands.len() - 1;
                for (i, c) in cands.iter().enumerate() {
                    ticket -= c.wgt;
                    if ticket <= 0.0 {
                        pick = i;
                        break;
                    }
                }
                let Cand { ti, qabs, qw, .. } = cands[pick];
                used_now[ti] = true;
                assign.push((h, qabs, qw));
            }

            let mut others: Vec<usize> = unplaced
                .iter()
                .copied()
                .filter(|&c| Some(c) != hooked)
                .collect();
            others.sort_unstable();
            for c in others {
                let mut best: Option<(F, usize, [i64; 3], [i64; 3])> = None;
                for (ti, &t) in T_STEPS.iter().enumerate() {
                    if used_now[ti] {
                        continue;
                    }
                    let qabs = abs_step(lat, pabs, t);
                    let qw = lat.wrap(qabs);
                    if !field.free_for(lat, qw, Some(pw), guard) {
                        continue;
                    }
                    let step = [qabs[0] - pabs[0], qabs[1] - pabs[1], qabs[2] - pabs[2]];
                    let st = if let (Some(u), Some(v)) = (ggp_step, gp_step) {
                        g_state(u, v, step)
                    } else {
                        0
                    };
                    let want = follows[c].map(|(_, off)| off).unwrap_or(0.0);
                    let err = wrap_pi(slot_angle(st) - want).abs();
                    if best.map(|(e, _, _, _)| err < e).unwrap_or(true) {
                        best = Some((err, ti, qabs, qw));
                    }
                }
                let Some((_, ti, qabs, qw)) = best else {
                    failed = true;
                    break 'place;
                };
                used_now[ti] = true;
                assign.push((c, qabs, qw));
            }

            if failed {
                break;
            }
            for (c, qabs, qw) in assign {
                abs[c] = Some(qabs);
                field.occ.insert(qw, chain_id);
                occupied.push(qw);
            }
        }

        if !failed && abs.iter().all(|s| s.is_some()) {
            return Some(Walk {
                sites: abs.into_iter().map(|s| s.unwrap()).collect(),
            });
        }
        clear(field, &occupied);
        if backtracks > max_backtrack {
            continue;
        }
    }
    None
}

/// Last-resort escape: a tetrahedral embedding of the tree from a random
/// seed, ignoring occupancy (still recorded). Honours `field.blocked`
/// (Region ∩ lattice). Returns `None` when no allowed tetrahedral embedding
/// exists inside the region.
#[allow(clippy::needless_range_loop)]
pub(crate) fn forced_zigzag(
    lat: &DiamondLattice,
    field: &mut SawField,
    chain_id: u32,
    parent: &[Option<usize>],
    children: &[Vec<usize>],
    follows: &[Option<(usize, F)>],
    rng: &mut SmallRng,
) -> Option<Walk> {
    let n = parent.len();
    let mut abs = vec![[0i64; 3]; n];
    let mut seed = None;
    for _ in 0..2000 {
        let p0 = lat.random_site(rng);
        if field.blocked.contains(&p0) {
            continue;
        }
        seed = Some(p0);
        break;
    }
    if seed.is_none() {
        seed = lat
            .iter_sites()
            .find(|p| DiamondLattice::is_a(*p) && !field.blocked.contains(p));
    }
    abs[0] = seed?;
    let mut t0 = None;
    for k in 0..4 {
        let ti = (k + (uniform01(rng) * 4.0) as usize) % 4;
        let q = abs_step(lat, abs[0], T_STEPS[ti]);
        if !field.blocked.contains(&lat.wrap(q)) {
            t0 = Some(ti);
            break;
        }
    }
    abs[1] = abs_step(lat, abs[0], T_STEPS[t0?]);
    // Atom 2: prefer trans of the 0→1 bond (the continuation that repeats
    // the previous step after the A/B sign flip is "the other" of a pair).
    let mut placed = vec![false; n];
    placed[0] = true;
    placed[1] = true;

    let mut used1 = [false; 4];
    if let Some(idx) = t_index(lat, abs[1], abs[0]) {
        used1[idx] = true;
    }
    let mut got2 = false;
    for (ti, &t) in T_STEPS.iter().enumerate() {
        if used1[ti] {
            continue;
        }
        let q = abs_step(lat, abs[1], t);
        if field.blocked.contains(&lat.wrap(q)) {
            continue;
        }
        abs[2] = q;
        placed[2] = true;
        got2 = true;
        break;
    }
    if !got2 {
        return None;
    }

    for p in 0..n {
        let unplaced: Vec<usize> = children[p]
            .iter()
            .copied()
            .filter(|&c| !placed[c])
            .collect();
        if unplaced.is_empty() {
            continue;
        }
        let mut used = [false; 4];
        if let Some(gp) = parent[p]
            && let Some(idx) = t_index(lat, abs[p], abs[gp])
        {
            used[idx] = true;
        }
        for &c in &children[p] {
            if placed[c]
                && let Some(idx) = t_index(lat, abs[p], abs[c])
            {
                used[idx] = true;
            }
        }
        let hooked = unplaced
            .iter()
            .copied()
            .filter(|&c| follows[c].is_some())
            .min();
        let gp_step = parent[p].map(|gp| {
            [
                abs[p][0] - abs[gp][0],
                abs[p][1] - abs[gp][1],
                abs[p][2] - abs[gp][2],
            ]
        });
        let ggp_step = parent[p].and_then(|gp| {
            parent[gp].map(|ggp| {
                [
                    abs[gp][0] - abs[ggp][0],
                    abs[gp][1] - abs[ggp][1],
                    abs[gp][2] - abs[ggp][2],
                ]
            })
        });
        if let Some(h) = hooked {
            let mut trans = None;
            let mut fallback = None;
            for (ti, &t) in T_STEPS.iter().enumerate() {
                if used[ti] {
                    continue;
                }
                let q = abs_step(lat, abs[p], t);
                if field.blocked.contains(&lat.wrap(q)) {
                    continue;
                }
                let step = [q[0] - abs[p][0], q[1] - abs[p][1], q[2] - abs[p][2]];
                let is_trans = match (ggp_step, gp_step) {
                    (Some(u), Some(v)) => g_state(u, v, step) == 0,
                    _ => false,
                };
                if is_trans {
                    trans = Some((ti, q));
                    break;
                }
                if fallback.is_none() {
                    fallback = Some((ti, q));
                }
            }
            let (ti, q) = trans.or(fallback)?;
            used[ti] = true;
            abs[h] = q;
            placed[h] = true;
        }
        let mut others: Vec<usize> = unplaced
            .iter()
            .copied()
            .filter(|&c| Some(c) != hooked)
            .collect();
        others.sort_unstable();
        for c in others {
            let mut best: Option<(F, usize, [i64; 3])> = None;
            for (ti, &t) in T_STEPS.iter().enumerate() {
                if used[ti] {
                    continue;
                }
                let q = abs_step(lat, abs[p], t);
                if field.blocked.contains(&lat.wrap(q)) {
                    continue;
                }
                let step = [q[0] - abs[p][0], q[1] - abs[p][1], q[2] - abs[p][2]];
                let st = match (ggp_step, gp_step) {
                    (Some(u), Some(v)) => g_state(u, v, step),
                    _ => 0,
                };
                let want = follows[c].map(|(_, off)| off).unwrap_or(0.0);
                let err = wrap_pi(slot_angle(st) - want).abs();
                if best.map(|(e, _, _)| err < e).unwrap_or(true) {
                    best = Some((err, ti, q));
                }
            }
            let (_, ti, q) = best?;
            used[ti] = true;
            abs[c] = q;
            placed[c] = true;
        }
    }

    for p in &abs {
        field.occ.insert(lat.wrap(*p), chain_id);
    }
    Some(Walk { sites: abs })
}

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
