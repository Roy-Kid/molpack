//! Tree embeddings on the diamond lattice: the guarded walk and the
//! last-resort zigzag. Both score a side-branch continuation the same way.

use molrs::types::F;
use rand::rngs::SmallRng;

use super::super::config::SEED_TRIES;
use super::{
    DiamondLattice, RisWeights, SawField, T_STEPS, Walk, abs_step, g_state, slot_angle, t_index,
};
use crate::grow::internal::wrap_pi;
use crate::random::uniform01;

/// Best unused tetrahedral continuation of `parent`, scored by how close its
/// RIS slot sits to `want`. `reject` drops a candidate absolute site; the
/// caller decides what "occupied" means and wraps the site itself.
fn best_continuation(
    lat: &DiamondLattice,
    parent: [i64; 3],
    used: &[bool],
    ggp_step: Option<[i64; 3]>,
    gp_step: Option<[i64; 3]>,
    want: F,
    mut reject: impl FnMut([i64; 3]) -> bool,
) -> Option<(usize, [i64; 3])> {
    let mut best: Option<(F, usize, [i64; 3])> = None;
    for (ti, &t) in T_STEPS.iter().enumerate() {
        if used[ti] {
            continue;
        }
        let qabs = abs_step(lat, parent, t);
        if reject(qabs) {
            continue;
        }
        let step = [
            qabs[0] - parent[0],
            qabs[1] - parent[1],
            qabs[2] - parent[2],
        ];
        let st = match (ggp_step, gp_step) {
            (Some(u), Some(v)) => g_state(u, v, step),
            _ => 0,
        };
        let err = wrap_pi(slot_angle(st) - want).abs();
        if best.map(|(e, _, _)| err < e).unwrap_or(true) {
            best = Some((err, ti, qabs));
        }
    }
    best.map(|(_, ti, q)| (ti, q))
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
        // Placement stack: (atom, wrapped site) in the order atoms were
        // placed. Every atom is placed after its parent, so any prefix is a
        // parent-closed partial walk and a recoil pops a suffix.
        let mut stack: Vec<(usize, [i64; 3])> = Vec::new();

        let clear = |field: &mut SawField, stack: &[(usize, [i64; 3])]| {
            for (_, p) in stack {
                field.occ.remove(p);
            }
        };

        let mut seeded = None;
        for _ in 0..SEED_TRIES {
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
        stack.push((0, p0));
        stack.push((1, p1w));

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
            stack.push((2, qw));
            placed2 = true;
            break;
        }
        if !placed2 {
            clear(field, &stack);
            continue;
        }

        let bond_into = |abs: &[Option<[i64; 3]>], j: usize| -> Option<[i64; 3]> {
            let p = parent[j]?;
            let a = abs[j]?;
            let b = abs[p]?;
            Some([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
        };

        // Recoil: a dead end pops the most recent placements off the stack
        // and regrows from there, instead of throwing the whole chain away.
        // The depth doubles every two consecutive dead ends that get no
        // further than the deepest one so far; `max_backtrack` recoils buy
        // one reseed.
        let mut backtracks = 0usize;
        let mut streak = 0usize;
        let mut wall = 0usize;
        let mut failed = false;
        let mut p = 0usize;
        'place: while p < n {
            let unplaced: Vec<usize> = children[p]
                .iter()
                .copied()
                .filter(|&c| abs[c].is_none())
                .collect();
            if unplaced.is_empty() {
                p += 1;
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
            let mut dead_end = false;

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
                    dead_end = true;
                }
                if !dead_end {
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
            }

            let mut others: Vec<usize> = unplaced
                .iter()
                .copied()
                .filter(|&c| Some(c) != hooked)
                .collect();
            others.sort_unstable();
            if dead_end {
                others.clear();
            }
            for c in others {
                let want = follows[c].map(|(_, off)| off).unwrap_or(0.0);
                let Some((ti, qabs)) =
                    best_continuation(lat, pabs, &used_now, ggp_step, gp_step, want, |q| {
                        !field.free_for(lat, lat.wrap(q), Some(pw), guard)
                    })
                else {
                    dead_end = true;
                    break;
                };
                used_now[ti] = true;
                assign.push((c, qabs, lat.wrap(qabs)));
            }

            if dead_end {
                // Nothing of `assign` is committed yet, so only the stack
                // beyond the seed triple can be popped.
                if backtracks >= max_backtrack {
                    failed = true;
                    break 'place;
                }
                backtracks += 1;
                if stack.len() > wall {
                    wall = stack.len();
                    streak = 0;
                } else {
                    streak += 1;
                }
                let depth = (1usize << (streak / 2).min(16)).min(stack.len() - 3);
                if depth == 0 {
                    failed = true;
                    break 'place;
                }
                let mut resume = n;
                for _ in 0..depth {
                    let (a, w) = stack.pop().expect("depth never reaches the seed");
                    field.occ.remove(&w);
                    abs[a] = None;
                    resume = resume.min(parent[a].expect("only the root has no parent"));
                }
                p = resume;
                continue 'place;
            }
            for (c, qabs, qw) in assign {
                abs[c] = Some(qabs);
                field.occ.insert(qw, chain_id);
                stack.push((c, qw));
            }
            p += 1;
        }

        if !failed && abs.iter().all(|s| s.is_some()) {
            return Some(Walk {
                sites: abs.into_iter().map(|s| s.unwrap()).collect(),
            });
        }
        clear(field, &stack);
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
    for _ in 0..SEED_TRIES {
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
            let want = follows[c].map(|(_, off)| off).unwrap_or(0.0);
            let (ti, q) = best_continuation(lat, abs[p], &used, ggp_step, gp_step, want, |q| {
                field.blocked.contains(&lat.wrap(q))
            })?;
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
