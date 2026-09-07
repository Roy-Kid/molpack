//! PEO melt evaluation: the chain-growth solver against the rigid-body path.
//!
//! Produces the quantified report backing the chain-growth spec's scientific
//! and performance acceptance criteria (ac-006 / ac-010): timing, convergence,
//! `degraded`, chain Rg statistics vs the Flory unperturbed value, chain-order
//! bias, the internal-distance curve ⟨R²(s)⟩/s, density homogeneity E(d), and
//! the minimum inter-molecular distance.
//!
//! ```sh
//! cargo run --release --example pack_peo --features io -- grow 25 200 0.5 42 80 24
//! cargo run --release --example pack_peo --features io -- rigid 200 25 0.3 42
//! cargo run --release --example pack_peo --features io -- grow 200 25 1.0 42 chain.pdb
//! ```
//!
//! One distance knob: `PEO_TOLERANCE` (default 2.0 Å). Every atom's packing
//! radius is `tolerance/2`; contact is `tolerance`. The growth ladder is
//! pinned at full contact (`min_hard_scale = 1`). Optional `PEO_H_RADIUS`
//! still exists but is off unless set.
//!
//! Without a PDB argument the program synthesizes an all-atom PEO chain
//! H-(CH₂-CH₂-O)ₙ-H itself (1.53/1.43/1.10 Å bonds, tetrahedral angles,
//! all-trans start) — dp = 200 gives the spec's 1402 atoms per chain. The
//! template's *shape* is irrelevant to growth (only its chemistry is used),
//! so the synthetic template and an AmberPolymerBuilder one measure the same
//! thing.

use molpack::grow::{GrowConfig, TorsionPrior};
use molpack::{
    CbmcGrow, F, GenCanPack, Handler, LatticeGrow, PackContext, PackEngine, State, StepInfo, Target,
};
use molrs::store::block::Block;
use molrs::store::frame::Frame;
use ndarray::Array1;
use std::sync::{Arc, Mutex};
use std::time::Instant;

/// Live growth ladder: round count and whether the hard-core scale dropped.
struct Ladder {
    start: Instant,
    rounds: usize,
    rad_min: F,
    first_drop: Option<(usize, F)>,
    last_hist: Vec<usize>,
}

struct LadderHandler {
    state: Arc<Mutex<Ladder>>,
    trace: Option<std::io::BufWriter<std::fs::File>>,
}

impl LadderHandler {
    fn new(state: Arc<Mutex<Ladder>>, trace_path: Option<&str>) -> Self {
        let trace = trace_path.map(|p| {
            let f = std::fs::File::create(p).expect("PEO_TRACE");
            let mut w = std::io::BufWriter::new(f);
            use std::io::Write;
            let _ = writeln!(
                w,
                "round,elapsed_s,n_placed,frac,n_done,placed_min,placed_med,placed_max,min_inter,max_overlap,handler_fdist"
            );
            w
        });
        Self { state, trace }
    }
}

impl Handler for LadderHandler {
    fn on_step(&mut self, info: &StepInfo, sys: &PackContext) {
        {
            let mut s = self.state.lock().expect("ladder");
            s.rounds = info.loop_idx;
            if info.radscale < 1.0 && s.first_drop.is_none() {
                s.first_drop = Some((info.loop_idx, info.radscale));
            }
            if info.radscale < s.rad_min {
                s.rad_min = info.radscale;
            }
            if info.loop_idx == 1 || info.loop_idx.is_multiple_of(200) {
                eprintln!(
                    "    round {:>5}  radscale={:.3}  ({:.1}s)",
                    info.loop_idx,
                    info.radscale,
                    s.start.elapsed().as_secs_f64()
                );
            }
        }
        let every = if info.loop_idx < 100 { 1 } else { 20 };
        if self.trace.is_some() && (info.loop_idx == 1 || info.loop_idx.is_multiple_of(every)) {
            let snap = placed_snapshot(sys);
            if let Some(w) = self.trace.as_mut() {
                use std::io::Write;
                let elapsed = self
                    .state
                    .lock()
                    .expect("ladder")
                    .start
                    .elapsed()
                    .as_secs_f64();
                let _ = writeln!(
                    w,
                    "{},{:.3},{},{:.4},{},{},{},{},{:.4},{:.4},{:.4e}",
                    info.loop_idx,
                    elapsed,
                    snap.n_placed,
                    snap.frac,
                    snap.n_done,
                    snap.placed_min,
                    snap.placed_med,
                    snap.placed_max,
                    snap.min_inter,
                    snap.max_overlap,
                    info.fdist
                );
            }
            self.state.lock().expect("ladder").last_hist = snap.per_chain;
        }
    }
}

struct PlacedSnap {
    n_placed: usize,
    frac: F,
    n_done: usize,
    placed_min: usize,
    placed_med: usize,
    placed_max: usize,
    min_inter: F,
    max_overlap: F,
    per_chain: Vec<usize>,
}

fn is_placed(p: [F; 3]) -> bool {
    p[0] * p[0] + p[1] * p[1] + p[2] * p[2] > 1e-16
}

fn placed_snapshot(sys: &PackContext) -> PlacedSnap {
    let na = sys.natoms.first().copied().unwrap_or(0);
    let n_chains = sys.nmols.first().copied().unwrap_or(0);
    let base0 = sys.idfirst.first().copied().unwrap_or(0);
    let mut per_chain = vec![0usize; n_chains];
    let mut placed_idx = Vec::new();
    for (c, count) in per_chain.iter_mut().enumerate() {
        let base = base0 + c * na;
        for a in 0..na {
            if is_placed(sys.xcart[base + a]) {
                *count += 1;
                placed_idx.push(base + a);
            }
        }
    }
    let n_placed = placed_idx.len();
    let n_done = per_chain.iter().filter(|&&n| n == na).count();
    let mut sorted = per_chain.clone();
    sorted.sort_unstable();
    let placed_min = sorted.first().copied().unwrap_or(0);
    let placed_max = sorted.last().copied().unwrap_or(0);
    let placed_med = if sorted.is_empty() {
        0
    } else {
        sorted[sorted.len() / 2]
    };
    let l = sys.simbox.lengths()[0];
    let (min_inter, max_overlap) = min_inter_placed(sys, &placed_idx, l);
    PlacedSnap {
        n_placed,
        frac: n_placed as F / (na * n_chains).max(1) as F,
        n_done,
        placed_min,
        placed_med,
        placed_max,
        min_inter,
        max_overlap,
        per_chain,
    }
}

/// Depth 3 exempts hops 1–3 (1-2 / 1-3 / 1-4). Hop 4 is 1-5 and is scored.
fn print_closest_scored(pos: &[[F; 3]], bonds: &[(usize, usize)], elem: &[String], na: usize) {
    let mut adj = vec![Vec::new(); na];
    for &(i, j) in bonds {
        if i < na && j < na {
            adj[i].push(j);
            adj[j].push(i);
        }
    }
    let hops_and_parent = |src: usize| -> (Vec<i32>, Vec<i32>) {
        let mut d = vec![-1i32; na];
        let mut parent = vec![-1i32; na];
        d[src] = 0;
        let mut q = vec![src];
        let mut qi = 0;
        while qi < q.len() {
            let a = q[qi];
            qi += 1;
            for &b in &adj[a] {
                if d[b] < 0 {
                    d[b] = d[a] + 1;
                    parent[b] = a as i32;
                    q.push(b);
                }
            }
        }
        (d, parent)
    };
    let path_str = |i: usize, j: usize, parent: &[i32]| -> String {
        let mut nodes = vec![j];
        let mut cur = j as i32;
        while cur >= 0 && cur as usize != i {
            cur = parent[cur as usize];
            if cur < 0 {
                break;
            }
            nodes.push(cur as usize);
        }
        nodes.reverse();
        nodes
            .iter()
            .map(|&k| format!("{}{}", elem.get(k).map(String::as_str).unwrap_or("?"), k))
            .collect::<Vec<_>>()
            .join("-")
    };
    let mut pairs: Vec<(F, usize, usize, i32, String)> = Vec::new();
    for i in 0..na {
        let (d, parent) = hops_and_parent(i);
        for j in (i + 1)..na {
            if d[j] >= 4 {
                let dr = vsub(pos[i], pos[j]);
                let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                pairs.push((r, i, j, d[j], path_str(i, j, &parent)));
            }
        }
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    println!("  closest scored intra (1-5+; depth-3 skips only 1-2/1-3/1-4):");
    for (r, i, j, hop, path) in pairs.iter().take(5) {
        println!(
            "    {}{}–{}{}  1-{}  {:.3} Å   {}",
            elem.get(*i).map(String::as_str).unwrap_or("?"),
            i,
            elem.get(*j).map(String::as_str).unwrap_or("?"),
            j,
            hop + 1,
            r,
            path
        );
    }
}

fn min_inter_placed(sys: &PackContext, placed: &[usize], l: F) -> (F, F) {
    if placed.len() < 2 {
        return (F::INFINITY, 0.0);
    }
    let na = sys.natoms.first().copied().unwrap_or(1);
    let cut = 4.0 as F;
    let nc = ((l / cut).floor() as usize).max(1);
    let cell = |p: [F; 3]| -> (usize, usize, usize) {
        let f = |v: F| (((v / l).rem_euclid(1.0) * nc as F) as usize).min(nc - 1);
        (f(p[0]), f(p[1]), f(p[2]))
    };
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); nc * nc * nc];
    let flat = |c: (usize, usize, usize)| (c.0 * nc + c.1) * nc + c.2;
    for &i in placed {
        bins[flat(cell(sys.xcart[i]))].push(i);
    }
    let mut best = F::INFINITY;
    let mut max_ov = 0.0 as F;
    let nci = nc as isize;
    for cx in 0..nc {
        for cy in 0..nc {
            for cz in 0..nc {
                let list = &bins[flat((cx, cy, cz))];
                for dx in -1..=1isize {
                    for dy in -1..=1isize {
                        for dz in -1..=1isize {
                            let nb = (
                                (cx as isize + dx).rem_euclid(nci) as usize,
                                (cy as isize + dy).rem_euclid(nci) as usize,
                                (cz as isize + dz).rem_euclid(nci) as usize,
                            );
                            for &i in list {
                                for &j in &bins[flat(nb)] {
                                    if j <= i || i / na == j / na {
                                        continue;
                                    }
                                    let d = min_image(vsub(sys.xcart[i], sys.xcart[j]), l);
                                    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                                    if r < best {
                                        best = r;
                                    }
                                    let need = sys.radius[i] + sys.radius[j];
                                    let ov = need - r;
                                    if ov > max_ov {
                                        max_ov = ov;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    (best, max_ov.max(0.0))
}

const PEO_C_INF: F = 5.5;
const TET: F = 1.910_633_2; // 109.4712° in radians
const FLORY_R2_PER_M: F = 0.805; // ⟨R²⟩₀/M, Å²·mol/g (Fetters via Everaers 2020)

// ── synthetic all-atom PEO ─────────────────────────────────────────────────

fn vadd(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn vsub(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn vscale(a: [F; 3], s: F) -> [F; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn vnorm(a: [F; 3]) -> F {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn vunit(a: [F; 3]) -> [F; 3] {
    vscale(a, 1.0 / vnorm(a))
}
fn vcross(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// All-atom H-(CH₂-CH₂-O)ₙ-H in an all-trans zigzag: 7n + 2 atoms.
fn synthesize_peo(dp: usize) -> Frame {
    let alpha = (std::f64::consts::PI as F - TET) / 2.0;
    let nb = 3 * dp;
    let bond = |i: usize| -> F {
        // Bond i connects backbone atom i-1 → i: C-C = 1.53, C-O/O-C = 1.43.
        match i % 3 {
            1 => 1.53,
            _ => 1.43,
        }
    };
    let mut bb = Vec::with_capacity(nb);
    bb.push([0.0 as F, 0.0, 0.0]);
    for i in 1..nb {
        let b = bond(i);
        let dir = if i % 2 == 1 {
            [alpha.cos(), 0.0, alpha.sin()]
        } else {
            [alpha.cos(), 0.0, -alpha.sin()]
        };
        bb.push(vadd(bb[i - 1], vscale(dir, b)));
    }

    let mut pos: Vec<[F; 3]> = Vec::with_capacity(7 * dp + 2);
    let mut elem: Vec<&str> = Vec::with_capacity(7 * dp + 2);
    let mut bonds: Vec<(u32, u32)> = Vec::new();
    let mut bb_index = vec![0u32; nb];

    // Terminal H on the first carbon, pointing back along -x.
    pos.push(vadd(bb[0], [-1.10, 0.0, 0.0]));
    elem.push("H");

    let beta = (109.47_f64.to_radians() as F) / 2.0;
    for (i, &p) in bb.iter().enumerate() {
        let idx = pos.len() as u32;
        bb_index[i] = idx;
        let is_o = i % 3 == 2;
        pos.push(p);
        elem.push(if is_o { "O" } else { "C" });
        if i == 0 {
            bonds.push((0, idx));
        } else {
            bonds.push((bb_index[i - 1], idx));
        }
        if !is_o {
            // Two hydrogens, tetrahedral off the backbone plane.
            let prev = if i == 0 {
                vsub(p, [1.0, 0.0, 0.0])
            } else {
                bb[i - 1]
            };
            let next = if i + 1 < nb {
                bb[i + 1]
            } else {
                vadd(p, [1.0, 0.0, 0.0])
            };
            let d1 = vunit(vsub(prev, p));
            let d2 = vunit(vsub(next, p));
            let bis = vunit(vscale(vadd(d1, d2), -1.0));
            let n = vunit(vcross(d1, d2));
            for sgn in [1.0 as F, -1.0] {
                let dir = vunit(vadd(vscale(bis, beta.cos()), vscale(n, sgn * beta.sin())));
                let h = pos.len() as u32;
                pos.push(vadd(p, vscale(dir, 1.10)));
                elem.push("H");
                bonds.push((idx, h));
            }
        }
    }
    // Terminal H on the last oxygen.
    let last_o = bb_index[nb - 1];
    let h = pos.len() as u32;
    pos.push(vadd(bb[nb - 1], [1.0, 0.0, 0.0]));
    elem.push("H");
    bonds.push((last_o, h));

    let n = pos.len();
    let mut atoms = Block::new();
    atoms
        .insert(
            "x",
            Array1::from_vec(pos.iter().map(|p| p[0]).collect()).into_dyn(),
        )
        .unwrap();
    atoms
        .insert(
            "y",
            Array1::from_vec(pos.iter().map(|p| p[1]).collect()).into_dyn(),
        )
        .unwrap();
    atoms
        .insert(
            "z",
            Array1::from_vec(pos.iter().map(|p| p[2]).collect()).into_dyn(),
        )
        .unwrap();
    atoms
        .insert(
            "element",
            Array1::from_vec(elem.iter().map(|s| s.to_string()).collect()).into_dyn(),
        )
        .unwrap();
    let mut bblock = Block::new();
    bblock
        .insert(
            "atomi",
            Array1::from_vec(bonds.iter().map(|b| b.0).collect()).into_dyn(),
        )
        .unwrap();
    bblock
        .insert(
            "atomj",
            Array1::from_vec(bonds.iter().map(|b| b.1).collect()).into_dyn(),
        )
        .unwrap();
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame.insert("bonds", bblock);
    assert_eq!(n, 7 * dp + 2, "PEO synthesis atom count");
    frame
}

// ── report metrics ─────────────────────────────────────────────────────────

fn min_image(d: [F; 3], l: F) -> [F; 3] {
    std::array::from_fn(|k| d[k] - (d[k] / l).round() * l)
}

/// Bond-graph unwrap of one chain under the minimum image.
fn unwrap_chain(xyz: &[[F; 3]], bonds: &[(usize, usize)], l: F) -> Vec<[F; 3]> {
    let n = xyz.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(i, j) in bonds {
        adj[i].push(j);
        adj[j].push(i);
    }
    let mut un = vec![[F::NAN; 3]; n];
    let mut seen = vec![false; n];
    un[0] = xyz[0];
    seen[0] = true;
    let mut queue = vec![0usize];
    while let Some(a) = queue.pop() {
        for &b in &adj[a] {
            if !seen[b] {
                let d = min_image(vsub(xyz[b], un[a]), l);
                un[b] = vadd(un[a], d);
                seen[b] = true;
                queue.push(b);
            }
        }
    }
    un
}

fn radius_of_gyration(un: &[[F; 3]]) -> F {
    let n = un.len() as F;
    let mut c = [0.0 as F; 3];
    for p in un {
        for k in 0..3 {
            c[k] += p[k];
        }
    }
    for v in c.iter_mut() {
        *v /= n;
    }
    let mut s = 0.0;
    for p in un {
        for k in 0..3 {
            s += (p[k] - c[k]).powi(2);
        }
    }
    (s / n).sqrt()
}

/// Minimum inter-molecular distance via a uniform grid (minimum image).
fn min_inter_distance(pos: &[[F; 3]], na: usize, l: F) -> F {
    let cut = 4.0 as F;
    let nc = ((l / cut).floor() as usize).max(1);
    let cell = |p: [F; 3]| -> (usize, usize, usize) {
        let f = |v: F| (((v / l).rem_euclid(1.0) * nc as F) as usize).min(nc - 1);
        (f(p[0]), f(p[1]), f(p[2]))
    };
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); nc * nc * nc];
    let flat = |c: (usize, usize, usize)| (c.0 * nc + c.1) * nc + c.2;
    for (i, &p) in pos.iter().enumerate() {
        bins[flat(cell(p))].push(i);
    }
    let mut best = F::INFINITY;
    let nci = nc as isize;
    for cx in 0..nc {
        for cy in 0..nc {
            for cz in 0..nc {
                let list = &bins[flat((cx, cy, cz))];
                for dx in -1..=1isize {
                    for dy in -1..=1isize {
                        for dz in -1..=1isize {
                            let nb = (
                                (cx as isize + dx).rem_euclid(nci) as usize,
                                (cy as isize + dy).rem_euclid(nci) as usize,
                                (cz as isize + dz).rem_euclid(nci) as usize,
                            );
                            for &i in list {
                                for &j in &bins[flat(nb)] {
                                    if j <= i || i / na == j / na {
                                        continue;
                                    }
                                    let d = min_image(vsub(pos[i], pos[j]), l);
                                    let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                                    if r2 < best {
                                        best = r2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    best.sqrt()
}

/// Density-fluctuation merit: Var(n)/⟨n⟩ over random spheres of radius `d`
/// (ideal-gas/random placement ⇒ ≈ 1; smaller = more homogeneous).
fn density_fluctuation(pos: &[[F; 3]], l: F, d: F, seed: u64) -> F {
    let mut state = seed | 1;
    let mut rng = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as F / (1u64 << 53) as F
    };
    let m = 400;
    let d2 = d * d;
    let (mut s1, mut s2) = (0.0 as F, 0.0 as F);
    for _ in 0..m {
        let c = [rng() * l, rng() * l, rng() * l];
        let mut count = 0usize;
        for &p in pos {
            let dd = min_image(vsub(p, c), l);
            if dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2] < d2 {
                count += 1;
            }
        }
        s1 += count as F;
        s2 += (count * count) as F;
    }
    let mean = s1 / m as F;
    let var = s2 / m as F - mean * mean;
    if mean > 0.0 { var / mean } else { 0.0 }
}

#[allow(clippy::too_many_arguments)]
fn report(
    result: &State,
    na: usize,
    n_chains: usize,
    bonds: &[(usize, usize)],
    l: F,
    rg_ideal: F,
    elapsed: F,
) {
    let pos = result.positions();
    println!("── result ─────────────────────────────────────");
    println!(
        "  elapsed      : {elapsed:.3} s   converged={}  fdist={:.4e}  frest={:.4e}  degraded={}",
        result.converged, result.fdist, result.frest, result.degraded
    );
    println!(
        "  intra        : scored {:.4} Å   exempted {:.4} Å",
        result.intra.scored, result.intra.exempted
    );
    if let Some(elem) = result
        .frame
        .get("atoms")
        .and_then(|a| a.get_string("element"))
    {
        let names: Vec<String> = elem.iter().take(na).map(|s| s.to_string()).collect();
        print_closest_scored(&pos[..na], bonds, &names, na);
    }
    let unwrapped: Vec<Vec<[F; 3]>> = (0..n_chains)
        .map(|c| unwrap_chain(&pos[c * na..(c + 1) * na], bonds, l))
        .collect();
    let rgs: Vec<F> = unwrapped.iter().map(|u| radius_of_gyration(u)).collect();
    let rg_mean = rgs.iter().sum::<F>() / rgs.len() as F;
    let rg_min = rgs.iter().cloned().fold(F::INFINITY, F::min);
    let rg_max = rgs.iter().cloned().fold(0.0, F::max);
    println!(
        "  chain Rg     : mean {rg_mean:.2} Å  (min {rg_min:.2}, max {rg_max:.2})   \
         Flory ideal {rg_ideal:.1} Å  → deviation {:+.1}%  (ac-006)",
        100.0 * (rg_mean - rg_ideal) / rg_ideal
    );
    let k = rgs.len().min(5);
    let front: F = rgs.iter().take(k).sum::<F>() / k as F;
    let back: F = rgs.iter().rev().take(k).sum::<F>() / k as F;
    println!(
        "  chain-order bias: first-{k} Rg {front:.2} vs last-{k} Rg {back:.2}  (Δ {:.1}%)",
        100.0 * (front - back).abs() / rg_mean
    );
    println!("  internal distances ⟨R²(s)⟩/s (atom separations along the chain):");
    for &s in &[10usize, 30, 100, 300, 700] {
        let mut acc = 0.0;
        let mut cnt = 0usize;
        for u in &unwrapped {
            if s + 1 < u.len() {
                for i in (0..u.len() - s).step_by(s.max(1)) {
                    let d = vsub(u[i + s], u[i]);
                    acc += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                    cnt += 1;
                }
            }
        }
        if cnt > 0 {
            println!("    s = {s:>4}: {:.3} Å²", acc / cnt as F / s as F);
        }
    }
    println!(
        "  density fluctuation Var(n)/⟨n⟩: d=2 Å: {:.3}   d=4 Å: {:.3}   (random ≈ 1)",
        density_fluctuation(&pos, l, 2.0, 9),
        density_fluctuation(&pos, l, 4.0, 11)
    );
    if pos.len() <= 60_000 {
        println!(
            "  min inter-molecular distance (PBC): {:.3} Å",
            min_inter_distance(&pos, na, l)
        );
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("grow");
    let dp: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(200);
    let n_chains: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(25);
    // Evaluation claim range / claim cap: ρ = 1.2 g/cm³ for PEO melt
    // benchmarks. Not a runtime clamp of `with_density`.
    let density: F = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let seed: u64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(42);
    // Optional knobs: max_loops, trials, retract, relax_every, relax_window.
    // Serial needs max_loops ≥ n_chains (one pass per chain); the default
    // below is applied after `serial` is known.
    let max_loops_arg: Option<usize> = args.get(6).and_then(|s| s.parse().ok());
    let trials: usize = args.get(7).and_then(|s| s.parse().ok()).unwrap_or(12);
    let retract: usize = args.get(8).and_then(|s| s.parse().ok()).unwrap_or(10);
    let relax_every: usize = args.get(9).and_then(|s| s.parse().ok()).unwrap_or(25);
    let relax_window: usize = args.get(10).and_then(|s| s.parse().ok()).unwrap_or(6);
    let template_path = args.get(11);

    let frame = match template_path {
        Some(p) => molrs::io::data::pdb::read_pdb_frame(p)?,
        None => synthesize_peo(dp),
    };
    let na = frame.get("atoms").and_then(|a| a.nrows()).unwrap_or(0);
    let bonds: Vec<(usize, usize)> = {
        let b = frame.get("bonds").expect("template bonds");
        let i = b.get_uint("atomi").expect("bonds.atomi");
        let j = b.get_uint("atomj").expect("bonds.atomj");
        i.iter()
            .zip(j.iter())
            .map(|(&a, &c)| (a as usize, c as usize))
            .collect()
    };
    let m_chain: F = frame
        .get("atoms")
        .and_then(|a| a.get_string("element"))
        .map(|e| {
            e.iter()
                .filter_map(|s| {
                    use std::str::FromStr;
                    molpack::Element::from_str(s.trim())
                        .ok()
                        .map(|el| el.atomic_mass() as F)
                })
                .sum()
        })
        .unwrap_or(0.0);
    let rg_ideal = (FLORY_R2_PER_M * m_chain / 6.0).sqrt();
    let total_g = n_chains as F * m_chain / 6.022_140_76e23;
    let l = (total_g / density * 1e24).cbrt();
    let tmpl_elem: Vec<String> = frame
        .get("atoms")
        .and_then(|a| a.get_string("element"))
        .map(|e| e.iter().map(|s| s.to_string()).collect())
        .unwrap_or_default();
    let tmpl_pos: Vec<[F; 3]> = {
        let a = frame.get("atoms").expect("atoms");
        let x = a.get_float("x").expect("x");
        let y = a.get_float("y").expect("y");
        let z = a.get_float("z").expect("z");
        x.iter()
            .zip(y.iter())
            .zip(z.iter())
            .map(|((&xi, &yi), &zi)| [xi as F, yi as F, zi as F])
            .collect()
    };

    println!("── system ─────────────────────────────────────");
    println!(
        "  template     : {} ({na} atoms, {} bonds)",
        template_path.map(String::as_str).unwrap_or("synthetic PEO"),
        bonds.len()
    );
    println!(
        "  chains       : {n_chains} × {m_chain:.1} amu   total atoms {}",
        na * n_chains
    );
    println!("  density      : {density} g/cm³  →  L = {l:.2} Å");
    println!("  Flory ideal Rg: {rg_ideal:.1} Å   mode: {mode}   seed: {seed}");

    let h_idx: Vec<usize> = frame
        .get("atoms")
        .and_then(|a| a.get_string("element"))
        .map(|e| {
            e.iter()
                .enumerate()
                .filter(|(_, s)| s.trim().eq_ignore_ascii_case("h"))
                .map(|(i, _)| i)
                .collect()
        })
        .unwrap_or_default();
    let h_radius = std::env::var("PEO_H_RADIUS")
        .ok()
        .map(|s| s.parse::<F>())
        .transpose()?;
    let tolerance: F = match std::env::var("PEO_TOLERANCE") {
        Ok(s) => s.parse()?,
        Err(_) => 2.0,
    };
    let void_bias = matches!(
        std::env::var("PEO_VOID_BIAS").as_deref(),
        Ok("1") | Ok("true") | Ok("TRUE")
    );
    // Default: one chain completes before the next starts. PEO_SERIAL=0
    // restores the crate default (every pending chain advances each round).
    let serial = match std::env::var("PEO_SERIAL") {
        Ok(s) => !matches!(s.as_str(), "0" | "false" | "FALSE" | "no"),
        Err(_) => false,
    };
    let max_loops: usize = max_loops_arg.unwrap_or(if serial {
        n_chains.saturating_mul(10).max(60)
    } else {
        60
    });
    let mut target = Target::new(frame, n_chains).with_name("peo");
    match h_radius {
        Some(r) if !h_idx.is_empty() => {
            target = target.with_atom_radius(&h_idx, r);
            println!(
                "  H radius     : {r} Å  ({} of {na} atoms; PEO_H_RADIUS override)",
                h_idx.len()
            );
        }
        _ => println!(
            "  packing R    : {:.3} Å for every atom (tolerance/2)",
            tolerance / 2.0
        ),
    }
    println!(
        "  tolerance    : {tolerance} Å  (contact = {tolerance} Å; min_hard_scale=1, no ladder)"
    );
    println!("  special_bonds: {:?}", target.special_bonds.as_slice());
    println!("  template (all-trans) scored pairs — chemically real geometry:");
    print_closest_scored(&tmpl_pos, &bonds, &tmpl_elem, na);

    let ladder = Arc::new(Mutex::new(Ladder {
        start: Instant::now(),
        rounds: 0,
        rad_min: 1.0,
        first_drop: None,
        last_hist: Vec::new(),
    }));
    let trace_path = std::env::var("PEO_TRACE").ok();

    let t0 = Instant::now();
    let result = match mode {
        "grow" => {
            let cfg = GrowConfig::new(TorsionPrior::three_state_from_c_inf(PEO_C_INF, TET))
                .with_trials(trials)
                .with_retract(retract)
                .with_relax(relax_every, relax_window)
                .with_void_bias(void_bias)
                .with_serial(serial)
                .with_min_hard_scale(1.0)
                .with_soften_after(1_000_000);
            println!(
                "  grow knobs   : trials={trials} retract={retract} relax=({relax_every},{relax_window}) void_bias={void_bias} serial={serial}"
            );
            println!(
                "  scheduling   : {}",
                if serial {
                    "1 chain at a time (GrowConfig::with_serial)"
                } else {
                    "all pending chains each round (crate default)"
                }
            );
            CbmcGrow::from_config(cfg)
                .with_seed(seed)
                .with_tolerance(tolerance)
                .with_density(density)
                .with_handler(Box::new(LadderHandler::new(
                    Arc::clone(&ladder),
                    trace_path.as_deref(),
                )))
                .run(&[target], max_loops)?
        }
        "rigid" => GenCanPack::new()
            .with_seed(seed)
            .with_tolerance(tolerance)
            .with_periodic_box([0.0; 3], [l, l, l], [true; 3])
            .run(&[target], max_loops)?,
        "auhl" => {
            // 1. Reduced-EV CBMC (all-parallel) then seeded GENCAN push-off.
            let grow_tol: F = std::env::var("PEO_GROW_TOL")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.6);
            let push_tol: F = 2.0;
            let h_push: F = h_radius.unwrap_or(0.05);
            println!(
                "  Auhl         : CbmcGrow all-parallel @ {grow_tol} Å → GenCanPack @ {push_tol} Å (H={h_push})"
            );
            let cfg = GrowConfig::new(TorsionPrior::three_state_from_c_inf(PEO_C_INF, TET))
                .with_trials(trials)
                .with_retract(retract)
                .with_relax(relax_every, relax_window)
                .with_serial(false)
                .with_min_hard_scale(1.0)
                .with_soften_after(1_000_000);
            let g0 = Instant::now();
            let grown = CbmcGrow::from_config(cfg)
                .with_seed(seed)
                .with_tolerance(grow_tol)
                .with_density(density)
                .run(&[target], max_loops)?;
            println!(
                "  grow         : {:.3} s  converged={}  fdist={:.4e}  degraded={}  intra scored {:.4} Å",
                g0.elapsed().as_secs_f64(),
                grown.converged,
                grown.fdist,
                grown.degraded,
                grown.intra.scored
            );
            let gp = grown.positions();
            println!(
                "  grow atom0   : [{:.4}, {:.4}, {:.4}]  (chain 0, not reshuffled if push-off)",
                gp[0][0], gp[0][1], gp[0][2]
            );
            let mut t_push = Target::new(synthesize_peo(dp), n_chains).with_name("peo");
            if !h_idx.is_empty() {
                t_push = t_push.with_atom_radius(&h_idx, h_push);
            }
            if matches!(
                std::env::var("PEO_PUSH").as_deref(),
                Ok("0") | Ok("false") | Ok("off")
            ) {
                grown
            } else {
                let p0 = Instant::now();
                let pushed = GenCanPack::new()
                    .with_restart(&grown)
                    .with_seed(seed)
                    .with_tolerance(push_tol)
                    .run(&[t_push], 30)?;
                let pp = pushed.positions();
                let d0 = vnorm(vsub(pp[0], gp[0]));
                println!(
                    "  push-off     : {:.3} s  converged={}  fdist={:.4e}  degraded={}  atom0 moved {:.3} Å",
                    p0.elapsed().as_secs_f64(),
                    pushed.converged,
                    pushed.fdist,
                    pushed.degraded,
                    d0
                );
                pushed
            }
        }
        "lattice" => {
            // 2. Diamond-lattice SAW then seeded GENCAN push-off.
            let push_tol: F = 2.0;
            let h_push: F = h_radius.unwrap_or(0.05);
            println!(
                "  lattice      : LatticeGrow occupancy-guard → GenCanPack @ {push_tol} Å (H={h_push})"
            );
            let g0 = Instant::now();
            let grown = LatticeGrow::new(TorsionPrior::three_state_from_c_inf(PEO_C_INF, TET))
                .with_seed(seed)
                .with_tolerance(push_tol)
                .with_density(density)
                .run(&[target], max_loops)?;
            println!(
                "  grow         : {:.3} s  converged={}  fdist={:.4e}  degraded={}  intra scored {:.4} Å",
                g0.elapsed().as_secs_f64(),
                grown.converged,
                grown.fdist,
                grown.degraded,
                grown.intra.scored
            );
            let mut t_push = Target::new(synthesize_peo(dp), n_chains).with_name("peo");
            if !h_idx.is_empty() {
                t_push = t_push.with_atom_radius(&h_idx, h_push);
            }
            if matches!(
                std::env::var("PEO_PUSH").as_deref(),
                Ok("0") | Ok("false") | Ok("off")
            ) {
                grown
            } else {
                let p0 = Instant::now();
                let pushed = GenCanPack::new()
                    .with_restart(&grown)
                    .with_seed(seed)
                    .with_tolerance(push_tol)
                    .run(&[t_push], 30)?;
                println!(
                    "  push-off     : {:.3} s  converged={}  fdist={:.4e}  degraded={}",
                    p0.elapsed().as_secs_f64(),
                    pushed.converged,
                    pushed.fdist,
                    pushed.degraded
                );
                pushed
            }
        }
        "grow1" => {
            // Grow one conformer, then rigid-pack N copies (not seeded push-off).
            let h_r: F = h_radius.unwrap_or(0.05);
            println!(
                "  grow1        : 1 chain (serial, H={h_r} Å, tol=2.0) → GenCanPack × {n_chains} @ {density} g/cm³"
            );
            let mut t1 = Target::new(synthesize_peo(dp), 1).with_name("peo");
            if !h_idx.is_empty() {
                t1 = t1.with_atom_radius(&h_idx, h_r);
            }
            let cfg = GrowConfig::new(TorsionPrior::three_state_from_c_inf(PEO_C_INF, TET))
                .with_trials(trials)
                .with_retract(retract)
                .with_relax(relax_every, relax_window)
                .with_serial(true)
                .with_min_hard_scale(1.0)
                .with_soften_after(1_000_000);
            let g0 = Instant::now();
            let grown = CbmcGrow::from_config(cfg)
                .with_seed(seed)
                .with_tolerance(2.0)
                .with_density(0.001)
                .run(&[t1], max_loops)?;
            println!(
                "  grow         : {:.3} s  converged={}  fdist={:.4e}  intra scored {:.4} Å",
                g0.elapsed().as_secs_f64(),
                grown.converged,
                grown.fdist,
                grown.intra.scored
            );
            let mut tn = Target::new(grown.frame.clone(), n_chains).with_name("peo");
            if !h_idx.is_empty() {
                tn = tn.with_atom_radius(&h_idx, h_r);
            }
            let p0 = Instant::now();
            let packed = GenCanPack::new()
                .with_seed(seed)
                .with_tolerance(2.0)
                .with_density(density)
                .run(&[tn], 60)?;
            println!(
                "  pack         : {:.3} s  converged={}  fdist={:.4e}",
                p0.elapsed().as_secs_f64(),
                packed.converged,
                packed.fdist
            );
            packed
        }
        other => {
            return Err(
                format!("unknown mode {other}: use grow | rigid | auhl | lattice | grow1").into(),
            );
        }
    };
    let elapsed = t0.elapsed().as_secs_f64();
    {
        let s = ladder.lock().expect("ladder");
        println!(
            "  ladder       : {} rounds   radscale min {:.3}{}",
            s.rounds,
            s.rad_min,
            match s.first_drop {
                Some((r, v)) => format!("   first drop at round {r} → {v:.3}"),
                None => "   no drop (hard core stayed at 1.0)".into(),
            }
        );
        if !s.last_hist.is_empty() {
            let n = s.last_hist.len();
            let done = s.last_hist.iter().filter(|&&k| k == na).count();
            let half = s
                .last_hist
                .iter()
                .filter(|&&k| k * 2 >= na && k < na)
                .count();
            let early = s.last_hist.iter().filter(|&&k| k * 2 < na).count();
            let mut v = s.last_hist.clone();
            v.sort_unstable();
            println!(
                "  chain stage  : done={done}/{n}  mid-unfinished={half}  early={early}  \
                 placed min/med/max={}/{}/{}",
                v.first().copied().unwrap_or(0),
                v[v.len() / 2],
                v.last().copied().unwrap_or(0)
            );
            if let Some(p) = trace_path.as_ref() {
                let hist = format!("{p}.chains");
                if let Ok(mut f) = std::fs::File::create(&hist) {
                    use std::io::Write;
                    for (i, k) in s.last_hist.iter().enumerate() {
                        let _ = writeln!(f, "{i},{k},{na}");
                    }
                    println!("  wrote        : {p}  and  {hist}");
                }
            }
        }
    }
    report(&result, na, n_chains, &bonds, l, rg_ideal, elapsed);
    Ok(())
}
