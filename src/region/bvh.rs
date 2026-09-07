//! Bounding-volume hierarchy over a triangle soup.
//!
//! [`StlRegion`](super::StlRegion) asks two questions of its mesh — closest
//! point, and how many times a ray crosses it — and both were a scan over
//! every triangle. That is fine for a twelve-triangle box and hopeless at melt
//! scale, where the lattice mask alone asks about a million sites: the scan,
//! not the packing, becomes the run.
//!
//! One tree serves both queries. The closest-point descent prunes on the
//! squared distance to a node box and visits the nearer child first; the ray
//! walk prunes on a slab test. Both are exact — the BVH changes what is
//! *visited*, never what is *answered*, which is what makes it swappable under
//! the existing goldens.

use molrs::types::F;

use super::stl::{dot, sub};

/// Triangles per leaf. A leaf scan is a handful of closest-point evaluations,
/// cheaper than the box tests that would separate them further.
const LEAF_SIZE: usize = 4;

#[derive(Debug, Clone, Default)]
struct Node {
    min: [F; 3],
    max: [F; 3],
    /// Interior node: index of the left child, with the right at `first + 1`.
    /// Leaf: index into [`Bvh::order`] of its first triangle.
    first: u32,
    /// Triangle count for a leaf; `0` marks an interior node.
    count: u32,
}

/// Median-split BVH. Empty only for an empty soup, which the region rejects
/// before it gets here.
#[derive(Debug, Clone)]
pub(super) struct Bvh {
    nodes: Vec<Node>,
    /// Triangle indices permuted so that every leaf owns a contiguous run.
    order: Vec<u32>,
}

impl Bvh {
    pub(super) fn build(triangles: &[[[F; 3]; 3]]) -> Self {
        let mut order: Vec<u32> = (0..triangles.len() as u32).collect();
        let centroids: Vec<[F; 3]> = triangles.iter().map(centroid).collect();
        let mut nodes = vec![Node::default()];
        if !order.is_empty() {
            let n = order.len();
            build_node(&mut nodes, &mut order, triangles, &centroids, 0, 0, n);
        }
        Self { nodes, order }
    }

    /// Closest point on the soup to `p`, as `(distance², point)`.
    ///
    /// `closest_on` returns the closest point on one triangle; the tree only
    /// decides which triangles are worth asking about.
    pub(super) fn nearest<G>(&self, p: &[F; 3], mut closest_on: G) -> (F, [F; 3])
    where
        G: FnMut(u32) -> [F; 3],
    {
        let mut best = (F::INFINITY, *p);
        if !self.order.is_empty() {
            self.nearest_in(0, p, &mut closest_on, &mut best);
        }
        best
    }

    fn nearest_in<G>(&self, n: usize, p: &[F; 3], closest_on: &mut G, best: &mut (F, [F; 3]))
    where
        G: FnMut(u32) -> [F; 3],
    {
        let node = &self.nodes[n];
        if node.count > 0 {
            let end = (node.first + node.count) as usize;
            for &t in &self.order[node.first as usize..end] {
                let c = closest_on(t);
                let d2 = dist2(*p, c);
                if d2 < best.0 {
                    *best = (d2, c);
                }
            }
            return;
        }
        let left = node.first as usize;
        let dl = box_dist2(&self.nodes[left], p);
        let dr = box_dist2(&self.nodes[left + 1], p);
        // Nearer child first: it tightens `best` before the other is tested,
        // which is what turns the second test into a prune.
        let (near, far, dfar) = if dl <= dr {
            (left, left + 1, dr)
        } else {
            (left + 1, left, dl)
        };
        self.nearest_in(near, p, closest_on, best);
        if dfar < best.0 {
            self.nearest_in(far, p, closest_on, best);
        }
    }

    /// Whether any triangle comes within `radius` of `p`.
    ///
    /// The same descent as [`nearest`](Self::nearest) without the bookkeeping:
    /// a membership test only needs to know whether the surface is within the
    /// boundary tolerance, and answering that as a radius query stops at the
    /// first hit instead of finding the closest one.
    pub(super) fn any_within<G>(&self, p: &[F; 3], radius: F, mut closest_on: G) -> bool
    where
        G: FnMut(u32) -> [F; 3],
    {
        !self.order.is_empty() && self.any_within_in(0, p, radius * radius, &mut closest_on)
    }

    fn any_within_in<G>(&self, n: usize, p: &[F; 3], r2: F, closest_on: &mut G) -> bool
    where
        G: FnMut(u32) -> [F; 3],
    {
        let node = &self.nodes[n];
        if box_dist2(node, p) > r2 {
            return false;
        }
        if node.count > 0 {
            let end = (node.first + node.count) as usize;
            return self.order[node.first as usize..end]
                .iter()
                .any(|&t| dist2(*p, closest_on(t)) < r2);
        }
        let left = node.first as usize;
        self.any_within_in(left, p, r2, closest_on)
            || self.any_within_in(left + 1, p, r2, closest_on)
    }

    /// How many triangles `hit` reports along the ray from `origin`.
    ///
    /// `inv` is the componentwise reciprocal of the direction; the caller owns
    /// the direction, so it also owns keeping it free of zero components.
    pub(super) fn count_hits<H>(&self, origin: &[F; 3], inv: &[F; 3], mut hit: H) -> u32
    where
        H: FnMut(u32) -> bool,
    {
        if self.order.is_empty() {
            return 0;
        }
        self.count_hits_in(0, origin, inv, &mut hit)
    }

    fn count_hits_in<H>(&self, n: usize, origin: &[F; 3], inv: &[F; 3], hit: &mut H) -> u32
    where
        H: FnMut(u32) -> bool,
    {
        let node = &self.nodes[n];
        if !ray_hits_box(node, origin, inv) {
            return 0;
        }
        if node.count > 0 {
            let end = (node.first + node.count) as usize;
            return self.order[node.first as usize..end]
                .iter()
                .filter(|&&t| hit(t))
                .count() as u32;
        }
        let left = node.first as usize;
        self.count_hits_in(left, origin, inv, hit) + self.count_hits_in(left + 1, origin, inv, hit)
    }
}

/// Fill node `n` with the bounds of `order[start..start + count]`, splitting
/// until a leaf is small enough.
fn build_node(
    nodes: &mut Vec<Node>,
    order: &mut Vec<u32>,
    triangles: &[[[F; 3]; 3]],
    centroids: &[[F; 3]],
    n: usize,
    start: usize,
    count: usize,
) {
    let (min, max) = bounds(triangles, &order[start..start + count]);
    nodes[n].min = min;
    nodes[n].max = max;
    if count <= LEAF_SIZE {
        nodes[n].first = start as u32;
        nodes[n].count = count as u32;
        return;
    }
    // Median split on the widest spread of centroids: no surface-area
    // heuristic, because a mesh from a mesher is already spatially coherent
    // and the extra build cost buys nothing measurable here.
    let axis = widest_axis(centroids, &order[start..start + count]);
    let half = count / 2;
    order[start..start + count].select_nth_unstable_by(half, |&a, &b| {
        centroids[a as usize][axis].total_cmp(&centroids[b as usize][axis])
    });
    let left = nodes.len();
    nodes.push(Node::default());
    nodes.push(Node::default());
    nodes[n].first = left as u32;
    nodes[n].count = 0;
    build_node(nodes, order, triangles, centroids, left, start, half);
    build_node(
        nodes,
        order,
        triangles,
        centroids,
        left + 1,
        start + half,
        count - half,
    );
}

fn centroid(t: &[[F; 3]; 3]) -> [F; 3] {
    [
        (t[0][0] + t[1][0] + t[2][0]) / 3.0,
        (t[0][1] + t[1][1] + t[2][1]) / 3.0,
        (t[0][2] + t[1][2] + t[2][2]) / 3.0,
    ]
}

fn bounds(triangles: &[[[F; 3]; 3]], idx: &[u32]) -> ([F; 3], [F; 3]) {
    let mut min = [F::INFINITY; 3];
    let mut max = [F::NEG_INFINITY; 3];
    for &i in idx {
        for v in &triangles[i as usize] {
            for k in 0..3 {
                min[k] = min[k].min(v[k]);
                max[k] = max[k].max(v[k]);
            }
        }
    }
    (min, max)
}

fn widest_axis(centroids: &[[F; 3]], idx: &[u32]) -> usize {
    let mut min = [F::INFINITY; 3];
    let mut max = [F::NEG_INFINITY; 3];
    for &i in idx {
        let c = centroids[i as usize];
        for k in 0..3 {
            min[k] = min[k].min(c[k]);
            max[k] = max[k].max(c[k]);
        }
    }
    let mut axis = 0;
    let mut span = max[0] - min[0];
    for (k, s) in (1..3).map(|k| (k, max[k] - min[k])) {
        if s > span {
            span = s;
            axis = k;
        }
    }
    axis
}

/// Squared distance from `p` to a node box; `0` when `p` is inside it.
fn box_dist2(node: &Node, p: &[F; 3]) -> F {
    let mut d2 = 0.0;
    for (k, &pk) in p.iter().enumerate() {
        let v = if pk < node.min[k] {
            node.min[k] - pk
        } else if pk > node.max[k] {
            pk - node.max[k]
        } else {
            continue;
        };
        d2 += v * v;
    }
    d2
}

/// Slab test for a ray that starts at `origin` and never turns back.
fn ray_hits_box(node: &Node, origin: &[F; 3], inv: &[F; 3]) -> bool {
    let mut near: F = 0.0;
    let mut far = F::INFINITY;
    for (k, (&o, &iv)) in origin.iter().zip(inv.iter()).enumerate() {
        let a = (node.min[k] - o) * iv;
        let b = (node.max[k] - o) * iv;
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        near = near.max(lo);
        far = far.min(hi);
        if far < near {
            return false;
        }
    }
    true
}

fn dist2(a: [F; 3], b: [F; 3]) -> F {
    let d = sub(a, b);
    dot(d, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic triangle soup — no rand dependency in a unit test.
    fn soup(n: usize) -> Vec<[[F; 3]; 3]> {
        let mut s: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as F / (1u64 << 53) as F
        };
        (0..n)
            .map(|_| {
                let o = [next() * 20.0, next() * 20.0, next() * 20.0];
                let mut v = || {
                    [
                        o[0] + next() * 3.0,
                        o[1] + next() * 3.0,
                        o[2] + next() * 3.0,
                    ]
                };
                [v(), v(), v()]
            })
            .collect()
    }

    /// Closest point on a triangle, by the same construction the region uses —
    /// duplicated here so the BVH test does not depend on the region.
    fn closest(p: [F; 3], t: &[[F; 3]; 3]) -> [F; 3] {
        let mut best = t[0];
        let mut best_d = dist2(p, t[0]);
        // Dense sampling of the triangle: crude, but this test only needs a
        // reference the tree must reproduce, not a fast one.
        let n = 24;
        for i in 0..=n {
            for j in 0..=(n - i) {
                let (u, v) = (i as F / n as F, j as F / n as F);
                let w = 1.0 - u - v;
                let q = [
                    u * t[0][0] + v * t[1][0] + w * t[2][0],
                    u * t[0][1] + v * t[1][1] + w * t[2][1],
                    u * t[0][2] + v * t[1][2] + w * t[2][2],
                ];
                let d = dist2(p, q);
                if d < best_d {
                    best_d = d;
                    best = q;
                }
            }
        }
        best
    }

    #[test]
    fn nearest_matches_the_linear_scan() {
        let tris = soup(300);
        let bvh = Bvh::build(&tris);
        for p in [
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 10.0],
            [-5.0, 22.0, 3.0],
            [21.5, 21.5, 21.5],
        ] {
            let (d2, _) = bvh.nearest(&p, |i| closest(p, &tris[i as usize]));
            let brute = tris
                .iter()
                .map(|t| dist2(p, closest(p, t)))
                .fold(F::INFINITY, F::min);
            assert!(
                (d2 - brute).abs() < 1e-12,
                "bvh {d2} vs scan {brute} at {p:?}"
            );
        }
    }

    /// Möller–Trumbore, forward hits only — the same shape of test the region
    /// runs, so the tree may only skip what this would have rejected anyway.
    fn ray_hits(o: [F; 3], d: [F; 3], t: &[[F; 3]; 3]) -> bool {
        let sub = |a: [F; 3], b: [F; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let cross = |a: [F; 3], b: [F; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let dot = |a: [F; 3], b: [F; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let (e1, e2) = (sub(t[1], t[0]), sub(t[2], t[0]));
        let h = cross(d, e2);
        let a = dot(e1, h);
        if a.abs() < 1e-12 {
            return false;
        }
        let f = 1.0 / a;
        let s = sub(o, t[0]);
        let u = f * dot(s, h);
        if !(0.0..=1.0).contains(&u) {
            return false;
        }
        let q = cross(s, e1);
        let v = f * dot(d, q);
        if v < 0.0 || u + v > 1.0 {
            return false;
        }
        f * dot(e2, q) > 1e-12
    }

    #[test]
    fn count_hits_matches_the_linear_scan() {
        let tris = soup(300);
        let bvh = Bvh::build(&tris);
        let dir = {
            let d = [1.0, std::f64::consts::SQRT_2, std::f64::consts::PI];
            let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            [d[0] / n, d[1] / n, d[2] / n]
        };
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        for origin in [
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 10.0],
            [-3.0, 7.5, 12.0],
            [19.0, 2.0, 8.0],
        ] {
            let counted =
                bvh.count_hits(&origin, &inv, |i| ray_hits(origin, dir, &tris[i as usize]));
            let brute = tris.iter().filter(|t| ray_hits(origin, dir, t)).count() as u32;
            assert_eq!(counted, brute, "at {origin:?}");
        }
    }

    #[test]
    fn a_single_triangle_is_one_leaf() {
        let tris = vec![[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]];
        let bvh = Bvh::build(&tris);
        let (d2, c) = bvh.nearest(&[0.0, 0.0, 2.0], |i| {
            closest([0.0, 0.0, 2.0], &tris[i as usize])
        });
        assert!((d2 - 4.0).abs() < 1e-12);
        assert!(c[2].abs() < 1e-12);
    }
}
