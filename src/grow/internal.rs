//! Internal-coordinate decomposition of a template molecule.
//!
//! Growth needs the template's *chemistry*, not its shape. This module walks
//! the template's bond graph breadth-first from a chain end and rewrites every
//! atom as `(bond, angle, torsion)` against three already-placed reference
//! atoms. Bond lengths, bond angles and non-rotatable torsions are copied
//! verbatim, so the rebuilt molecule carries the template's bonded geometry
//! exactly — a force field typified against the template stays valid. Only the
//! rotatable torsions become free variables, and those are what the packer
//! samples.
//!
//! The atoms are grouped into **steps**: one step per free torsion, holding the
//! representative atom of that torsion plus every following atom whose own
//! torsion is already determined. A step is the unit the growth driver places
//! and scores at once.

use std::collections::HashMap;

use molrs::store::frame::Frame;
use molrs::system::atomistic::Atomistic;
use molrs::system::bond::BondType;
use molrs::types::F;

use crate::grow::GrowError;
use crate::topology::{Topology, TopologyError};

/// One atom's placement recipe against three earlier atoms.
#[derive(Debug, Clone)]
pub struct Site {
    /// Template atom index this site places.
    pub atom: usize,
    /// Reference atoms `(p, g, gg)`, all placed before this site.
    pub refs: [usize; 3],
    /// `|atom - p|`, copied from the template.
    pub bond: F,
    /// `angle(g, p, atom)`, copied from the template.
    pub angle: F,
    /// `dihedral(gg, g, p, atom)` in the template.
    pub torsion: F,
    /// Free-torsion variable this site follows, with the offset to add.
    /// `None` ⇒ the torsion is rigid and `torsion` is used as-is.
    pub follows: Option<(usize, F)>,
}

/// A template molecule rewritten as a growable internal-coordinate tree.
#[derive(Debug, Clone)]
pub struct InternalTree {
    n_atoms: usize,
    /// The three atoms placed as one rigid seed, in placement order.
    seed: [usize; 3],
    /// Seed geometry in the template frame, centred on `seed[0]`.
    seed_local: [[F; 3]; 3],
    /// Placement recipes, in BFS order (excludes the seed).
    sites: Vec<Site>,
    /// `steps[k]` is the half-open range of `sites` placed by step `k`.
    steps: Vec<(usize, usize)>,
    /// Free-torsion variable owned by each step (`None` for the first step,
    /// whose sites are all rigid relative to the seed).
    step_var: Vec<Option<usize>>,
    n_vars: usize,
    /// Per-atom excluded partners: same-molecule atoms within the exclusion
    /// depth (in bonds) this tree was built with; they are governed by the
    /// template geometry / torsion prior and must not be scored. Each list is
    /// ascending and has its own root included, so it is the overlap field's
    /// skip set verbatim.
    exclusions: Vec<Vec<u32>>,
}

/// Default bond-graph exclusion depth: 3 bonds, i.e. 1-2 / 1-3 / 1-4 — the
/// all-atom convention. 1-2 and 1-3 distances are fixed by the template's
/// bond lengths and angles. 1-4 distances are **not** fixed — they swing with
/// the sampled torsion (butane: trans ≈ 3.9 Å vs cis ≈ 2.9 Å) — but they are
/// governed by the torsion prior, and a hard-core check would wrongly reject
/// legitimate gauche/cis conformers, so they are excluded too (the same
/// division of labor as force-field 1-4 scaling). 1-5 and beyond are what a
/// chain must not thread through itself, so they are scored.
///
/// Coarse-grained templates conventionally use a shallower depth (1-2 or
/// 1-3); it is a per-target parameter (`GrowConfig::exclusion_depth`), not a
/// constant — this default serves the all-atom case.
const DEFAULT_EXCLUDE_BONDS: usize = 3;

impl InternalTree {
    /// Decompose a template frame with the all-atom default exclusion depth
    /// (`DEFAULT_EXCLUDE_BONDS` = 3, i.e. 1-2/1-3/1-4). The frame must carry
    /// an `atoms` block with `x`/`y`/`z` and a `bonds` block with
    /// `atomi`/`atomj`.
    pub fn from_frame(frame: &Frame) -> Result<Self, GrowError> {
        Self::from_frame_with_depth(frame, DEFAULT_EXCLUDE_BONDS)
    }

    /// Decompose a template frame with an explicit intramolecular exclusion
    /// depth in bonds (`3` = exclude 1-2/1-3/1-4; CG templates typically use
    /// `1` or `2`).
    pub fn from_frame_with_depth(frame: &Frame, exclude_bonds: usize) -> Result<Self, GrowError> {
        let (topo, xyz) =
            Topology::from_frame_with_positions(frame).map_err(GrowError::Topology)?;
        let n = topo.natoms();
        if n < 3 {
            return Err(GrowError::TemplateTooSmall(n));
        }
        topo.require_connected().map_err(GrowError::Topology)?;

        let root = diameter_endpoint(&topo);
        let (order, parent) = bfs_order(&topo, root, n)?;

        // Rotatable bonds, as an unordered key set. PDB / GROMACS / XYZ
        // connectivity reads back as `BondType::Unknown`; for a conformer
        // search an unclassed bond is a rotatable single bond.
        let graph = Atomistic::from_frame(frame).map_err(|e| GrowError::Perceive(e.to_string()))?;
        let rotatable = rotatable_bond_keys(&graph);

        let seed = [order[0], order[1], order[2]];
        let seed_local = {
            let o = xyz[seed[0]];
            [[0.0, 0.0, 0.0], sub(xyz[seed[1]], o), sub(xyz[seed[2]], o)]
        };

        let mut placed = vec![false; n];
        for a in seed {
            placed[a] = true;
        }

        let mut sites: Vec<Site> = Vec::with_capacity(n - 3);
        // (p, g) key → (variable index, representative torsion value)
        let mut var_of_bond: HashMap<(usize, usize), (usize, F)> = HashMap::new();
        let mut n_vars = 0usize;

        for &i in &order[3..] {
            let p = parent[i].expect("BFS parent");
            let g = pick_ref(&topo, &placed, p, i, parent[p]).ok_or(GrowError::NoReference(i))?;
            let gg = pick_ref(&topo, &placed, g, p, parent[g])
                .filter(|&c| c != i)
                .or_else(|| {
                    // Fallback reference: pick the placed atom whose frame is
                    // best conditioned — a gg (near-)collinear with g→p
                    // degenerates NeRF's frame and corrupts the placement
                    // silently, so maximize the cross-product norm instead of
                    // taking the first index.
                    let gp = sub(xyz[p], xyz[g]);
                    (0..n)
                        .filter(|&c| placed[c] && c != p && c != g)
                        .max_by(|&a, &b| {
                            let na = norm(cross(sub(xyz[g], xyz[a]), gp));
                            let nb = norm(cross(sub(xyz[g], xyz[b]), gp));
                            na.total_cmp(&nb)
                        })
                })
                .ok_or(GrowError::NoReference(i))?;

            let bond = norm(sub(xyz[i], xyz[p]));
            let angle = angle_at(xyz[g], xyz[p], xyz[i]);
            let torsion = dihedral(xyz[gg], xyz[g], xyz[p], xyz[i]);

            let key = (p.min(g), p.max(g));
            let follows = if rotatable.contains(&key) {
                match var_of_bond.get(&key) {
                    // A later child of the same rotatable bond: it keeps its
                    // template offset from the representative, so the local
                    // geometry (and chirality) around `p` is preserved exactly.
                    Some(&(v, rep_torsion)) => Some((v, wrap_pi(torsion - rep_torsion))),
                    None => {
                        let v = n_vars;
                        n_vars += 1;
                        var_of_bond.insert(key, (v, torsion));
                        Some((v, 0.0))
                    }
                }
            } else {
                None
            };

            sites.push(Site {
                atom: i,
                refs: [p, g, gg],
                bond,
                angle,
                torsion,
                follows,
            });
            placed[i] = true;
        }

        // Group into steps: a new step opens at each site that *introduces* a
        // variable (offset 0.0 and a variable index not yet seen).
        let mut steps: Vec<(usize, usize)> = Vec::new();
        let mut step_var: Vec<Option<usize>> = Vec::new();
        let mut seen_var = vec![false; n_vars];
        let mut start = 0usize;
        let mut current: Option<usize> = None;
        for (idx, s) in sites.iter().enumerate() {
            let opens = match s.follows {
                Some((v, _)) if !seen_var[v] => {
                    seen_var[v] = true;
                    Some(v)
                }
                _ => None,
            };
            if let Some(v) = opens {
                if idx > start || current.is_some() {
                    steps.push((start, idx));
                    step_var.push(current);
                }
                start = idx;
                current = Some(v);
            }
        }
        steps.push((start, sites.len()));
        step_var.push(current);

        let exclusions = topo.exclusions(exclude_bonds);

        Ok(Self {
            n_atoms: n,
            seed,
            seed_local,
            sites,
            steps,
            step_var,
            n_vars,
            exclusions,
        })
    }

    pub fn n_atoms(&self) -> usize {
        self.n_atoms
    }
    pub fn n_steps(&self) -> usize {
        self.steps.len()
    }
    pub fn n_vars(&self) -> usize {
        self.n_vars
    }
    pub fn seed_atoms(&self) -> [usize; 3] {
        self.seed
    }
    /// Free-torsion variable sampled at step `k`, if any.
    pub fn step_var(&self, k: usize) -> Option<usize> {
        self.step_var[k]
    }
    /// Same-molecule partners of `atom` within this tree's exclusion depth,
    /// ascending and **including `atom` itself** — the skip set the overlap
    /// field expects, since an atom must not be scored against its own
    /// position. Produced by [`Topology::exclusions`](crate::Topology::exclusions).
    pub fn exclusions(&self, atom: usize) -> &[u32] {
        &self.exclusions[atom]
    }

    /// The placement recipes of step `k`'s sites, in site order (lattice
    /// decoration reads refs/offsets to solve torsion variables exactly).
    pub(crate) fn step_sites(&self, k: usize) -> &[Site] {
        let (a, b) = self.steps[k];
        &self.sites[a..b]
    }

    /// Template atom indices placed by step `k`.
    pub fn step_atoms(&self, k: usize) -> impl Iterator<Item = usize> + '_ {
        let (a, b) = self.steps[k];
        self.sites[a..b].iter().map(|s| s.atom)
    }

    /// Number of atoms step `k` places.
    pub fn step_len(&self, k: usize) -> usize {
        let (a, b) = self.steps[k];
        b - a
    }

    /// Place the rigid seed at `origin` under rotation `rot` (row-major).
    pub fn place_seed(&self, origin: [F; 3], rot: &[[F; 3]; 3], coords: &mut [[F; 3]]) {
        for (k, &a) in self.seed.iter().enumerate() {
            let v = self.seed_local[k];
            coords[a] = [
                origin[0] + rot[0][0] * v[0] + rot[0][1] * v[1] + rot[0][2] * v[2],
                origin[1] + rot[1][0] * v[0] + rot[1][1] * v[1] + rot[1][2] * v[2],
                origin[2] + rot[2][0] * v[0] + rot[2][1] * v[1] + rot[2][2] * v[2],
            ];
        }
    }

    /// Place every atom of step `k`. `vars` must already hold the value of the
    /// variable this step owns (and of every earlier variable).
    pub fn place_step(&self, k: usize, vars: &[F], coords: &mut [[F; 3]]) {
        let (a, b) = self.steps[k];
        for s in &self.sites[a..b] {
            let phi = match s.follows {
                Some((v, offset)) => vars[v] + offset,
                None => s.torsion,
            };
            coords[s.atom] = nerf(
                coords[s.refs[2]],
                coords[s.refs[1]],
                coords[s.refs[0]],
                s.bond,
                s.angle,
                phi,
            );
        }
    }

    /// [`place_step`][Self::place_step] with per-site placement angles
    /// overriding the template's (`angles[i]` pairs with the step's `i`-th
    /// site). The CG angle-prior path — an all-atom template keeps its stiff
    /// template angles via `place_step` instead.
    pub fn place_step_with_angles(
        &self,
        k: usize,
        vars: &[F],
        angles: &[F],
        coords: &mut [[F; 3]],
    ) {
        let (a, b) = self.steps[k];
        debug_assert_eq!(angles.len(), b - a);
        for (s, &angle) in self.sites[a..b].iter().zip(angles) {
            let phi = match s.follows {
                Some((v, offset)) => vars[v] + offset,
                None => s.torsion,
            };
            coords[s.atom] = nerf(
                coords[s.refs[2]],
                coords[s.refs[1]],
                coords[s.refs[0]],
                s.bond,
                angle,
                phi,
            );
        }
    }

    /// Template placement angles of step `k`'s sites, in site order.
    pub fn step_angles(&self, k: usize) -> impl Iterator<Item = F> + '_ {
        let (a, b) = self.steps[k];
        self.sites[a..b].iter().map(|s| s.angle)
    }

    /// The template's own value for the variable step `k` owns.
    ///
    /// # Panics
    ///
    /// Panics when step `k` owns no variable ([`step_var`][Self::step_var]
    /// returns `None` — the seed-adjacent rigid step, or any step of a
    /// molecule without rotatable bonds). Callers gate on `step_var` first.
    pub fn template_var(&self, k: usize) -> F {
        assert!(
            self.step_var[k].is_some(),
            "template_var on step {k}, which owns no free variable"
        );
        let (a, _) = self.steps[k];
        self.sites[a].torsion
    }
}

// ── graph helpers ──────────────────────────────────────────────────────────

/// Copy of `graph` with every unclassed bond re-classed as `Single`, then the
/// perceived rotatable bonds as an unordered index-pair set.
///
/// Formats that carry connectivity without orders — PDB `CONECT`, GROMACS
/// `.top`, XYZ `Connct` — read back [`BondType::Unknown`], because molrs
/// reports what the file said rather than guessing. Perception only accepts
/// [`BondType::Single`], so without this fallback a PDB template would yield
/// zero rotatable bonds and grow as a rigid body.
pub(crate) fn rotatable_bond_keys(graph: &Atomistic) -> std::collections::HashSet<(usize, usize)> {
    use molrs::perceive::rotatable::detect_rotatable_bonds_with_downstream;

    let mut perceived = graph.clone();
    let unclassed: Vec<_> = perceived
        .bonds()
        .filter(|(id, _)| perceived.bond_type(*id) == BondType::Unknown)
        .map(|(id, _)| id)
        .collect();
    for id in unclassed {
        let _ = perceived.set_bond_type(id, BondType::Single);
    }
    // `RotatableBond.j` / `.k` are positional indices in `Atomistic::atoms`
    // order, which `Atomistic::from_frame` builds in frame row order — the
    // same index space as the template coordinates.
    detect_rotatable_bonds_with_downstream(&perceived)
        .iter()
        .map(|b| (b.j.min(b.k), b.j.max(b.k)))
        .collect()
}

/// One endpoint of a longest shortest-path in the graph — a chain end for a
/// linear polymer, so growth runs along the backbone instead of starting in
/// the middle and having to grow two ways at once.
fn diameter_endpoint(topo: &Topology) -> usize {
    let far = |from: usize| -> usize {
        let n = topo.natoms();
        let mut dist = vec![usize::MAX; n];
        let mut queue = std::collections::VecDeque::new();
        dist[from] = 0;
        queue.push_back(from);
        let mut best = from;
        while let Some(a) = queue.pop_front() {
            if dist[a] > dist[best] {
                best = a;
            }
            for &b in topo.neighbors(a) {
                let b = b as usize;
                if dist[b] == usize::MAX {
                    dist[b] = dist[a] + 1;
                    queue.push_back(b);
                }
            }
        }
        best
    };
    far(far(0))
}

fn bfs_order(
    topo: &Topology,
    root: usize,
    n: usize,
) -> Result<(Vec<usize>, Vec<Option<usize>>), GrowError> {
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut seen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut queue = std::collections::VecDeque::new();
    seen[root] = true;
    queue.push_back(root);
    while let Some(a) = queue.pop_front() {
        order.push(a);
        for &b in topo.neighbors(a) {
            let b = b as usize;
            if !seen[b] {
                seen[b] = true;
                parent[b] = Some(a);
                queue.push_back(b);
            }
        }
    }
    if order.len() != n {
        return Err(GrowError::Topology(TopologyError::Disconnected));
    }
    Ok((order, parent))
}

/// A placed neighbour of `at`, preferring its BFS parent, never `exclude`.
fn pick_ref(
    topo: &Topology,
    placed: &[bool],
    at: usize,
    exclude: usize,
    prefer: Option<usize>,
) -> Option<usize> {
    if let Some(p) = prefer
        && p != exclude
        && placed[p]
    {
        return Some(p);
    }
    topo.neighbors(at)
        .iter()
        .map(|&c| c as usize)
        .find(|&c| c != exclude && placed[c])
}

// ── geometry ───────────────────────────────────────────────────────────────

#[inline]
fn sub(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot(a: [F; 3], b: [F; 3]) -> F {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn norm(a: [F; 3]) -> F {
    dot(a, a).sqrt()
}
#[inline]
fn unit(a: [F; 3]) -> [F; 3] {
    let n = norm(a).max(crate::numerics::near_zero_norm_floor());
    [a[0] / n, a[1] / n, a[2] / n]
}

fn angle_at(a: [F; 3], vertex: [F; 3], c: [F; 3]) -> F {
    let u = unit(sub(a, vertex));
    let v = unit(sub(c, vertex));
    dot(u, v).clamp(-1.0, 1.0).acos()
}

/// Dihedral `a-b-c-d`, in `(-π, π]`.
pub(crate) fn dihedral(a: [F; 3], b: [F; 3], c: [F; 3], d: [F; 3]) -> F {
    let b1 = sub(b, a);
    let b2 = sub(c, b);
    let b3 = sub(d, c);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m1 = cross(n1, unit(b2));
    let x = dot(n1, n2);
    let y = dot(m1, n2);
    let d = (-y).atan2(x);
    // atan2 can return exactly -π on planar input with a signed zero;
    // fold it onto +π so the documented (-π, π] range holds everywhere.
    if d == -(std::f64::consts::PI as F) {
        std::f64::consts::PI as F
    } else {
        d
    }
}

/// Natural-extension reference frame: place `d` given `a-b-c` and `(bond,
/// angle, torsion)`, the exact inverse of [`dihedral`] / [`angle_at`].
pub(crate) fn nerf(a: [F; 3], b: [F; 3], c: [F; 3], bond: F, angle: F, torsion: F) -> [F; 3] {
    let bc = unit(sub(c, b));
    let n = unit(cross(sub(b, a), bc));
    let m = cross(n, bc);
    let (st, ct) = angle.sin_cos();
    let (sp, cp) = torsion.sin_cos();
    let d2 = [-bond * ct, bond * st * cp, bond * st * sp];
    [
        c[0] + bc[0] * d2[0] + m[0] * d2[1] + n[0] * d2[2],
        c[1] + bc[1] * d2[0] + m[1] * d2[1] + n[1] * d2[2],
        c[2] + bc[2] * d2[0] + m[2] * d2[1] + n[2] * d2[2],
    ]
}

#[inline]
fn wrap_pi(x: F) -> F {
    use std::f64::consts::PI;
    let two_pi = 2.0 * PI as F;
    let mut v = x % two_pi;
    if v > PI as F {
        v -= two_pi;
    } else if v <= -(PI as F) {
        v += two_pi;
    }
    v
}
