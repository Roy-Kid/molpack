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

use molrs::op::types::F;
use molrs::store::Frame;
use molrs::system::Atomistic;

use crate::grow::GrowError;

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

impl InternalTree {
    /// Decompose a template frame using `weights` as the intramolecular skip
    /// table. A non-0/1 weight is [`GrowError::NonBinarySpecialBond`] before
    /// exclusions are compiled. The all-atom convention is
    /// [`molrs::system::BondDistanceWeights::from_exclusion_depth`]`(3)`.
    ///
    /// Coordinates are in Å. The frame must carry an `atoms` block with
    /// `x`/`y`/`z` and a connected acyclic bond graph of at least 3 atoms.
    pub fn from_frame(
        frame: &Frame,
        weights: &molrs::system::BondDistanceWeights,
    ) -> Result<Self, GrowError> {
        if let Some((index, weight)) = super::config::binary_violation(weights) {
            return Err(GrowError::NonBinarySpecialBond { index, weight });
        }
        let (topo, xyz) = super::topology_for_growth(frame)?;
        let n = topo.n_atoms();

        let root = diameter_endpoint(&topo);
        let (order, parent) = bfs_order(&topo, root, n);

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
            let angle = angle(xyz[g], xyz[p], xyz[i]);
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

        let exclusions = topo
            .exclusions(weights)
            .into_iter()
            .map(|row| row.into_iter().map(|i| i as u32).collect())
            .collect();

        Ok(Self {
            seed,
            seed_local,
            sites,
            steps,
            step_var,
            n_vars,
            exclusions,
        })
    }

    /// Atoms in the tree: the three-atom seed plus one site per other atom.
    #[cfg(test)]
    pub fn n_atoms(&self) -> usize {
        3 + self.sites.len()
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
    /// position. Produced by [`molrs::system::Topology::exclusions`].
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

/// The rotatable bonds of `graph` (see [`crate::template::rotatable_bonds`])
/// as an unordered index-pair set.
pub(crate) fn rotatable_bond_keys(graph: &Atomistic) -> std::collections::HashSet<(usize, usize)> {
    // `RotatableBond.j` / `.k` are positional indices in `Atomistic::atoms`
    // order, which `Atomistic::from_frame` builds in frame row order — the
    // same index space as the template coordinates.
    crate::template::rotatable_bonds(graph)
        .iter()
        .map(|b| (b.j.min(b.k), b.j.max(b.k)))
        .collect()
}

/// One endpoint of a longest shortest-path in the graph — a chain end for a
/// linear polymer, so growth runs along the backbone instead of starting in
/// the middle and having to grow two ways at once.
fn diameter_endpoint(topo: &molrs::system::Topology) -> usize {
    // The farthest atom from `from`; ties go to the lowest index.
    let far = |from: usize| -> usize {
        let dist = topo.distances(from);
        (0..dist.len()).fold(from, |best, a| if dist[a] > dist[best] { a } else { best })
    };
    far(far(0))
}

fn bfs_order(
    topo: &molrs::system::Topology,
    root: usize,
    n: usize,
) -> (Vec<usize>, Vec<Option<usize>>) {
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut seen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut queue = std::collections::VecDeque::new();
    seen[root] = true;
    queue.push_back(root);
    while let Some(a) = queue.pop_front() {
        order.push(a);
        for b in topo.neighbors(a) {
            if !seen[b] {
                seen[b] = true;
                parent[b] = Some(a);
                queue.push_back(b);
            }
        }
    }
    debug_assert_eq!(order.len(), n);
    (order, parent)
}

/// A placed neighbour of `at`, preferring its BFS parent, never `exclude`.
fn pick_ref(
    topo: &molrs::system::Topology,
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
        .into_iter()
        .find(|&c| c != exclude && placed[c])
}

// ── geometry ───────────────────────────────────────────────────────────────

// Internal coordinates are molrs's: `op::vec3::{angle, dihedral}` read them
// off a template and `op::rigid::nerf` is their exact inverse.
use molrs::op::rigid::nerf;
use molrs::op::vec3::{angle, cross, dihedral, norm, sub};

/// Wrap an angle (radians) into `(-π, π]`.
#[inline]
pub(crate) fn wrap_pi(x: F) -> F {
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
