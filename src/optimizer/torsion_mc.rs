//! Torsion Monte-Carlo optimizer implementing [`molrs::optimize::Optimizer`].
//!
//! Lives in molpack (not molrs): uses packer-local geometry helpers and
//! self-avoidance scoring on a Frame assembled by the packer.

use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

use molrs::core::Frame;
use molrs::core::{Atomistic, NodeId};
use molrs::core::{BondDistanceWeights, Topology};
use molrs::op::F;
use molrs::op::centroid;
use molrs::op::{axis_angle, rotation_about, transform_point, vec3};
use molrs::optimize::{OptimizationReport, Optimizer};
use molrs::perceive::RotatableBond;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::random::uniform01_core;

/// Monte-Carlo torsion-angle optimizer for flexible molecules.
///
/// Largest torsion change one proposal makes (radians).
const MAX_DELTA: F = (PI / 6.0) as F;

/// Implements [`Optimizer`]: each `minimize` proposes rotations about rotatable
/// bonds on the Frame's free atoms and accepts against self-avoidance energy
/// (plus optional soft contact with fixed environment atoms present in the
/// Frame). Packing non-harm is enforced by the packer after write-back.
#[derive(Debug, Clone)]
pub struct TorsionMcOptimizer {
    bonds: Vec<RotatableBond>,
    steps: usize,
    temperature: F,
    self_avoidance_radius: F,
    /// The molecule's bond graph, kept so the exclusions can be re-derived
    /// when [`with_special_bonds`](Self::with_special_bonds) changes the table.
    topology: Topology,
    excluded_pairs: HashSet<(usize, usize)>,
    seed: u64,
}

impl TorsionMcOptimizer {
    /// Build from a molecule's topology, perceiving its rotatable bonds.
    ///
    /// A bond whose class the input left unstated counts as a rotatable
    /// single bond — molpack's one policy for every solver. Self-avoidance
    /// skips pairs at bond distance 1–3 (the all-atom default); a
    /// coarse-grained model states its own table with
    /// [`with_special_bonds`](Self::with_special_bonds).
    pub fn new(graph: &Atomistic) -> Self {
        // Atoms by their position in `Atomistic::atoms`, the order the
        // rotatable bonds' indices use.
        let id_to_idx: HashMap<NodeId, usize> = graph
            .atoms()
            .enumerate()
            .map(|(idx, (id, _))| (id, idx))
            .collect();
        let edges: Vec<[usize; 2]> = graph
            .bonds()
            .map(|(_, b)| [id_to_idx[&b.nodes[0]], id_to_idx[&b.nodes[1]]])
            .collect();
        let topology = Topology::from_edges(id_to_idx.len(), &edges);
        Self {
            bonds: crate::template::rotatable_bonds(graph),
            steps: 10,
            temperature: 1.0,
            self_avoidance_radius: 0.0,
            excluded_pairs: excluded_pairs(
                &topology,
                &BondDistanceWeights::from_exclusion_depth(3),
            ),
            topology,
            seed: 1,
        }
    }

    /// Which intramolecular pairs self-avoidance skips: every pair whose
    /// bond-distance weight is `0` (the same table
    /// [`Target::with_special_bonds`](crate::Target) takes).
    pub fn with_special_bonds(mut self, weights: BondDistanceWeights) -> Self {
        self.excluded_pairs = excluded_pairs(&self.topology, &weights);
        self
    }

    /// Number of rotatable bonds perceived for this molecule, after the
    /// unclassed-bond fallback described on [`new`](Self::new).
    ///
    /// Zero means every proposal this optimizer makes is the identity — worth
    /// asserting on when wiring one up, since a no-op optimizer is otherwise
    /// indistinguishable from one that simply never improves anything.
    pub fn rotatable_bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn with_temperature(mut self, t: F) -> Self {
        self.temperature = t;
        self
    }
    pub fn with_steps(mut self, n: usize) -> Self {
        self.steps = n;
        self
    }
    pub fn with_self_avoidance(mut self, radius: F) -> Self {
        self.self_avoidance_radius = radius;
        self
    }
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }
}

impl Optimizer for TorsionMcOptimizer {
    fn minimize(&mut self, frame: &mut Frame) -> Result<OptimizationReport, String> {
        if self.bonds.is_empty() {
            return Ok(OptimizationReport {
                converged: true,
                n_steps: 0,
                final_energy: 0.0,
                final_fmax: 0.0,
                final_grad_rms: 0.0,
            });
        }
        let xyz = frame.coords().map_err(|e| e.to_string())?;
        let n = xyz.nrows();
        let coords = crate::template::coord_rows(&xyz);

        // Free mask: only free atoms may be torsion-rotated (environment fixed).
        let free: Vec<bool> = match frame
            .get("atoms")
            .and_then(|a| a.get("free"))
            .and_then(molrs::core::Column::as_bool)
        {
            Some(col) if col.len() == n => col.iter().copied().collect(),
            _ => vec![true; n],
        };

        let mut rng = SmallRng::seed_from_u64(self.seed);
        self.seed = self.seed.wrapping_add(1);

        let use_sa = self.self_avoidance_radius > 0.0;
        let mut best = coords.clone();
        let mut best_e = energy(
            &best,
            use_sa,
            self.self_avoidance_radius,
            &self.excluded_pairs,
        );
        let mut trial = best.clone();
        let mut accepts = 0usize;

        for _ in 0..self.steps {
            if self.bonds.is_empty() {
                break;
            }
            let bond_idx = (rng.next_u32() as usize) % self.bonds.len();
            let bond = &self.bonds[bond_idx];
            // Skip moves that would rotate only fixed atoms.
            if !bond
                .downstream
                .iter()
                .any(|&i| free.get(i).copied().unwrap_or(false))
                && !free.get(bond.j).copied().unwrap_or(false)
            {
                continue;
            }
            let delta = (uniform01_core(&mut rng) * 2.0 - 1.0) * MAX_DELTA;
            trial.copy_from_slice(&best);
            rotate_around_bond(&mut trial, bond, delta);
            // Restore fixed atoms.
            for i in 0..n {
                if !free[i] {
                    trial[i] = best[i];
                }
            }
            recenter_free(&mut trial, &free);
            let e = energy(
                &trial,
                use_sa,
                self.self_avoidance_radius,
                &self.excluded_pairs,
            );
            if metropolis_accept(e, best_e, self.temperature, &mut rng) {
                std::mem::swap(&mut best, &mut trial);
                best_e = e;
                accepts += 1;
            }
        }

        let out = ndarray::Array2::from(best);
        frame.set_coords(out.view()).map_err(|e| e.to_string())?;
        Ok(OptimizationReport {
            converged: accepts > 0 || self.steps == 0,
            n_steps: self.steps,
            final_energy: best_e,
            final_fmax: 0.0,
            final_grad_rms: 0.0,
        })
    }
}

fn energy(coords: &[[F; 3]], use_sa: bool, radius: F, excluded: &HashSet<(usize, usize)>) -> F {
    if !use_sa {
        return 0.0;
    }
    self_avoidance_penalty(coords, radius, excluded)
}

fn self_avoidance_penalty(coords: &[[F; 3]], radius: F, excluded: &HashSet<(usize, usize)>) -> F {
    let cutoff = 2.0 * radius;
    let cutoff_sq = cutoff * cutoff;
    let n = coords.len();
    let mut penalty: F = 0.0;
    for i in 0..n {
        let ci = coords[i];
        for (j, cj) in coords.iter().enumerate().skip(i + 1) {
            let dx = ci[0] - cj[0];
            let dy = ci[1] - cj[1];
            let dz = ci[2] - cj[2];
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq < cutoff_sq && !excluded.contains(&(i, j)) {
                let gap = dist_sq - cutoff_sq;
                penalty += gap * gap;
            }
        }
    }
    penalty
}

/// Rotate the downstream side of `bond` by `angle` radians about the `j → k`
/// axis. A degenerate bond (ends closer than molrs's
/// `MIN_DIRECTION_LENGTH`, so not a direction) is left alone.
fn rotate_around_bond(coords: &mut [[F; 3]], bond: &RotatableBond, angle: F) {
    let (j, k) = (coords[bond.j], coords[bond.k]);
    let axis = vec3::sub(k, j);
    if vec3::normalize(axis).is_none() {
        return;
    }
    let Some(rotation) = axis_angle(axis, angle) else {
        return;
    };
    let motion = rotation_about(rotation, j);
    for &idx in &bond.downstream {
        coords[idx] = transform_point(&motion, coords[idx]);
    }
}

/// Shift the free atoms so their centroid (molrs's
/// [`centroid`](molrs::op::centroid), unit weights) sits at the
/// origin. No free atom: nothing moves.
fn recenter_free(coords: &mut [[F; 3]], free: &[bool]) {
    let points: Vec<[F; 3]> = coords
        .iter()
        .zip(free)
        .filter_map(|(p, &f)| f.then_some(*p))
        .collect();
    let Some(c) = centroid(&points, &vec![1.0; points.len()]) else {
        return;
    };
    for (i, p) in coords.iter_mut().enumerate() {
        if free[i] {
            p[0] -= c[0];
            p[1] -= c[1];
            p[2] -= c[2];
        }
    }
}

fn metropolis_accept(f_trial: F, f_current: F, temperature: F, rng: &mut dyn Rng) -> bool {
    if f_trial <= f_current {
        return true;
    }
    if temperature <= 0.0 {
        return false;
    }
    let delta = (f_trial - f_current) / temperature;
    uniform01_core(rng) < (-delta).exp()
}

/// Unordered `(i, j)`, `i < j`, pairs the table exempts, from molrs's
/// bond-distance walk.
fn excluded_pairs(topology: &Topology, weights: &BondDistanceWeights) -> HashSet<(usize, usize)> {
    topology
        .exclusions(weights)
        .into_iter()
        .enumerate()
        .flat_map(|(i, partners)| {
            partners
                .into_iter()
                .filter(move |&j| j > i)
                .map(move |j| (i, j))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::chain_graph;

    /// The all-atom default skips 1-2, 1-3 and 1-4 pairs and scores 1-5.
    #[test]
    fn default_exclusions_reach_bond_distance_three() {
        let opt = TorsionMcOptimizer::new(&chain_graph(5));
        assert!(opt.excluded_pairs.contains(&(0, 3)));
        assert!(!opt.excluded_pairs.contains(&(0, 4)));
    }

    /// A coarse-grained table is per-target data: excluding only bonded
    /// neighbours leaves 1-3 pairs to self-avoidance.
    #[test]
    fn special_bonds_set_the_exclusion_depth() {
        let opt = TorsionMcOptimizer::new(&chain_graph(5))
            .with_special_bonds(BondDistanceWeights::from_exclusion_depth(1));
        assert!(opt.excluded_pairs.contains(&(0, 1)));
        assert!(!opt.excluded_pairs.contains(&(0, 2)));
    }
}
