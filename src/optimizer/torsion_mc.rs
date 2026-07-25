//! Torsion Monte-Carlo optimizer implementing [`molrs::optimize::Optimizer`].
//!
//! Lives in molpack (not molrs): uses packer-local geometry helpers and
//! self-avoidance scoring on a Frame assembled by the packer.

#![cfg(feature = "ff")]

use std::collections::HashSet;
use std::f64::consts::PI;

use molrs::ff::potential::{extract_coords, write_coords};
use molrs::optimize::{OptReport, Optimizer};
use molrs::perceive::rotatable::{
    RotatableBond, atom_id_to_index, detect_rotatable_bonds_with_downstream,
};
use molrs::store::frame::Frame;
use molrs::system::atomistic::Atomistic;
use molrs::types::F;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::numerics::near_zero_norm_floor;
use crate::random::uniform01_core;

/// Monte-Carlo torsion-angle optimizer for flexible molecules.
///
/// Implements [`Optimizer`]: each `run` proposes rotations about rotatable
/// bonds on the Frame's free atoms and accepts against self-avoidance energy
/// (plus optional soft contact with fixed environment atoms present in the
/// Frame). Packing non-harm is enforced by the packer after write-back.
#[derive(Debug, Clone)]
pub struct TorsionMcOptimizer {
    bonds: Vec<RotatableBond>,
    max_delta: F,
    steps: usize,
    temperature: F,
    self_avoidance_radius: F,
    excluded_pairs: HashSet<(usize, usize)>,
    seed: u64,
}

impl TorsionMcOptimizer {
    pub fn new(graph: &Atomistic) -> Self {
        let bonds = detect_rotatable_bonds_with_downstream(graph);
        let excluded_pairs = compute_excluded_pairs(graph);
        Self {
            bonds,
            max_delta: (PI / 6.0) as F,
            steps: 10,
            temperature: 1.0,
            self_avoidance_radius: 0.0,
            excluded_pairs,
            seed: 1,
        }
    }

    pub fn with_temperature(mut self, t: F) -> Self {
        self.temperature = t;
        self
    }
    pub fn with_steps(mut self, n: usize) -> Self {
        self.steps = n;
        self
    }
    pub fn with_max_delta(mut self, rad: F) -> Self {
        self.max_delta = rad;
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
    fn run(&mut self, frame: &mut Frame) -> Result<OptReport, String> {
        if self.bonds.is_empty() {
            return Ok(OptReport {
                converged: true,
                n_steps: 0,
                final_energy: 0.0,
                final_fmax: 0.0,
            });
        }
        let flat = extract_coords(frame)?;
        let n = flat.len() / 3;
        let mut coords: Vec<[F; 3]> = (0..n)
            .map(|i| [flat[3 * i], flat[3 * i + 1], flat[3 * i + 2]])
            .collect();

        // Free mask: only free atoms may be torsion-rotated (environment fixed).
        let free: Vec<bool> = match frame.get("atoms").and_then(|a| a.get_bool("free")) {
            Some(col) if col.len() == n => col.iter().copied().collect(),
            _ => vec![true; n],
        };

        let mut rng = SmallRng::seed_from_u64(self.seed);
        self.seed = self.seed.wrapping_add(1);

        let use_sa = self.self_avoidance_radius > 0.0;
        let mut best = coords.clone();
        let mut best_e = energy(&best, use_sa, self.self_avoidance_radius, &self.excluded_pairs);
        let mut trial = best.clone();
        let mut accepts = 0usize;

        for _ in 0..self.steps {
            if self.bonds.is_empty() {
                break;
            }
            let bond_idx = (rng.next_u32() as usize) % self.bonds.len();
            let bond = &self.bonds[bond_idx];
            // Skip moves that would rotate only fixed atoms.
            if !bond.downstream.iter().any(|&i| free.get(i).copied().unwrap_or(false))
                && !free.get(bond.j).copied().unwrap_or(false)
            {
                continue;
            }
            let delta = (uniform01_core(&mut rng) * 2.0 - 1.0) * self.max_delta;
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

        let mut out = Vec::with_capacity(n * 3);
        for p in &best {
            out.extend_from_slice(p);
        }
        write_coords(frame, &out)?;
        Ok(OptReport {
            converged: accepts > 0 || self.steps == 0,
            n_steps: self.steps,
            final_energy: best_e,
            final_fmax: 0.0,
        })
    }
}

fn energy(
    coords: &[[F; 3]],
    use_sa: bool,
    radius: F,
    excluded: &HashSet<(usize, usize)>,
) -> F {
    if !use_sa {
        return 0.0;
    }
    self_avoidance_penalty(coords, radius, excluded)
}

fn self_avoidance_penalty(
    coords: &[[F; 3]],
    radius: F,
    excluded: &HashSet<(usize, usize)>,
) -> F {
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

fn rotate_around_bond(coords: &mut [[F; 3]], bond: &RotatableBond, angle: F) {
    let j = coords[bond.j];
    let k = coords[bond.k];
    let mut u = [k[0] - j[0], k[1] - j[1], k[2] - j[2]];
    let norm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    if norm < near_zero_norm_floor() {
        return;
    }
    u[0] /= norm;
    u[1] /= norm;
    u[2] /= norm;
    let origin = j;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    for &idx in &bond.downstream {
        let p0 = coords[idx];
        let p = [p0[0] - origin[0], p0[1] - origin[1], p0[2] - origin[2]];
        let udotp = u[0] * p[0] + u[1] * p[1] + u[2] * p[2];
        let cross = [
            u[1] * p[2] - u[2] * p[1],
            u[2] * p[0] - u[0] * p[2],
            u[0] * p[1] - u[1] * p[0],
        ];
        coords[idx] = [
            p[0] * cos_a + cross[0] * sin_a + u[0] * udotp * (1.0 - cos_a) + origin[0],
            p[1] * cos_a + cross[1] * sin_a + u[1] * udotp * (1.0 - cos_a) + origin[1],
            p[2] * cos_a + cross[2] * sin_a + u[2] * udotp * (1.0 - cos_a) + origin[2],
        ];
    }
}

fn recenter_free(coords: &mut [[F; 3]], free: &[bool]) {
    let mut n = 0.0 as F;
    let mut c = [0.0 as F; 3];
    for (i, p) in coords.iter().enumerate() {
        if free[i] {
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
            n += 1.0;
        }
    }
    if n < 1.0 {
        return;
    }
    c[0] /= n;
    c[1] /= n;
    c[2] /= n;
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

fn compute_excluded_pairs(graph: &Atomistic) -> HashSet<(usize, usize)> {
    let id_to_idx = atom_id_to_index(graph);
    let atom_ids: Vec<_> = graph.atoms().map(|(id, _)| id).collect();
    let mut adj: std::collections::HashMap<_, Vec<_>> = std::collections::HashMap::new();
    for &id in &atom_ids {
        adj.insert(id, graph.neighbors(id).collect());
    }
    let mut excluded = HashSet::new();
    for &root in &atom_ids {
        let root_idx = id_to_idx[&root];
        for &n1 in adj.get(&root).unwrap_or(&Vec::new()) {
            let n1_idx = id_to_idx[&n1];
            excluded.insert((root_idx.min(n1_idx), root_idx.max(n1_idx)));
            for &n2 in adj.get(&n1).unwrap_or(&Vec::new()) {
                if n2 == root {
                    continue;
                }
                let n2_idx = id_to_idx[&n2];
                excluded.insert((root_idx.min(n2_idx), root_idx.max(n2_idx)));
                for &n3 in adj.get(&n2).unwrap_or(&Vec::new()) {
                    if n3 == root || n3 == n1 {
                        continue;
                    }
                    let n3_idx = id_to_idx[&n3];
                    excluded.insert((root_idx.min(n3_idx), root_idx.max(n3_idx)));
                }
            }
        }
    }
    excluded
}
