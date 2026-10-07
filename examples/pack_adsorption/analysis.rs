//! Train / loop / tail analysis of an adsorbed chain layer.
//!
//! The classical Scheutjens–Fleer decomposition: a monomer within the surface
//! layer is *adsorbed*; a maximal run of adsorbed monomers is a **train**, a
//! maximal run of desorbed monomers between two trains is a **loop**, and a
//! desorbed run at either chain end is a **tail**.
//!
//! This is the observable the example exists to produce. An extended chain
//! cannot show it — only conformational relaxation during packing can put
//! every k-th bead on the surface while the segments between them arch away.

use molrs::op::F;

/// Train / loop / tail counts and lengths over one chain layer.
#[derive(Debug, Default)]
pub struct LayerStats {
    pub trains: Vec<usize>,
    pub loops: Vec<usize>,
    pub tails: Vec<usize>,
    /// Sticky beads inside the adsorption slab, over all sticky beads.
    pub sticky_adsorbed: usize,
    pub sticky_total: usize,
}

impl LayerStats {
    pub fn sticky_fraction(&self) -> F {
        if self.sticky_total == 0 {
            0.0
        } else {
            self.sticky_adsorbed as F / self.sticky_total as F
        }
    }

    pub fn mean(v: &[usize]) -> F {
        if v.is_empty() {
            0.0
        } else {
            v.iter().sum::<usize>() as F / v.len() as F
        }
    }
}

/// Split one chain's adsorbed/desorbed mask into trains, loops and tails.
///
/// Runs are classified by position: a desorbed run touching either chain end is
/// a tail, one bounded by trains on both sides is a loop.
fn classify(mask: &[bool], stats: &mut LayerStats) {
    let n = mask.len();
    let mut i = 0;
    while i < n {
        let value = mask[i];
        let start = i;
        while i < n && mask[i] == value {
            i += 1;
        }
        let len = i - start;
        if value {
            stats.trains.push(len);
        } else if start == 0 || i == n {
            stats.tails.push(len);
        } else {
            stats.loops.push(len);
        }
    }
}

/// Analyse `n_chains` chains of `n_beads` laid out copy-by-copy, bead-by-bead
/// from `chain_positions[0]`.
///
/// `adsorbed_below` is the z below which a bead counts as being in the surface
/// layer; `sticky` holds the 0-based sticky indices within one chain.
pub fn layer_stats(
    chain_positions: &[[F; 3]],
    n_chains: usize,
    n_beads: usize,
    sticky: &[usize],
    adsorbed_below: F,
) -> LayerStats {
    let mut stats = LayerStats::default();
    for c in 0..n_chains {
        let base = c * n_beads;
        let mask: Vec<bool> = (0..n_beads)
            .map(|i| chain_positions[base + i][2] < adsorbed_below)
            .collect();
        classify(&mask, &mut stats);
        for &s in sticky {
            stats.sticky_total += 1;
            if mask[s] {
                stats.sticky_adsorbed += 1;
            }
        }
    }
    stats
}

/// A z-density histogram, as (bin centre, count) pairs.
pub fn z_profile(positions: &[[F; 3]], z_min: F, z_max: F, nbins: usize) -> Vec<(F, usize)> {
    let width = (z_max - z_min) / nbins as F;
    let mut counts = vec![0usize; nbins];
    for p in positions {
        if p[2] < z_min || p[2] >= z_max {
            continue;
        }
        let b = ((p[2] - z_min) / width) as usize;
        if b < nbins {
            counts[b] += 1;
        }
    }
    counts
        .into_iter()
        .enumerate()
        .map(|(b, c)| (z_min + (b as F + 0.5) * width, c))
        .collect()
}

/// Render a histogram as a fixed-width bar chart.
pub fn render_profile(profile: &[(F, usize)], width: usize) -> String {
    let peak = profile.iter().map(|(_, c)| *c).max().unwrap_or(0).max(1);
    profile
        .iter()
        .map(|(z, c)| {
            let bar = c * width / peak;
            format!("  z={z:6.1}  {:<width$} {c}\n", "█".repeat(bar))
        })
        .collect()
}
