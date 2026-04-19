//! Reusable temporary buffers for objective/gradient and movebad paths.

use molrs::types::F;
#[cfg(feature = "rayon")]
use rayon::prelude::*;

/// Crossover where the parallel gradient-merge starts paying off.
///
/// The merge sums `n_threads` partial buffers of `ntotat` triples into
/// `gxcar`, i.e. roughly `n_threads × ntotat × 3` scalar adds. Rayon's
/// `par_chunks_mut` task-dispatch overhead is ~1 µs per chunk on
/// current hardware; at `n = 4096` triples the serial merge itself is
/// on the order of 10 µs, so any smaller `ntotat` goes the serial
/// route. See the `partial_gradient_merge` bench for the data behind
/// this threshold.
#[cfg(feature = "rayon")]
pub const MERGE_PARALLEL_THRESHOLD_NTOTAT: usize = 4096;

/// Rayon chunk size (in atoms) for the parallel merge. ~256-atom
/// chunks give each thread 256 × 3 = 768 scalar adds per partial,
/// which is long enough to amortize the per-task overhead but short
/// enough that all cores get roughly equal work on typical
/// `n_threads × ntotat` ratios.
#[cfg(feature = "rayon")]
pub const MERGE_PARALLEL_CHUNK_ATOMS: usize = 256;
/// Reusable mutable buffers shared across packing iterations.
pub struct WorkBuffers {
    /// Cartesian gradient accumulator used by objective gradient evaluation.
    pub gxcar: Vec<[F; 3]>,
    /// Per-thread partial gradient buffers for the rayon-parallel
    /// accumulators. Allocated lazily on first use; each buffer is sized
    /// to `ntotat`. Kept across iterations to amortize allocation cost,
    /// which dominates the parallel pair kernel at large `ntotat`.
    #[cfg(feature = "rayon")]
    pub partial_gxcar: Vec<Vec<[F; 3]>>,
    /// Temporary radius backup used by movebad/radius scaling paths.
    pub radiuswork: Vec<F>,
    /// Per-molecule score buffer used by flashsort/movebad ranking.
    pub fmol: Vec<F>,
    /// Index permutation buffer reused by flashsort in movebad.
    pub flash_ind: Vec<usize>,
    /// Histogram bucket buffer reused by flashsort.
    pub flash_l: Vec<usize>,
    /// Last x-vector whose expanded Cartesian geometry is still resident in `PackContext`.
    pub cached_x: Vec<F>,
    /// Active-type mask associated with `cached_x`.
    pub cached_comptype: Vec<bool>,
    /// Whether the cached geometry was built in init1 mode.
    pub cached_init1: bool,
    /// Cell grid signature for the cached geometry.
    pub cached_ncells: [usize; 3],
    pub cached_cell_length: [F; 3],
    pub cached_pbc_min: [F; 3],
    pub cached_pbc_length: [F; 3],
    pub cached_pbc_periodic: [bool; 3],
    /// Whether the cached geometry metadata is valid.
    pub cached_geometry_valid: bool,
}

impl WorkBuffers {
    pub fn new(ntotat: usize) -> Self {
        Self {
            gxcar: vec![[0.0; 3]; ntotat],
            #[cfg(feature = "rayon")]
            partial_gxcar: Vec::new(),
            radiuswork: vec![0.0; ntotat],
            fmol: Vec::new(),
            flash_ind: Vec::new(),
            flash_l: Vec::new(),
            cached_x: Vec::new(),
            cached_comptype: Vec::new(),
            cached_init1: false,
            cached_ncells: [0; 3],
            cached_cell_length: [0.0; 3],
            cached_pbc_min: [0.0; 3],
            cached_pbc_length: [0.0; 3],
            cached_pbc_periodic: [false; 3],
            cached_geometry_valid: false,
        }
    }

    /// Resize the per-thread partial gradient buffers to `n_threads`
    /// slots of `ntotat` entries each, filling them with zeros so callers
    /// can accumulate into them directly. Reuses existing allocations
    /// when sizes match — the parallel pair kernel is called repeatedly.
    #[cfg(feature = "rayon")]
    pub fn reset_partial_gxcar(&mut self, n_threads: usize, ntotat: usize) {
        if self.partial_gxcar.len() != n_threads {
            self.partial_gxcar = (0..n_threads).map(|_| vec![[0.0; 3]; ntotat]).collect();
            return;
        }
        for buf in &mut self.partial_gxcar {
            if buf.len() != ntotat {
                buf.clear();
                buf.resize(ntotat, [0.0; 3]);
            } else {
                buf.fill([0.0; 3]);
            }
        }
    }

    /// Sum every `partials[t][i]` triple into `out[i]`. Serial version:
    /// each partial buffer walked in order, one scalar-add column at a
    /// time. Complexity: O(`n_threads` × `ntotat`).
    ///
    /// Used for small `ntotat` where the rayon dispatch overhead of
    /// [`Self::merge_partials_parallel`] would dominate.
    pub fn merge_partials_serial(partials: &[Vec<[F; 3]>], out: &mut [[F; 3]]) {
        let ntotat = out.len();
        for partial in partials {
            debug_assert_eq!(partial.len(), ntotat);
            for i in 0..ntotat {
                out[i][0] += partial[i][0];
                out[i][1] += partial[i][1];
                out[i][2] += partial[i][2];
            }
        }
    }

    /// Parallel counterpart to [`Self::merge_partials_serial`]. Splits `out`
    /// into disjoint atom chunks of [`MERGE_PARALLEL_CHUNK_ATOMS`] and
    /// has each rayon task sum every partial buffer into its chunk.
    ///
    /// Gated at the call site on
    /// `ntotat >= MERGE_PARALLEL_THRESHOLD_NTOTAT`.
    #[cfg(feature = "rayon")]
    pub fn merge_partials_parallel(partials: &[Vec<[F; 3]>], out: &mut [[F; 3]]) {
        out.par_chunks_mut(MERGE_PARALLEL_CHUNK_ATOMS)
            .enumerate()
            .for_each(|(chunk_idx, out_chunk)| {
                let start = chunk_idx * MERGE_PARALLEL_CHUNK_ATOMS;
                for partial in partials {
                    // SAFETY-ish: the serial-side `reset_partial_gxcar`
                    // always sizes each partial to `ntotat`, and `out`
                    // is `sys.work.gxcar` which is also sized to
                    // `ntotat`. `start + out_chunk.len() <= ntotat`
                    // therefore holds by construction — the debug
                    // assertion pins the invariant.
                    debug_assert!(start + out_chunk.len() <= partial.len());
                    let src = &partial[start..start + out_chunk.len()];
                    for (dst, s) in out_chunk.iter_mut().zip(src) {
                        dst[0] += s[0];
                        dst[1] += s[1];
                        dst[2] += s[2];
                    }
                }
            });
    }

    pub fn ensure_atom_capacity(&mut self, ntotat: usize) {
        if self.gxcar.len() != ntotat {
            self.gxcar.resize(ntotat, [0.0; 3]);
        }
        if self.radiuswork.len() != ntotat {
            self.radiuswork.resize(ntotat, 0.0);
        }
        self.cached_geometry_valid = false;
    }

    #[allow(clippy::too_many_arguments)]
    pub fn matches_cached_geometry(
        &self,
        x: &[F],
        comptype: &[bool],
        init1: bool,
        ncells: [usize; 3],
        cell_length: [F; 3],
        pbc_min: [F; 3],
        pbc_length: [F; 3],
        pbc_periodic: [bool; 3],
    ) -> bool {
        self.cached_geometry_valid
            && self.cached_init1 == init1
            && self.cached_ncells == ncells
            && self.cached_cell_length == cell_length
            && self.cached_pbc_min == pbc_min
            && self.cached_pbc_length == pbc_length
            && self.cached_pbc_periodic == pbc_periodic
            && self.cached_x == x
            && self.cached_comptype == comptype
    }

    #[allow(clippy::too_many_arguments)]
    pub fn update_cached_geometry(
        &mut self,
        x: &[F],
        comptype: &[bool],
        init1: bool,
        ncells: [usize; 3],
        cell_length: [F; 3],
        pbc_min: [F; 3],
        pbc_length: [F; 3],
        pbc_periodic: [bool; 3],
    ) {
        self.cached_x.clear();
        self.cached_x.extend_from_slice(x);
        self.cached_comptype.clear();
        self.cached_comptype.extend_from_slice(comptype);
        self.cached_init1 = init1;
        self.cached_ncells = ncells;
        self.cached_cell_length = cell_length;
        self.cached_pbc_min = pbc_min;
        self.cached_pbc_length = pbc_length;
        self.cached_pbc_periodic = pbc_periodic;
        self.cached_geometry_valid = true;
    }
}

#[cfg(test)]
mod merge_tests {
    use super::*;

    fn make_partials(n_threads: usize, ntotat: usize) -> Vec<Vec<[F; 3]>> {
        (0..n_threads)
            .map(|t| {
                (0..ntotat)
                    .map(|i| {
                        let v = (t * 37 + i * 13) as F;
                        [v * 0.125, -v * 0.0625, v * 0.03125]
                    })
                    .collect()
            })
            .collect()
    }

    fn reference_sum(partials: &[Vec<[F; 3]>], ntotat: usize) -> Vec<[F; 3]> {
        let mut out = vec![[0.0 as F; 3]; ntotat];
        for partial in partials {
            for i in 0..ntotat {
                out[i][0] += partial[i][0];
                out[i][1] += partial[i][1];
                out[i][2] += partial[i][2];
            }
        }
        out
    }

    #[test]
    fn merge_serial_matches_reference() {
        let partials = make_partials(8, 1023);
        let expected = reference_sum(&partials, 1023);
        let mut got = vec![[0.0 as F; 3]; 1023];
        WorkBuffers::merge_partials_serial(&partials, &mut got);
        assert_eq!(got, expected);
    }

    /// Parallel merge must be bit-identical to the serial path: the
    /// partitioning is contiguous and each atom is touched by exactly
    /// one thread, so there is no cross-thread accumulation order.
    #[cfg(feature = "rayon")]
    #[test]
    fn merge_parallel_matches_serial_bitwise() {
        // Sizes chosen to straddle both the threshold (4096) and the
        // parallel chunk size (256), exercising full chunks, partial
        // chunks, and one chunk-per-thread edge case.
        for &ntotat in &[1usize, 255, 256, 257, 4_095, 4_096, 8_193] {
            let partials = make_partials(6, ntotat);
            let mut got_serial = vec![[0.0 as F; 3]; ntotat];
            let mut got_parallel = vec![[0.0 as F; 3]; ntotat];
            WorkBuffers::merge_partials_serial(&partials, &mut got_serial);
            WorkBuffers::merge_partials_parallel(&partials, &mut got_parallel);
            assert_eq!(
                got_serial, got_parallel,
                "parallel merge diverged at ntotat={ntotat}"
            );
        }
    }

    /// An empty partial list is a no-op: `out` stays as the caller
    /// initialized it. Guards against a panic in the parallel iterator
    /// when the inner `for partial in partials` loop sees zero items.
    #[test]
    fn merge_empty_partials_leaves_out_unchanged() {
        let partials: Vec<Vec<[F; 3]>> = vec![];
        let mut out = vec![[1.0, 2.0, 3.0]; 128];
        let expected = out.clone();
        WorkBuffers::merge_partials_serial(&partials, &mut out);
        assert_eq!(out, expected);
        #[cfg(feature = "rayon")]
        {
            let mut out = vec![[1.0, 2.0, 3.0]; 128];
            WorkBuffers::merge_partials_parallel(&partials, &mut out);
            assert_eq!(out, expected);
        }
    }
}
