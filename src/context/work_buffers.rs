//! Reusable temporary buffers for objective/gradient and movebad paths.

use super::geometry::GeometryKey;
use molrs::types::F;

/// Reusable mutable buffers shared across packing iterations.
pub struct WorkBuffers {
    /// Cartesian gradient accumulator used by objective gradient evaluation.
    pub gxcar: Vec<[F; 3]>,
    /// Per-rayon-worker scratch gradient buffers for the parallel pair-gradient
    /// path. Flat, thread-major: worker `t` owns the contiguous region
    /// `[t*ntotat .. (t+1)*ntotat)`. Each worker accumulates its half-stencil
    /// pair forces here race-free (its region is touched by no other
    /// concurrently-running task), then the regions are reduced into [`gxcar`].
    /// Reusing one persistent buffer keeps the hot path allocation-free; it is
    /// zeroed (in parallel) at the start of each gradient evaluation. Sized
    /// lazily to `nthreads * ntotat`.
    ///
    /// [`gxcar`]: Self::gxcar
    #[cfg(feature = "rayon")]
    pub grad_partials: Vec<[F; 3]>,
    /// Reused per-active-molecule descriptor list `(itype, icart0, ilubar,
    /// ilugan)` for the phase-structured `expand_molecules` /
    /// `project_cartesian_gradient` passes. Both rebuild it from the current
    /// `comptype` each call (cheap index arithmetic, no trig); persisting the
    /// `Vec` keeps that hot rebuild allocation-free.
    pub mol_descs: Vec<(usize, usize, usize, usize)>,
    /// Reused `(icart0, len, natoms_per_copy)` span list for the collective
    /// restraints, rebuilt from the current `comptype` on every evaluation.
    /// Persisting the `Vec` keeps `accumulate_collective_fg` allocation-free
    /// on a path GENCAN calls thousands of times.
    pub collective_spans: Vec<(usize, usize, usize)>,
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
    /// Geometry the cached cell assignment was built for; `None` invalidates
    /// the cache.
    pub cached_geometry: Option<GeometryKey>,
}

impl WorkBuffers {
    pub fn new(ntotat: usize) -> Self {
        Self {
            gxcar: vec![[0.0; 3]; ntotat],
            #[cfg(feature = "rayon")]
            grad_partials: Vec::new(),
            mol_descs: Vec::new(),
            collective_spans: Vec::new(),
            radiuswork: vec![0.0; ntotat],
            fmol: Vec::new(),
            flash_ind: Vec::new(),
            flash_l: Vec::new(),
            cached_x: Vec::new(),
            cached_comptype: Vec::new(),
            cached_init1: false,
            cached_geometry: None,
        }
    }

    pub fn matches_cached_geometry(
        &self,
        x: &[F],
        comptype: &[bool],
        init1: bool,
        geometry: GeometryKey,
    ) -> bool {
        self.cached_geometry == Some(geometry)
            && self.cached_init1 == init1
            && self.cached_x == x
            && self.cached_comptype == comptype
    }

    pub fn update_cached_geometry(
        &mut self,
        x: &[F],
        comptype: &[bool],
        init1: bool,
        geometry: GeometryKey,
    ) {
        self.cached_x.clear();
        self.cached_x.extend_from_slice(x);
        self.cached_comptype.clear();
        self.cached_comptype.extend_from_slice(comptype);
        self.cached_init1 = init1;
        self.cached_geometry = Some(geometry);
    }
}
