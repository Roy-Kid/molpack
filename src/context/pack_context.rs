//! Packing runtime context — mirrors Packmol `compute_data` behavior.

use std::sync::Arc;

use crate::restraint::{AtomRestraint, Restraint};
use molrs::core::CellGrid;
use molrs::core::Element;
use molrs::core::SimBox;
use molrs::op::F;
use ndarray::array;

use super::geometry::GeometryKey;
use super::work_buffers::WorkBuffers;

/// `flags` bit for a fixed-structure atom inside [`AtomProps`].
pub const ATOM_FLAG_FIXED: u32 = 1 << 0;
/// `flags` bit for an atom that participates in the optional short-radius
/// penalty inside [`AtomProps`].
pub const ATOM_FLAG_SHORT: u32 = 1 << 1;

/// Sentinel for "no entry" in the linked-cell lists (`latomfirst`,
/// `latomnext`, `latomfix`, `lcellfirst`, `lcellnext`). Using `u32::MAX`
/// instead of `Option<usize>` shrinks each slot from 16 B (unoptimized
/// `Option<usize>`) to 4 B, which dominates the pair-kernel cache
/// footprint where these arrays are traversed per atom visit. Valid
/// indices must satisfy `idx < NONE_IDX`; `PackContext::new` debug-asserts
/// this for `ntotat` and `ncell_total`.
pub const NONE_IDX: u32 = u32::MAX;

/// Compact AoS view of every atom's **hot** pair-kernel inputs.
///
/// The pair kernel in `objective::fparc` / `gparc` / `fgparc` reads the
/// molecule identity, radii, `fscale`, and a flag byte every time it
/// looks at another atom. Scattering those across seven separate `Vec<_>`s
/// turns each visit into seven independent cache-line fetches. Packing
/// them into one 40-byte struct (8-byte aligned) collapses that to at
/// most two cache-line fetches, with a matching load for `xcart[jcart]`
/// kept separate because positions change every evaluation.
///
/// Layout (40 bytes, natural 8-byte alignment — two atoms fit 80 bytes,
/// i.e. ~1.25 cache lines so most pair visits hit a single line):
/// - `ibmol, ibtype` — 2 × u32 (8 bytes) — same-molecule skip.
/// - `fscale, radius, radius_ini` — 3 × F (24 bytes) — distance kernel.
/// - `flags` — u32 (4 bytes) — `ATOM_FLAG_FIXED`, `ATOM_FLAG_SHORT`.
/// - private padding — u32 (4 bytes) to keep the struct multiple of 8.
///
/// Cold-path fields (`short_radius`, `short_radius_scale`) stay in
/// separate `Vec<F>`s on `PackContext` so the common "no short radius"
/// workload does not pay to load them.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct AtomProps {
    pub ibmol: u32,
    pub ibtype: u32,
    pub fscale: F,
    pub radius: F,
    pub radius_ini: F,
    pub flags: u32,
    /// 8-byte alignment padding. Private so the struct's size remains a
    /// layout detail: flipping `F` or adding fields recomputes size at
    /// compile time (see the size assertion below) without churning the
    /// public API.
    _padding: u32,
}

// Compile-time assertion that `AtomProps` stays 40 bytes (`F = f64`).
// Shrinking it again without noticing would regress the pair-kernel
// cache-line budget; growing it past 48 bytes would cost an extra
// fetch per atom visit. The check is always compiled (not
// `#[cfg(test)]`) so release builds catch layout drift too.
const _: () = assert!(
    std::mem::size_of::<AtomProps>() == 40,
    "AtomProps must stay 40 bytes"
);

/// Default quadratic-penalty scale (`scale2`) applied when no caller override
/// is supplied — the packer seeds [`PackContext`] with it, and post-pack
/// validation scores penalties on the same scale.
pub(crate) const DEFAULT_SCALE2: F = 0.01;

/// Full runtime context for one packing execution.
/// All arrays are 0-based; Fortran 1-based arrays are shifted by -1.
pub struct PackContext {
    // ---- Atom Cartesian coordinates (updated every function evaluation) ----
    /// Current Cartesian positions: `xcart[icart]` = `[x, y, z]`. Size: ntotat.
    pub xcart: Vec<[F; 3]>,
    /// Element per atom: `elements[icart]`. Size: ntotat. `None` means unknown/"X".
    pub elements: Vec<Option<Element>>,

    // ---- Reference (centered) coordinates ----
    /// Reference conformer per **copy**: `coor[icart]` = `[x, y, z]`.
    ///
    /// Shares `xcart`'s index space exactly (type-major, copy-major,
    /// atom-minor, free types then fixed types), so `icart` addresses both.
    /// Copies of one type start identical; in-loop optimizers
    /// (module `crate::optimizer`) relax each copy
    /// independently, after which they diverge. Size: `ntotat`.
    pub coor: Vec<[F; 3]>,

    // ---- Radii ----
    /// Current radii (may be scaled): `radius[icart]`. Size: ntotat.
    pub radius: Vec<F>,
    /// Original (unscaled) radii: `radius_ini[icart]`. Size: ntotat.
    pub radius_ini: Vec<F>,
    /// Function scaling per atom: `fscale[icart]`. Size: ntotat.
    pub fscale: Vec<F>,

    // ---- Short radius (optional secondary penalty) ----
    pub use_short_radius: Vec<bool>,
    pub short_radius: Vec<F>,
    pub short_radius_scale: Vec<F>,
    /// Summary flag — `true` iff any atom has `use_short_radius` set. Lets
    /// the objective hot loop skip the short-radius branch entirely for the
    /// common case (no short-radius usage). Maintained incrementally via
    /// setters + re-synced by [`Self::sync_atom_props`].
    pub any_short_radius: bool,
    /// Summary flag — `true` iff any atom is a fixed-structure atom. Lets
    /// the objective hot loop skip the `fixedatom[i] && fixedatom[j]`
    /// short-circuit when there are no fixed atoms. Maintained
    /// incrementally via setters + re-synced by [`Self::sync_atom_props`].
    pub any_fixed_atoms: bool,
    /// Incremental counter driving `any_fixed_atoms`. Private because
    /// the `any_*` flag is the observable contract; setters keep both
    /// counters and flag consistent.
    n_fixed_atoms: usize,
    /// Incremental counter driving `any_short_radius`.
    n_short_radius: usize,

    /// AoS mirror of the frequently-read per-atom fields. Kept in sync
    /// with the individual `Vec<_>`s by [`Self::sync_atom_props`] plus the
    /// per-field setters (`set_radius`, `set_fscale`, `set_fixed_atom`,
    /// `set_use_short_radius`, `set_ibmol`, `set_ibtype`). The objective
    /// hot kernels read exclusively from here; callers mutating the
    /// underlying `Vec<_>`s directly must call [`Self::sync_atom_props`]
    /// before the next `evaluate()`, or the debug-build invariant in
    /// [`Self::debug_assert_atom_props_sync`] will catch the drift.
    pub atom_props: Vec<AtomProps>,

    // ---- Objective function accumulators ----
    /// Maximum inter-molecular distance violation (fdist in Fortran).
    pub fdist: F,
    /// Maximum constraint violation (frest in Fortran).
    pub frest: F,
    /// Per-atom distance violation (for movebad).
    pub fdist_atom: Vec<F>,
    /// Per-atom constraint violation (for movebad).
    pub frest_atom: Vec<F>,

    // ---- Molecule topology ----
    /// Number of molecules per type: `nmols[itype]`. 0-based type index.
    pub nmols: Vec<usize>,
    /// Number of atoms per type: `natoms[itype]`. 0-based type index.
    pub natoms: Vec<usize>,
    /// First atom index (0-based) of each type's first copy: `idfirst[itype]`.
    /// Base into both [`Self::coor`] and [`Self::xcart`]; copy `imol` of type
    /// `itype` starts at `idfirst[itype] + imol * natoms[itype]`.
    pub idfirst: Vec<usize>,
    /// Total number of types (free).
    pub ntype: usize,
    /// Total number of types including fixed types.
    pub ntype_with_fixed: usize,
    /// Total number of free molecules.
    pub ntotmol: usize,
    /// Total number of atoms (free + fixed).
    pub ntotat: usize,
    /// Number of fixed atoms.
    pub nfixedat: usize,

    // ---- Rotation constraints (Packmol constrain_rotation) ----
    /// Rotation constraint flags per free type in Euler variable order
    /// [beta(y), gama(z), teta(x)].
    pub constrain_rot: Vec<[bool; 3]>,
    /// Rotation bounds per free type and Euler variable:
    /// [center_rad, half_width_rad].
    pub rot_bound: Vec<[[F; 2]; 3]>,

    // ---- Restraints ----
    /// All restraints pool: `restraints[irest]`.
    pub restraints: Vec<Arc<dyn AtomRestraint>>,
    /// CSR offsets for per-atom restraint indices:
    /// restraints of atom `icart` are in `iratom_data[iratom_offsets[icart]..iratom_offsets[icart+1]]`.
    pub iratom_offsets: Vec<usize>,
    /// Flattened per-atom restraint indices.
    pub iratom_data: Vec<usize>,
    /// Group-level restraints, paired with the (0-based) type they act on:
    /// `(itype, restraint)`. Evaluated once per group in the objective with the
    /// coordinates of all copies of `itype`; the coupled gradient is scattered
    /// back into `work.gxcar`. Empty in the common (no collective restraint) case.
    pub collective: Vec<(usize, Arc<dyn Restraint>)>,

    // ---- Cell list bookkeeping ----
    /// Type index per atom: `ibtype[icart]` (0-based type index).
    pub ibtype: Vec<usize>,
    /// Molecule index within its type: `ibmol[icart]` (0-based).
    pub ibmol: Vec<usize>,
    /// Is this a fixed atom?
    pub fixedatom: Vec<bool>,
    /// Is this type being computed in the current iteration?
    pub comptype: Vec<bool>,

    // ---- Cell geometry ----
    /// The packing cell: lattice matrix, origin and per-axis periodicity.
    ///
    /// Replaces the axis-aligned `pbc_min` / `pbc_length` pair, so the cell may
    /// be hexagonal, monoclinic or fully triclinic — the geometry Packmol
    /// cannot express at all. It is also the single source of truth for the
    /// minimum image: the pair kernel calls
    /// [`SimBox::shortest_vector_impl`], which honours `pbc` per axis.
    pub simbox: SimBox,
    /// Partition of [`simbox`](Self::simbox) into cells, in fractional space.
    ///
    /// Wraps on periodic axes and clamps on non-periodic ones, so an atom
    /// pushed outside the cell mid-optimisation lands in the nearest edge cell
    /// instead of on the opposite face.
    pub grid: CellGrid,

    // ---- Linked cell lists ----
    /// `latomfirst[icell]` = first atom index in cell, `NONE_IDX` if empty.
    /// Stored as flat Vec indexed by `index_cell`.
    pub latomfirst: Vec<u32>,
    /// `latomnext[icart]` = next atom in the same cell (`NONE_IDX` = end).
    pub latomnext: Vec<u32>,
    /// Fixed atom list per cell (permanent), `NONE_IDX` if cell has no fixed atoms.
    pub latomfix: Vec<u32>,
    /// Occupied cell linked list: first cell (`NONE_IDX` if none occupied).
    pub lcellfirst: u32,
    /// `lcellnext[icell]` = next occupied cell (`NONE_IDX` = end).
    pub lcellnext: Vec<u32>,
    /// Is cell empty?
    pub empty_cell: Vec<bool>,
    /// Cells that contain fixed atoms and must be restored on every reset.
    pub fixed_cells: Vec<usize>,
    /// Cells touched during the previous objective/gradient evaluation.
    pub active_cells: Vec<usize>,
    /// Forward-neighbour cells per cell, flattened (CSR): cell `i` owns
    /// `neighbor_cells[neighbor_start[i]..neighbor_start[i + 1]]`. Read it
    /// through [`neighbors`](Self::neighbors).
    ///
    /// Every unordered pair of adjacent cells appears exactly once across a
    /// full sweep, so both the objective and the gradient walk this one
    /// half-stencil and neither needs a full 26-neighbour list.
    ///
    /// Storage is variable-length because forwardness is decided by cell
    /// **index** (`nc > i`), not by a fixed set of 13 offset directions. The
    /// direction-based scheme Packmol uses double-counts as soon as a periodic
    /// axis holds two cells — cell 0's `+1` neighbour is cell 1, and cell 1's
    /// `+1` wraps back to cell 0 — which is common in a thin slab or a
    /// flat triclinic cell. Per-cell counts then vary from 0 to 26 while the
    /// total stays at 13 per cell.
    pub neighbor_cells: Vec<u32>,
    /// CSR offsets into [`neighbor_cells`](Self::neighbor_cells); length is
    /// `n_cells + 1`.
    pub neighbor_start: Vec<u32>,

    // ---- State flags ----
    /// If true, skip pair-distance computations (constraints only during init).
    pub init1: bool,
    /// If true, accumulate per-atom fdist/frest (movebad mode).
    pub move_flag: bool,
    /// Run the pair-kernel reductions (`accumulate_pair_f`,
    /// `accumulate_pair_fg`) on rayon. Off by default — parallelism is
    /// an explicit opt-in via [`PackEngine::with_parallel_eval`](crate::PackEngine::with_parallel_eval) because the
    /// crossover is workload-shaped and can't be inferred reliably from
    /// `active_cells.len()`. The flag is stored regardless of the
    /// `rayon` feature so the engine API stays the same; when the
    /// feature is off the field is read but the parallel path doesn't
    /// exist and the serial branch runs unconditionally.
    pub parallel_pair_eval: bool,

    // ---- Algorithm parameters ----
    pub scale: F,
    pub scale2: F,

    // ---- Bounding box ----
    pub sizemin: [F; 3],
    pub sizemax: [F; 3],

    // ---- Maximum internal distances per type ----
    pub dmax: Vec<F>,

    // ---- Work buffers ----
    pub work: WorkBuffers,

    // ---- Debug: call counters (zeroed per pgencan call) ----
    ncf: usize,
    ncg: usize,
}

impl PackContext {
    /// Allocate and zero-initialize all arrays.
    pub fn new(ntotat: usize, ntotmol: usize, ntype: usize) -> Self {
        let simbox =
            SimBox::cube(1.0, array![0.0, 0.0, 0.0], [false; 3]).expect("unit placeholder cell");
        let grid = CellGrid::with_dims([1; 3], [false; 3]);
        let ncell_total = grid.n_cells();
        debug_assert!(
            ntotat < NONE_IDX as usize,
            "ntotat={ntotat} must fit in u32 (< NONE_IDX)"
        );
        Self {
            xcart: vec![[0.0; 3]; ntotat],
            elements: vec![None; ntotat],
            coor: Vec::new(),
            radius: vec![0.0; ntotat],
            radius_ini: vec![0.0; ntotat],
            fscale: vec![1.0; ntotat],
            use_short_radius: vec![false; ntotat],
            short_radius: vec![0.0; ntotat],
            short_radius_scale: vec![0.0; ntotat],
            any_short_radius: false,
            any_fixed_atoms: false,
            n_fixed_atoms: 0,
            n_short_radius: 0,
            atom_props: vec![AtomProps::default(); ntotat],
            fdist: 0.0,
            frest: 0.0,
            fdist_atom: vec![0.0; ntotat],
            frest_atom: vec![0.0; ntotat],
            nmols: Vec::new(),
            natoms: Vec::new(),
            idfirst: Vec::new(),
            ntype,
            ntype_with_fixed: ntype,
            ntotmol,
            ntotat,
            nfixedat: 0,
            constrain_rot: vec![[false; 3]; ntype],
            rot_bound: vec![[[0.0; 2]; 3]; ntype],
            restraints: Vec::new(),
            iratom_offsets: vec![0; ntotat + 1],
            iratom_data: Vec::new(),
            collective: Vec::new(),
            ibtype: vec![0; ntotat],
            ibmol: vec![0; ntotat],
            fixedatom: vec![false; ntotat],
            comptype: vec![true; ntype],
            simbox,
            grid,
            latomfirst: vec![NONE_IDX; ncell_total],
            latomnext: vec![NONE_IDX; ntotat],
            latomfix: vec![NONE_IDX; ncell_total],
            lcellfirst: NONE_IDX,
            lcellnext: vec![NONE_IDX; ncell_total],
            empty_cell: vec![true; ncell_total],
            fixed_cells: Vec::new(),
            active_cells: Vec::new(),
            neighbor_cells: Vec::new(),
            neighbor_start: vec![0; ncell_total + 1],
            init1: false,
            move_flag: false,
            parallel_pair_eval: false,
            scale: 1.0,
            scale2: DEFAULT_SCALE2,
            sizemin: [0.0; 3],
            sizemax: [0.0; 3],
            dmax: vec![0.0; ntype],
            work: WorkBuffers::new(ntotat),
            ncf: 0,
            ncg: 0,
        }
    }

    /// Resize cell list arrays after ncells is set.
    pub fn resize_cell_arrays(&mut self) {
        let nc = self.grid.n_cells();
        debug_assert!(
            nc < NONE_IDX as usize,
            "ncell_total={nc} must fit in u32 (< NONE_IDX)"
        );
        self.latomfirst = vec![NONE_IDX; nc];
        self.latomfix = vec![NONE_IDX; nc];
        self.lcellnext = vec![NONE_IDX; nc];
        self.empty_cell = vec![true; nc];
        self.fixed_cells.clear();
        self.active_cells.clear();
        self.rebuild_neighbor_cells();
    }

    /// Reset cell lists (called at start of each compute_f/compute_g).
    /// Port of `resetcells.f90`.
    pub fn resetcells(&mut self) {
        self.lcellfirst = NONE_IDX;
        for &icell in &self.active_cells {
            self.latomfirst[icell] = NONE_IDX;
            self.lcellnext[icell] = NONE_IDX;
            self.empty_cell[icell] = true;
        }
        self.active_cells.clear();

        for &icell in &self.fixed_cells {
            self.latomfirst[icell] = self.latomfix[icell];
            self.empty_cell[icell] = false;
            self.lcellnext[icell] = self.lcellfirst;
            self.lcellfirst = icell as u32;
            self.active_cells.push(icell);
        }

        // Reset latomnext for free atoms only
        let free_atoms = self.ntotat - self.nfixedat;
        self.latomnext[..free_atoms].fill(NONE_IDX);
    }

    #[inline]
    pub fn reset_eval_counters(&mut self) {
        self.ncf = 0;
        self.ncg = 0;
    }

    /// Rebuild `atom_props` from the individual per-atom `Vec<_>`s. Call
    /// once after packer setup has populated every per-atom field, and
    /// whenever a field has been mutated via the underlying `Vec<_>`
    /// directly rather than through a setter.
    ///
    /// Also refreshes the summary flags (`any_fixed_atoms`,
    /// `any_short_radius`) and their backing counters.
    pub fn sync_atom_props(&mut self) {
        let n = self.ntotat;
        if self.atom_props.len() != n {
            self.atom_props.resize(n, AtomProps::default());
        }
        let mut n_fixed = 0usize;
        let mut n_short = 0usize;
        for i in 0..n {
            let fixed = self.fixedatom[i];
            let use_short = self.use_short_radius[i];
            if fixed {
                n_fixed += 1;
            }
            if use_short {
                n_short += 1;
            }
            let mut flags = 0u32;
            if fixed {
                flags |= ATOM_FLAG_FIXED;
            }
            if use_short {
                flags |= ATOM_FLAG_SHORT;
            }
            self.atom_props[i] = AtomProps {
                ibmol: self.ibmol[i] as u32,
                ibtype: self.ibtype[i] as u32,
                flags,
                _padding: 0,
                fscale: self.fscale[i],
                radius: self.radius[i],
                radius_ini: self.radius_ini[i],
            };
        }
        self.n_fixed_atoms = n_fixed;
        self.n_short_radius = n_short;
        self.any_fixed_atoms = n_fixed > 0;
        self.any_short_radius = n_short > 0;
    }

    /// Update atom `i`'s live radius on both the `Vec<F>` and the AoS
    /// mirror.  Preferred over writing `sys.radius[i]` directly: the
    /// hot-loop kernel reads `atom_props[i].radius`, so a raw write
    /// would silently desynchronize.
    #[inline]
    pub fn set_radius(&mut self, i: usize, value: F) {
        self.radius[i] = value;
        if i < self.atom_props.len() {
            self.atom_props[i].radius = value;
        }
    }

    /// Update atom `i`'s live `fscale` on both storages. Used by the
    /// scaling-phase paths that modulate per-atom weight.
    #[inline]
    pub fn set_fscale(&mut self, i: usize, value: F) {
        self.fscale[i] = value;
        if i < self.atom_props.len() {
            self.atom_props[i].fscale = value;
        }
    }

    /// Drop the cached Cartesian expansion so the next evaluation rebuilds it.
    ///
    /// The cache is keyed on `x` (COM / Euler), `comptype` and the cell
    /// geometry — **not** on [`coor`](Self::coor). Anything that mutates the
    /// reference conformers while leaving `x` alone (an in-loop optimizer, say)
    /// must call this, or the next `evaluate` at the same `x` returns the
    /// pre-mutation objective and any comparison against it is meaningless.
    #[inline]
    pub fn invalidate_geometry_cache(&mut self) {
        self.work.cached_geometry = None;
    }

    /// Toggle atom `i`'s fixed-structure flag and keep the `atom_props`
    /// mirror, the `any_fixed_atoms` summary flag, and the private
    /// counter in lock-step.
    #[inline]
    pub fn set_fixed_atom(&mut self, i: usize, is_fixed: bool) {
        let was_fixed = self.fixedatom[i];
        if was_fixed == is_fixed {
            return;
        }
        self.fixedatom[i] = is_fixed;
        if i < self.atom_props.len() {
            let flags = &mut self.atom_props[i].flags;
            if is_fixed {
                *flags |= ATOM_FLAG_FIXED;
            } else {
                *flags &= !ATOM_FLAG_FIXED;
            }
        }
        if is_fixed {
            self.n_fixed_atoms += 1;
        } else {
            self.n_fixed_atoms -= 1;
        }
        self.any_fixed_atoms = self.n_fixed_atoms > 0;
    }

    /// Toggle atom `i`'s `use_short_radius` flag, maintaining the mirror,
    /// summary flag, and counter.
    #[inline]
    pub fn set_use_short_radius(&mut self, i: usize, use_short: bool) {
        let was = self.use_short_radius[i];
        if was == use_short {
            return;
        }
        self.use_short_radius[i] = use_short;
        if i < self.atom_props.len() {
            let flags = &mut self.atom_props[i].flags;
            if use_short {
                *flags |= ATOM_FLAG_SHORT;
            } else {
                *flags &= !ATOM_FLAG_SHORT;
            }
        }
        if use_short {
            self.n_short_radius += 1;
        } else {
            self.n_short_radius -= 1;
        }
        self.any_short_radius = self.n_short_radius > 0;
    }

    /// Update `ibmol[i]` and the matching mirror field.
    #[inline]
    pub fn set_ibmol(&mut self, i: usize, value: usize) {
        self.ibmol[i] = value;
        if i < self.atom_props.len() {
            self.atom_props[i].ibmol = value as u32;
        }
    }

    /// Update `ibtype[i]` and the matching mirror field.
    #[inline]
    pub fn set_ibtype(&mut self, i: usize, value: usize) {
        self.ibtype[i] = value;
        if i < self.atom_props.len() {
            self.atom_props[i].ibtype = value as u32;
        }
    }

    /// Debug-only invariant: every `atom_props[i]` matches the state
    /// derivable from the individual per-atom `Vec<_>`s, and the
    /// summary counters / flags agree with reality.
    ///
    /// O(ntotat) per call in debug builds; compiled to a single
    /// early-return in release builds (see the `cfg!(debug_assertions)`
    /// gate below — with `#[inline(always)]` the release body DCE's).
    /// The objective hot loop calls this at the entry of `compute_f` /
    /// `compute_g` / `compute_fg` so direct-write drift fires at the
    /// next evaluate instead of silently producing wrong energies.
    #[inline(always)]
    pub fn debug_assert_atom_props_sync(&self) {
        if !cfg!(debug_assertions) {
            return;
        }
        let n = self.ntotat;
        assert_eq!(
            self.atom_props.len(),
            n,
            "atom_props length {} != ntotat {} — call sync_atom_props after a resize",
            self.atom_props.len(),
            n
        );
        let mut n_fixed = 0usize;
        let mut n_short = 0usize;
        for i in 0..n {
            let ap = &self.atom_props[i];
            let expected_fixed = self.fixedatom[i];
            let expected_short = self.use_short_radius[i];
            let expected_flags = if expected_fixed { ATOM_FLAG_FIXED } else { 0 }
                | if expected_short { ATOM_FLAG_SHORT } else { 0 };
            if expected_fixed {
                n_fixed += 1;
            }
            if expected_short {
                n_short += 1;
            }
            assert_eq!(
                ap.ibmol, self.ibmol[i] as u32,
                "atom_props[{i}].ibmol drift: mirror={} vec={}",
                ap.ibmol, self.ibmol[i]
            );
            assert_eq!(
                ap.ibtype, self.ibtype[i] as u32,
                "atom_props[{i}].ibtype drift"
            );
            assert_eq!(
                ap.fscale, self.fscale[i],
                "atom_props[{i}].fscale drift — did you write sys.fscale[{i}] directly?"
            );
            assert_eq!(
                ap.radius, self.radius[i],
                "atom_props[{i}].radius drift — use set_radius()"
            );
            assert_eq!(
                ap.radius_ini, self.radius_ini[i],
                "atom_props[{i}].radius_ini drift"
            );
            assert_eq!(
                ap.flags, expected_flags,
                "atom_props[{i}].flags drift — did you write sys.fixedatom/use_short_radius directly?"
            );
        }
        assert_eq!(
            self.n_fixed_atoms, n_fixed,
            "n_fixed_atoms counter drift: stored={} derived={}",
            self.n_fixed_atoms, n_fixed
        );
        assert_eq!(
            self.n_short_radius, n_short,
            "n_short_radius counter drift: stored={} derived={}",
            self.n_short_radius, n_short
        );
        assert_eq!(self.any_fixed_atoms, n_fixed > 0);
        assert_eq!(self.any_short_radius, n_short > 0);
    }

    #[inline]
    pub fn increment_ncf(&mut self) {
        self.ncf += 1;
    }

    #[inline]
    pub fn increment_ncg(&mut self) {
        self.ncg += 1;
    }

    #[inline]
    pub fn ncf(&self) -> usize {
        self.ncf
    }

    #[inline]
    pub fn ncg(&self) -> usize {
        self.ncg
    }

    /// Recompute the forward-neighbour table from the cell partition.
    ///
    /// Delegates the stencil to [`CellGrid::stencil_forward`], so periodicity,
    /// small-`celldim` aliasing and deduplication are decided in one place
    /// rather than re-derived here.
    fn rebuild_neighbor_cells(&mut self) {
        let nc = self.grid.n_cells();
        self.neighbor_start.clear();
        self.neighbor_start.reserve(nc + 1);
        self.neighbor_cells.clear();

        let mut buf = [0usize; 27];
        for icell in 0..nc {
            self.neighbor_start.push(self.neighbor_cells.len() as u32);
            let n = self.grid.stencil_forward(icell, &mut buf);
            self.neighbor_cells
                .extend(buf[..n].iter().map(|&c| c as u32));
        }
        self.neighbor_start.push(self.neighbor_cells.len() as u32);
    }

    /// Forward neighbours of `icell` — see [`neighbor_cells`](Self::neighbor_cells).
    #[inline(always)]
    pub fn neighbors(&self, icell: usize) -> &[u32] {
        let lo = self.neighbor_start[icell] as usize;
        let hi = self.neighbor_start[icell + 1] as usize;
        &self.neighbor_cells[lo..hi]
    }

    /// Forward neighbours of `icell` copied into a caller-owned buffer.
    ///
    /// The serial pair loops mutate the context while walking the neighbour
    /// list, so they cannot hold a borrow of it. Copying into a stack array
    /// keeps that allocation-free — the same thing the fixed `[usize; 13]`
    /// table gave for free when it was `Copy`.
    #[inline(always)]
    pub fn copy_neighbors(&self, icell: usize, out: &mut [u32; 27]) -> usize {
        let nbs = self.neighbors(icell);
        out[..nbs.len()].copy_from_slice(nbs);
        nbs.len()
    }

    /// Number of cells along each lattice direction.
    #[inline(always)]
    pub fn ncells(&self) -> [usize; 3] {
        self.grid.celldim().map(|d| d as usize)
    }

    /// Compact identity of the packing geometry, for the evaluation cache.
    ///
    /// Everything the cell list depends on: the partition and the lattice it
    /// partitions. Comparing this is what lets a repeated evaluation at the
    /// same coordinates reuse the previous cell assignment.
    pub fn geometry_key(&self) -> GeometryKey {
        let h = self.simbox.h_view();
        let o = self.simbox.origin_view();
        GeometryKey {
            celldim: self.grid.celldim(),
            pbc: self.grid.pbc(),
            h: [
                h[[0, 0]],
                h[[0, 1]],
                h[[0, 2]],
                h[[1, 0]],
                h[[1, 1]],
                h[[1, 2]],
                h[[2, 0]],
                h[[2, 1]],
                h[[2, 2]],
            ],
            origin: [o[0], o[1], o[2]],
        }
    }
}

#[cfg(test)]
mod atom_props_tests {
    use super::*;

    fn tiny_ctx(ntotat: usize) -> PackContext {
        let mut sys = PackContext::new(ntotat, ntotat, 1);
        for i in 0..ntotat {
            sys.ibmol[i] = i;
            sys.ibtype[i] = 0;
            sys.radius[i] = 1.0;
            sys.radius_ini[i] = 1.0;
            sys.fscale[i] = 1.0;
        }
        sys.sync_atom_props();
        sys
    }

    #[test]
    fn atom_props_size_is_40_bytes_on_f64() {
        // Runtime echo of the compile-time size assertion — cheap and also
        // readable as failing test output.
        assert_eq!(std::mem::size_of::<AtomProps>(), 40);
        assert_eq!(std::mem::align_of::<AtomProps>(), 8);
    }

    #[test]
    fn sync_atom_props_populates_mirror_and_flags() {
        let mut sys = PackContext::new(3, 3, 1);
        sys.ibmol = vec![10, 20, 30];
        sys.ibtype = vec![1, 2, 3];
        sys.fscale = vec![0.5, 0.25, 0.125];
        sys.radius = vec![1.1, 2.2, 3.3];
        sys.radius_ini = vec![1.0, 2.0, 3.0];
        sys.fixedatom = vec![false, true, false];
        sys.use_short_radius = vec![false, false, true];
        sys.sync_atom_props();

        for i in 0..3 {
            assert_eq!(sys.atom_props[i].ibmol, sys.ibmol[i] as u32);
            assert_eq!(sys.atom_props[i].ibtype, sys.ibtype[i] as u32);
            assert_eq!(sys.atom_props[i].fscale, sys.fscale[i]);
            assert_eq!(sys.atom_props[i].radius, sys.radius[i]);
            assert_eq!(sys.atom_props[i].radius_ini, sys.radius_ini[i]);
        }
        assert_eq!(sys.atom_props[0].flags, 0);
        assert_eq!(sys.atom_props[1].flags, ATOM_FLAG_FIXED);
        assert_eq!(sys.atom_props[2].flags, ATOM_FLAG_SHORT);
        assert!(sys.any_fixed_atoms);
        assert!(sys.any_short_radius);
        sys.debug_assert_atom_props_sync();
    }

    #[test]
    fn set_radius_keeps_mirror_in_sync() {
        let mut sys = tiny_ctx(4);
        sys.set_radius(2, 7.25);
        assert_eq!(sys.radius[2], 7.25);
        assert_eq!(sys.atom_props[2].radius, 7.25);
        // Other atoms unchanged.
        assert_eq!(sys.atom_props[0].radius, 1.0);
        assert_eq!(sys.atom_props[3].radius, 1.0);
        sys.debug_assert_atom_props_sync();
    }

    #[test]
    fn set_fscale_keeps_mirror_in_sync() {
        let mut sys = tiny_ctx(4);
        sys.set_fscale(1, 0.125);
        assert_eq!(sys.fscale[1], 0.125);
        assert_eq!(sys.atom_props[1].fscale, 0.125);
        sys.debug_assert_atom_props_sync();
    }

    #[test]
    fn set_fixed_atom_updates_mirror_flag_counter_and_summary() {
        let mut sys = tiny_ctx(3);
        assert!(!sys.any_fixed_atoms);
        assert_eq!(sys.n_fixed_atoms, 0);

        sys.set_fixed_atom(1, true);
        assert!(sys.fixedatom[1]);
        assert_eq!(sys.atom_props[1].flags & ATOM_FLAG_FIXED, ATOM_FLAG_FIXED);
        assert_eq!(sys.n_fixed_atoms, 1);
        assert!(sys.any_fixed_atoms);
        sys.debug_assert_atom_props_sync();

        // Setting the same value again is a no-op — counter must not re-increment.
        sys.set_fixed_atom(1, true);
        assert_eq!(sys.n_fixed_atoms, 1);

        sys.set_fixed_atom(0, true);
        assert_eq!(sys.n_fixed_atoms, 2);
        sys.set_fixed_atom(1, false);
        assert_eq!(sys.n_fixed_atoms, 1);
        assert!(sys.any_fixed_atoms);
        sys.set_fixed_atom(0, false);
        assert_eq!(sys.n_fixed_atoms, 0);
        assert!(!sys.any_fixed_atoms);
        sys.debug_assert_atom_props_sync();
    }

    #[test]
    fn set_use_short_radius_updates_mirror_flag_counter_and_summary() {
        let mut sys = tiny_ctx(3);
        sys.set_use_short_radius(2, true);
        assert!(sys.use_short_radius[2]);
        assert_eq!(sys.atom_props[2].flags & ATOM_FLAG_SHORT, ATOM_FLAG_SHORT);
        assert_eq!(sys.n_short_radius, 1);
        assert!(sys.any_short_radius);
        sys.debug_assert_atom_props_sync();

        sys.set_use_short_radius(2, false);
        assert_eq!(sys.n_short_radius, 0);
        assert!(!sys.any_short_radius);
        assert_eq!(sys.atom_props[2].flags & ATOM_FLAG_SHORT, 0);
        sys.debug_assert_atom_props_sync();
    }

    #[test]
    fn set_fixed_and_short_flags_coexist_on_same_atom() {
        let mut sys = tiny_ctx(2);
        sys.set_fixed_atom(0, true);
        sys.set_use_short_radius(0, true);
        assert_eq!(
            sys.atom_props[0].flags,
            ATOM_FLAG_FIXED | ATOM_FLAG_SHORT,
            "both flags must combine without clobbering each other"
        );
        sys.set_fixed_atom(0, false);
        assert_eq!(
            sys.atom_props[0].flags, ATOM_FLAG_SHORT,
            "clearing FIXED must leave SHORT intact"
        );
        sys.debug_assert_atom_props_sync();
    }

    #[test]
    fn set_ibmol_and_set_ibtype_keep_mirror_in_sync() {
        let mut sys = tiny_ctx(3);
        sys.set_ibmol(1, 42);
        assert_eq!(sys.ibmol[1], 42);
        assert_eq!(sys.atom_props[1].ibmol, 42);
        sys.set_ibtype(2, 7);
        assert_eq!(sys.ibtype[2], 7);
        assert_eq!(sys.atom_props[2].ibtype, 7);
        sys.debug_assert_atom_props_sync();
    }

    /// Regression guard: the pre-Fix-1 state had two paths that wrote
    /// `sys.fixedatom[i]` directly (Packmol-faithful but redundant).
    /// Should one ever be reintroduced and fail to re-sync, the
    /// debug-build invariant must catch it. This test flips the flag
    /// directly on the `Vec<bool>` and confirms the assertion panics.
    ///
    /// Gated on `debug_assertions` because the underlying invariant is
    /// a `debug_assert!` — release builds compile it out and the
    /// `#[should_panic]` would never fire.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "atom_props")]
    fn debug_invariant_catches_direct_fixedatom_write() {
        let mut sys = tiny_ctx(2);
        sys.fixedatom[0] = true; // bypass setter — simulates the bug class
        sys.debug_assert_atom_props_sync();
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "atom_props")]
    fn debug_invariant_catches_direct_fscale_write() {
        let mut sys = tiny_ctx(2);
        sys.fscale[1] = 99.0;
        sys.debug_assert_atom_props_sync();
    }

    /// Counter consistency under mixed operations: `sync_atom_props`
    /// and the setters must produce identical `n_fixed_atoms` /
    /// `n_short_radius` values. This guards against a setter
    /// forgetting to increment/decrement.
    #[test]
    fn counters_match_sync_after_mixed_mutations() {
        let mut sys = tiny_ctx(10);
        for i in [0usize, 3, 7] {
            sys.set_fixed_atom(i, true);
        }
        for i in [2usize, 5] {
            sys.set_use_short_radius(i, true);
        }
        let pre_n_fixed = sys.n_fixed_atoms;
        let pre_n_short = sys.n_short_radius;

        // Re-sync from Vecs and compare — counters must match.
        sys.sync_atom_props();
        assert_eq!(sys.n_fixed_atoms, pre_n_fixed);
        assert_eq!(sys.n_short_radius, pre_n_short);
        assert_eq!(sys.n_fixed_atoms, 3);
        assert_eq!(sys.n_short_radius, 2);
    }
}

#[cfg(test)]
mod neighbor_table_tests {
    use super::*;

    fn ctx_with_grid(celldim: [u32; 3], pbc: [bool; 3]) -> PackContext {
        let mut sys = PackContext::new(1, 1, 1);
        sys.simbox = SimBox::cube(10.0, array![0.0, 0.0, 0.0], pbc).expect("cell");
        sys.grid = CellGrid::with_dims(celldim, pbc);
        sys.resize_cell_arrays();
        sys
    }

    /// Every unordered pair of adjacent cells appears exactly once, so the
    /// table holds 13 entries per cell on a fully periodic grid — the same
    /// total the fixed 13-offset table carried.
    #[test]
    fn periodic_grid_holds_thirteen_forward_neighbours_per_cell() {
        let sys = ctx_with_grid([4, 4, 4], [true; 3]);
        let n = sys.grid.n_cells();
        assert_eq!(sys.neighbor_cells.len(), 13 * n);
    }

    /// What the move from offset-based to index-based forwardness actually
    /// changed: the *distribution*. Cell 0 sees all 26 of its neighbours as
    /// forward, the last cell sees none. Total work is unchanged, but it is no
    /// longer flat across cells, which is what a rayon-over-cells traversal
    /// divides up.
    #[test]
    fn forward_counts_are_uneven_while_the_total_is_not() {
        let sys = ctx_with_grid([4, 4, 4], [true; 3]);
        let counts: Vec<usize> = (0..sys.grid.n_cells())
            .map(|i| sys.neighbors(i).len())
            .collect();
        assert_eq!(counts.iter().sum::<usize>(), 13 * counts.len());
        assert_eq!(*counts.iter().max().expect("non-empty"), 26);
        assert_eq!(*counts.iter().min().expect("non-empty"), 0);
    }

    /// A non-periodic axis has no wrap-around neighbours, so the table is
    /// smaller than the periodic case rather than padded with far-side cells
    /// the way an unconditional wrap would leave it.
    #[test]
    fn a_non_periodic_axis_drops_its_wrap_neighbours() {
        let periodic = ctx_with_grid([4, 4, 4], [true; 3]);
        let confined = ctx_with_grid([4, 4, 4], [true, true, false]);
        assert!(
            confined.neighbor_cells.len() < periodic.neighbor_cells.len(),
            "confining an axis must remove neighbour entries, got {} vs {}",
            confined.neighbor_cells.len(),
            periodic.neighbor_cells.len()
        );
    }

    /// Two cells on an axis: the aliasing case that double-counted under a
    /// fixed `{0, +1}` offset set.
    #[test]
    fn two_cells_on_an_axis_are_paired_once() {
        let sys = ctx_with_grid([2, 1, 1], [true; 3]);
        assert_eq!(sys.neighbors(0), &[1]);
        assert_eq!(sys.neighbors(1), &[] as &[u32]);
    }
}

#[cfg(test)]
mod geometry_cache_tests {
    #![allow(clippy::needless_range_loop)]
    //! Tests for the geometry cache fast path in `compute_f` / `compute_fg` /
    //! `compute_g`. The cache is hit when the caller evaluates at the same `x`
    //! (and identical comptype / cell grid) as the previous call — in that case
    //! the Cartesian expansion and cell-list rebuild are skipped and only the
    //! pair / constraint kernels re-run on the stored state.
    //!
    //! These tests assert the cache path produces bit-identical results to the
    //! fresh-rebuild path across several call sequences used by the packer.

    use std::sync::Arc;

    use crate::PackContext;
    use crate::objective::{compute_f, compute_fg};
    use crate::testutil::inside_box;
    use molrs::op::F;

    // ── molrs regions lifted to "stay inside" (the one geometric restraint) ─────

    // ── setup helpers (mirror restraint::geometric::tests::gradient patterns) ──────────────────────

    fn setup_cells(sys: &mut PackContext, cell_n: usize, cell_len: F) {
        let side = cell_len * cell_n as F;
        sys.simbox =
            molrs::core::SimBox::cube(side, molrs::op::F3::zeros(3), [false; 3]).expect("cell");
        sys.grid = molrs::core::CellGrid::with_dims([cell_n as u32; 3], [false; 3]);
        sys.resize_cell_arrays();
    }

    /// Three single-atom molecules inside a 5³ box with a pair-overlap setup.
    fn mixed_system() -> (PackContext, Vec<F>) {
        let mut sys = PackContext::new(3, 3, 1);
        sys.ntype_with_fixed = 1;
        sys.nmols = vec![3];
        sys.natoms = vec![1];
        sys.idfirst = vec![0];
        sys.comptype = vec![true];
        // `coor` holds one reference conformer **per copy**, sharing `xcart`'s
        // index space — three single-atom copies, so three entries.
        sys.coor = vec![[0.0, 0.0, 0.0]; 3];

        sys.radius = vec![1.0; 3];
        sys.radius_ini = vec![1.0; 3];
        sys.fscale = vec![1.0; 3];
        sys.ibmol = vec![0, 1, 2];
        sys.sync_atom_props();

        sys.restraints = vec![Arc::new(inside_box([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]))];
        sys.iratom_offsets = vec![0, 1, 2, 3];
        sys.iratom_data = vec![0, 0, 0];

        setup_cells(&mut sys, 1, 10.0);

        // x = [com0(3), com1(3), com2(3), euler0(3), euler1(3), euler2(3)]
        let x = vec![
            6.0, 2.0, 2.0, // com0: outside box on +x, forces restraint penalty
            3.0, 2.0, 2.0, // com1: close to com2, forces pair penalty
            3.5, 2.5, 2.0, // com2
            0.1, 0.2, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, 0.4,
        ];
        (sys, x)
    }

    fn force_cache_miss(sys: &mut PackContext) {
        sys.work.cached_geometry = None;
    }

    // ── compute_f cache ────────────────────────────────────────────────────────

    #[test]
    fn compute_f_cached_path_matches_fresh() {
        let (mut sys_a, x) = mixed_system();
        let (mut sys_b, _) = mixed_system();

        // Prime both systems: one cache, one fresh each call.
        let f_a1 = compute_f(&x, &mut sys_a);
        let fdist_a1 = sys_a.fdist;
        let frest_a1 = sys_a.frest;

        force_cache_miss(&mut sys_b);
        let f_b1 = compute_f(&x, &mut sys_b);
        let fdist_b1 = sys_b.fdist;
        let frest_b1 = sys_b.frest;

        assert_eq!(
            f_a1.to_bits(),
            f_b1.to_bits(),
            "compute_f first call must agree bitwise"
        );
        assert_eq!(fdist_a1.to_bits(), fdist_b1.to_bits());
        assert_eq!(frest_a1.to_bits(), frest_b1.to_bits());

        // Second call at same x: a hits cache, b forced miss.
        let f_a2 = compute_f(&x, &mut sys_a);
        force_cache_miss(&mut sys_b);
        let f_b2 = compute_f(&x, &mut sys_b);

        assert_eq!(f_a2.to_bits(), f_b2.to_bits());
        assert_eq!(f_a2.to_bits(), f_a1.to_bits(), "cache must be pure");
        assert_eq!(sys_a.fdist.to_bits(), sys_b.fdist.to_bits());
        assert_eq!(sys_a.frest.to_bits(), sys_b.frest.to_bits());
    }

    // ── compute_fg cache ───────────────────────────────────────────────────────

    #[test]
    fn compute_fg_cached_path_matches_fresh() {
        let (mut sys_a, x) = mixed_system();
        let (mut sys_b, _) = mixed_system();

        let mut g_a1 = vec![0.0; x.len()];
        let mut g_b1 = vec![0.0; x.len()];
        let f_a1 = compute_fg(&x, &mut sys_a, &mut g_a1);
        force_cache_miss(&mut sys_b);
        let f_b1 = compute_fg(&x, &mut sys_b, &mut g_b1);

        assert_eq!(f_a1.to_bits(), f_b1.to_bits());
        for i in 0..x.len() {
            assert_eq!(
                g_a1[i].to_bits(),
                g_b1[i].to_bits(),
                "compute_fg first call: g[{i}] mismatch {} vs {}",
                g_a1[i],
                g_b1[i]
            );
        }

        // Second call at same x — a hits cache, b forced miss.
        let mut g_a2 = vec![0.0; x.len()];
        let mut g_b2 = vec![0.0; x.len()];
        let f_a2 = compute_fg(&x, &mut sys_a, &mut g_a2);
        force_cache_miss(&mut sys_b);
        let f_b2 = compute_fg(&x, &mut sys_b, &mut g_b2);

        assert_eq!(f_a2.to_bits(), f_b2.to_bits());
        assert_eq!(f_a2.to_bits(), f_a1.to_bits(), "cache must be pure");
        for i in 0..x.len() {
            assert_eq!(g_a2[i].to_bits(), g_b2[i].to_bits());
            assert_eq!(g_a1[i].to_bits(), g_a2[i].to_bits());
        }
    }

    // ── cross-mode cache reuse (compute_fg → compute_f at same x) ──────────────

    #[test]
    fn compute_f_reuses_compute_fg_geometry() {
        let (mut sys_a, x) = mixed_system();
        let (mut sys_b, _) = mixed_system();

        // A: warm with compute_fg then call compute_f — cache hit expected.
        let mut g_a = vec![0.0; x.len()];
        let _ = compute_fg(&x, &mut sys_a, &mut g_a);
        let f_a = compute_f(&x, &mut sys_a);

        // B: always fresh.
        let mut g_b = vec![0.0; x.len()];
        force_cache_miss(&mut sys_b);
        let _ = compute_fg(&x, &mut sys_b, &mut g_b);
        force_cache_miss(&mut sys_b);
        let f_b = compute_f(&x, &mut sys_b);

        assert_eq!(f_a.to_bits(), f_b.to_bits());
        assert_eq!(sys_a.fdist.to_bits(), sys_b.fdist.to_bits());
        assert_eq!(sys_a.frest.to_bits(), sys_b.frest.to_bits());
    }

    // ── packer's "unscaled re-evaluation" pattern ─────────────────────────────
    //
    // After `pgencan` converges, `packer.rs` swaps `radius := radius_ini` and calls
    // `compute_f` at the same `x` to measure violations under the true (unscaled)
    // atomic radii. The cache key intentionally does not include `radius`, so this
    // pattern hits the cache — verify the result under a radius mutation between
    // the scaled and unscaled calls is identical to a fresh full evaluation.

    #[test]
    fn radii_swap_between_fg_and_f_cached_matches_fresh() {
        let (mut sys_a, x) = mixed_system();
        let (mut sys_b, _) = mixed_system();

        // Scaled radii (typical during packing: discale=1.2).
        let scaled: Vec<F> = sys_a.radius_ini.iter().map(|r| r * 1.2).collect();
        let unscaled = sys_a.radius_ini.clone();

        // --- A: cached path ---
        sys_a.radius = scaled.clone();
        sys_a.sync_atom_props();
        let mut g_a = vec![0.0; x.len()];
        let _ = compute_fg(&x, &mut sys_a, &mut g_a);

        sys_a.radius = unscaled.clone();
        sys_a.sync_atom_props();
        let f_a = compute_f(&x, &mut sys_a);

        // --- B: always rebuild ---
        sys_b.radius = scaled.clone();
        sys_b.sync_atom_props();
        let mut g_b = vec![0.0; x.len()];
        force_cache_miss(&mut sys_b);
        let _ = compute_fg(&x, &mut sys_b, &mut g_b);

        sys_b.radius = unscaled.clone();
        sys_b.sync_atom_props();
        force_cache_miss(&mut sys_b);
        let f_b = compute_f(&x, &mut sys_b);

        assert_eq!(f_a.to_bits(), f_b.to_bits());
        assert_eq!(sys_a.fdist.to_bits(), sys_b.fdist.to_bits());
        assert_eq!(sys_a.frest.to_bits(), sys_b.frest.to_bits());
    }

    // ── move_flag path must stay on the slow path ─────────────────────────────
    //
    // When `move_flag` is true (during movebad), per-atom `fdist_atom` /
    // `frest_atom` are accumulated inside the pair / constraint kernels. Running
    // the cache path a second time would double-count; so cache must be bypassed
    // whenever `move_flag` is set.

    #[test]
    fn move_flag_true_bypasses_cache() {
        let (mut sys, x) = mixed_system();

        // Warm the cache at normal (move_flag=false) state.
        let _ = compute_f(&x, &mut sys);
        assert!(sys.work.cached_geometry.is_some());

        // Now turn on move_flag and reset per-atom trackers.
        sys.move_flag = true;
        sys.fdist_atom.iter_mut().for_each(|v| *v = 0.0);
        sys.frest_atom.iter_mut().for_each(|v| *v = 0.0);
        let _ = compute_f(&x, &mut sys);
        let fdist_move_a = sys.fdist_atom.clone();
        let frest_move_a = sys.frest_atom.clone();

        // Fresh context, same sequence, cache forced off each call.
        let (mut sys2, _) = mixed_system();
        force_cache_miss(&mut sys2);
        let _ = compute_f(&x, &mut sys2);
        sys2.move_flag = true;
        sys2.fdist_atom.iter_mut().for_each(|v| *v = 0.0);
        sys2.frest_atom.iter_mut().for_each(|v| *v = 0.0);
        force_cache_miss(&mut sys2);
        let _ = compute_f(&x, &mut sys2);

        assert_eq!(
            fdist_move_a.len(),
            sys2.fdist_atom.len(),
            "fdist_atom shape mismatch"
        );
        for (i, (&a, &b)) in fdist_move_a.iter().zip(sys2.fdist_atom.iter()).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "fdist_atom[{i}] mismatch under move_flag: {a} vs {b}"
            );
        }
        for (i, (&a, &b)) in frest_move_a.iter().zip(sys2.frest_atom.iter()).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "frest_atom[{i}] mismatch under move_flag: {a} vs {b}"
            );
        }
    }
}
