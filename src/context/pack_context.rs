//! Packing runtime context — mirrors Packmol `compute_data` behavior.

use std::sync::Arc;

use crate::cell::{cell_ind, icell_to_cell, index_cell};
use crate::constraints::{Constraints, EvalMode, EvalOutput};
use crate::restraint::Restraint;
use molrs::Element;
use molrs::types::F;

use super::model::ModelData;
use super::state::{RuntimeState, RuntimeStateMut};
use super::work_buffers::WorkBuffers;

/// Index of a restraint assigned to a specific atom.
pub type RestraintRef = usize;

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
    /// compile time (see [`ATOM_PROPS_SIZE`] below) without churning the
    /// public API.
    _padding: u32,
}

/// Compile-time assertion that [`AtomProps`] stays 40 bytes (`F = f64`).
/// Shrinking it again without noticing would regress the pair-kernel
/// cache-line budget; growing it past 48 bytes would cost an extra
/// fetch per atom visit. The check is always compiled (not
/// `#[cfg(test)]`) so release builds catch layout drift too.
pub const ATOM_PROPS_SIZE: usize = 40;
const _ATOM_PROPS_IS_40_BYTES: [(); ATOM_PROPS_SIZE] = [(); std::mem::size_of::<AtomProps>()];

/// Neighbor offsets used by `computef.f90` (13 forward neighbors).
const NEIGHBOR_OFFSETS_F: [(isize, isize, isize); 13] = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, -1, 0),
    (1, 0, -1),
    (0, 1, -1),
    (0, 1, 1),
    (1, 1, 0),
    (1, 0, 1),
    (1, -1, -1),
    (1, -1, 1),
    (1, 1, -1),
    (1, 1, 1),
];

/// Neighbor offsets used by `computeg.f90` (13 forward neighbors, different order).
const NEIGHBOR_OFFSETS_G: [(isize, isize, isize); 13] = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (0, 1, 1),
    (0, 1, -1),
    (1, 1, 0),
    (1, 0, 1),
    (1, -1, 0),
    (1, 0, -1),
    (1, 1, 1),
    (1, 1, -1),
    (1, -1, 1),
    (1, -1, -1),
];

/// Full runtime context for one packing execution.
/// All arrays are 0-based; Fortran 1-based arrays are shifted by -1.
pub struct PackContext {
    // ---- Constraints facade ----
    pub constraints: Constraints,

    // ---- Atom Cartesian coordinates (updated every function evaluation) ----
    /// Current Cartesian positions: `xcart[icart]` = `[x, y, z]`. Size: ntotat.
    pub xcart: Vec<[F; 3]>,
    /// Element per atom: `elements[icart]`. Size: ntotat. `None` means unknown/"X".
    pub elements: Vec<Option<Element>>,

    // ---- Reference (centered) coordinates ----
    /// Reference coordinates `coor[idatom]` = `[x, y, z]`. Size: total atoms across all types.
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
    /// First datum atom index (0-based) for each type: `idfirst[itype]`.
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
    pub restraints: Vec<Arc<dyn Restraint>>,
    /// CSR offsets for per-atom restraint indices:
    /// restraints of atom `icart` are in `iratom_data[iratom_offsets[icart]..iratom_offsets[icart+1]]`.
    pub iratom_offsets: Vec<usize>,
    /// Flattened per-atom restraint indices.
    pub iratom_data: Vec<RestraintRef>,

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
    pub ncells: [usize; 3],
    pub cell_length: [F; 3],
    pub pbc_length: [F; 3],
    pub pbc_min: [F; 3],
    /// Per-axis periodicity flags. `pbc_periodic[k] == true` means axis
    /// `k` wraps in the pair-kernel minimum image and the cell list;
    /// `false` means the cell list clamps and no wrap is applied. Set
    /// from restraints that override `Restraint::periodic_box()`.
    pub pbc_periodic: [bool; 3],

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
    /// Precomputed 13 forward-neighbor cell indices per cell for `compute_f`.
    pub neighbor_cells_f: Vec<[usize; 13]>,
    /// Precomputed 13 forward-neighbor cell indices per cell for `compute_g`.
    pub neighbor_cells_g: Vec<[usize; 13]>,

    // ---- State flags ----
    /// If true, skip pair-distance computations (constraints only during init).
    pub init1: bool,
    /// If true, accumulate per-atom fdist/frest (movebad mode).
    pub move_flag: bool,
    /// Run the pair-kernel reductions (`accumulate_pair_f`,
    /// `accumulate_pair_fg`) on rayon. Off by default — parallelism is
    /// an explicit opt-in via [`Molpack::with_parallel_eval`](crate::Molpack::with_parallel_eval) because the
    /// crossover is workload-shaped and can't be inferred reliably from
    /// `active_cells.len()`. The flag is stored regardless of the
    /// `rayon` feature so the `Molpack` API stays the same; when the
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

    // ---- Output frame (owned, built incrementally) ----
    /// Frame that accumulates constant columns (element, mol_id) during init
    /// and receives position columns at the end of packing.
    pub frame: molrs::Frame,

    // ---- Debug: call counters (zeroed per pgencan call) ----
    ncf: usize,
    ncg: usize,
}

impl PackContext {
    /// Allocate and zero-initialize all arrays.
    pub fn new(ntotat: usize, ntotmol: usize, ntype: usize) -> Self {
        let ncells = [1, 1, 1];
        let ncell_total = ncells[0] * ncells[1] * ncells[2];
        debug_assert!(
            ntotat < NONE_IDX as usize,
            "ntotat={ntotat} must fit in u32 (< NONE_IDX)"
        );
        Self {
            constraints: Constraints,
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
            ibtype: vec![0; ntotat],
            ibmol: vec![0; ntotat],
            fixedatom: vec![false; ntotat],
            comptype: vec![true; ntype],
            ncells,
            cell_length: [1.0; 3],
            pbc_length: [1.0; 3],
            pbc_min: [0.0; 3],
            pbc_periodic: [false; 3],
            latomfirst: vec![NONE_IDX; ncell_total],
            latomnext: vec![NONE_IDX; ntotat],
            latomfix: vec![NONE_IDX; ncell_total],
            lcellfirst: NONE_IDX,
            lcellnext: vec![NONE_IDX; ncell_total],
            empty_cell: vec![true; ncell_total],
            fixed_cells: Vec::new(),
            active_cells: Vec::new(),
            neighbor_cells_f: vec![[0; 13]; ncell_total],
            neighbor_cells_g: vec![[0; 13]; ncell_total],
            init1: false,
            move_flag: false,
            parallel_pair_eval: false,
            scale: 1.0,
            scale2: 0.01,
            sizemin: [0.0; 3],
            sizemax: [0.0; 3],
            dmax: vec![0.0; ntype],
            work: WorkBuffers::new(ntotat),
            frame: molrs::Frame::new(),
            ncf: 0,
            ncg: 0,
        }
    }

    /// Context view for mostly static model data.
    #[inline]
    pub fn model(&self) -> ModelData<'_> {
        ModelData { ctx: self }
    }

    /// Read-only runtime state view.
    #[inline]
    pub fn runtime(&self) -> RuntimeState<'_> {
        RuntimeState { ctx: self }
    }

    /// Mutable runtime state view.
    #[inline]
    pub fn runtime_mut(&mut self) -> RuntimeStateMut<'_> {
        RuntimeStateMut { ctx: self }
    }

    /// Unified constraints evaluation entrypoint.
    #[inline]
    pub fn evaluate(&mut self, x: &[F], mode: EvalMode, gradient: Option<&mut [F]>) -> EvalOutput {
        let constraints = self.constraints;
        constraints.evaluate(x, self, mode, gradient)
    }

    /// Resize cell list arrays after ncells is set.
    pub fn resize_cell_arrays(&mut self) {
        let nc = self.ncells[0] * self.ncells[1] * self.ncells[2];
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
        self.neighbor_cells_f = vec![[0; 13]; nc];
        self.neighbor_cells_g = vec![[0; 13]; nc];
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

    fn rebuild_neighbor_cells(&mut self) {
        let (nx, ny, nz) = (self.ncells[0], self.ncells[1], self.ncells[2]);
        let nc = nx * ny * nz;
        for icell in 0..nc {
            let cell = icell_to_cell(icell, &self.ncells);
            let (ci, cj, ck) = (cell[0], cell[1], cell[2]);

            let mut nbs_f = [0usize; 13];
            for (idx, &(di, dj, dk)) in NEIGHBOR_OFFSETS_F.iter().enumerate() {
                let ncell = [
                    cell_ind(ci as isize + di, nx),
                    cell_ind(cj as isize + dj, ny),
                    cell_ind(ck as isize + dk, nz),
                ];
                nbs_f[idx] = index_cell(&ncell, &self.ncells);
            }
            self.neighbor_cells_f[icell] = nbs_f;

            let mut nbs_g = [0usize; 13];
            for (idx, &(di, dj, dk)) in NEIGHBOR_OFFSETS_G.iter().enumerate() {
                let ncell = [
                    cell_ind(ci as isize + di, nx),
                    cell_ind(cj as isize + dj, ny),
                    cell_ind(ck as isize + dk, nz),
                ];
                nbs_g[idx] = index_cell(&ncell, &self.ncells);
            }
            self.neighbor_cells_g[icell] = nbs_g;
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
        // Runtime echo of the compile-time `_ATOM_PROPS_IS_40_BYTES`
        // assertion — cheap and also readable as failing test output.
        assert_eq!(std::mem::size_of::<AtomProps>(), 40);
        assert_eq!(std::mem::align_of::<AtomProps>(), 8);
        assert_eq!(ATOM_PROPS_SIZE, 40);
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
