//! Main packing orchestration.
//! Port of the outer loop in `app/packmol.f90`.

use std::sync::Arc;

use molrs::Element;
use molrs::types::F;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::constraints::EvalMode;
use crate::context::PackContext;
use crate::error::PackError;
use crate::euler::{compcart, eulerfixed};
use crate::gencan::{GencanParams, GencanWorkspace, pgencan};
use crate::handler::{Handler, PhaseInfo, StepInfo};
use crate::initial::{SwapState, init_xcart_from_x, initial};
use crate::movebad::{MoveBadConfig, movebad};
use crate::numerics::objective_small_floor;
use crate::relaxer::RelaxerRunner;
use crate::restraint::Restraint;
use crate::target::{CenteringMode, Target};

/// Result of a packing run.
///
/// The `frame` contains an "atoms" block with x, y, z, element, mol_id
/// columns — moved from the packing context (zero-copy ownership transfer).
#[derive(Debug, Clone)]
pub struct PackResult {
    /// Atoms frame with x, y, z (f32), element (String), mol_id (i64).
    pub frame: molrs::Frame,
    /// Maximum inter-molecular distance violation at termination.
    pub fdist: F,
    /// Maximum constraint violation at termination.
    pub frest: F,
    /// Whether the packing converged (`fdist < precision && frest < precision`).
    pub converged: bool,
}

impl PackResult {
    /// Extract atom positions as `Vec<[F; 3]>` (SoA→AoS conversion).
    pub fn positions(&self) -> Vec<[F; 3]> {
        let atoms = self.frame.get("atoms").expect("frame has no 'atoms' block");
        let x = atoms.get_float("x").expect("no 'x' column");
        let y = atoms.get_float("y").expect("no 'y' column");
        let z = atoms.get_float("z").expect("no 'z' column");
        x.iter()
            .zip(y.iter())
            .zip(z.iter())
            .map(|((&xi, &yi), &zi)| [xi, yi, zi])
            .collect()
    }

    /// Number of atoms in the result.
    #[inline]
    pub fn natoms(&self) -> usize {
        self.frame.get("atoms").and_then(|b| b.nrows()).unwrap_or(0)
    }
}

/// Default Packmol parameters
const PRECISION: F = 0.01;
// Packmol default from getinp.f90: discale = 1.1d0
const DISCALE: F = 1.1;
/// Fixed GENCAN inner iteration limit (Packmol default maxit = 20).
const GENCAN_MAXIT: usize = 20;
/// Packmol default sidemax (getinp.f90).
const SIDEMAX: F = 1000.0;
/// Packmol default movefrac.
const MOVEFRAC: F = 0.05;
/// Default minimum atom-atom distance tolerance (Packmol's `dism` default = 2.0 Å).
/// Atom radii are set to `tolerance / 2` for all atoms, matching Packmol's
/// `radius(i) = dism/2.d0` (packmol.f90 line 283).
const DEFAULT_TOLERANCE: F = 2.0;

/// The packer.
pub struct Molpack {
    handlers: Vec<Box<dyn Handler>>,
    /// Global restraints — broadcast to every target at `pack()` time
    /// (semantic equivalence to calling `target.with_restraint(r.clone())`
    /// on every target; no separate global-storage code path).
    global_restraints: Vec<Arc<dyn Restraint>>,
    precision: F,
    discale: F,
    /// Minimum atom-atom distance (Packmol's `tolerance`/`dism`). Default 2.0 Å.
    /// Atom radii = `tolerance / 2`.
    tolerance: F,
    /// GENCAN inner iterations (`maxit` keyword).
    inner_iterations: usize,
    /// Initialization outer loops (`nloop0` keyword). `None` means Packmol default (20*ntype).
    init_passes: Option<usize>,
    /// Maximum system half-size used in initial restmol stage (`sidemax` keyword).
    init_box_half_size: F,
    /// Fraction of molecules perturbed when packing stalls (Packmol's `movefrac`).
    perturb_fraction: F,
    /// Randomize perturbation target selection (Packmol's `movebadrandom`).
    random_perturb: bool,
    /// Master switch for the stall-perturbation heuristic (inverts Packmol's
    /// `disable_movebad`: `true` = perturb enabled, `false` = disabled).
    perturb: bool,
    /// Seed for the internal RNG. Default `0` (deterministic).
    seed: u64,
    /// Run the pair-kernel reductions on rayon. Off by default: see
    /// [`with_parallel_eval`][Self::with_parallel_eval].
    parallel_eval: bool,
    /// Global periodic-boundary box, as set by
    /// [`with_periodic_box`][Self::with_periodic_box] (Packmol `pbc`
    /// keyword). When both this and a restraint-declared PBC are
    /// present, they must match exactly or `pack()` returns
    /// [`PackError::ConflictingPeriodicBoxes`].
    periodic_box: Option<PeriodicSpec>,
}

impl Default for Molpack {
    fn default() -> Self {
        Self::new()
    }
}

impl Molpack {
    /// Create a packer with default settings and no handlers.
    ///
    /// All tuning knobs are set via `with_*` methods below; they all have
    /// defensible defaults, so `Molpack::new().pack(&targets, max_loops)`
    /// is a valid invocation. The one argument that has no default is
    /// `max_loops` — it depends on system size and convergence difficulty,
    /// so it lives on the terminal [`pack`][Self::pack] call, not on the
    /// builder.
    pub fn new() -> Self {
        Self {
            handlers: Vec::new(),
            global_restraints: Vec::new(),
            precision: PRECISION,
            discale: DISCALE,
            tolerance: DEFAULT_TOLERANCE,
            inner_iterations: GENCAN_MAXIT,
            init_passes: None,
            init_box_half_size: SIDEMAX,
            perturb_fraction: MOVEFRAC,
            random_perturb: false,
            perturb: true,
            seed: 0,
            parallel_eval: false,
            periodic_box: None,
        }
    }

    /// Append a progress handler. Multiple handlers compose in call order.
    pub fn with_handler(mut self, h: impl Handler + 'static) -> Self {
        self.handlers.push(Box::new(h));
        self
    }

    /// Append a **global** restraint — applied to every atom of every
    /// target at `pack()` time.
    ///
    /// Semantic equivalence (scope law):
    /// ```text
    /// molpack.with_global_restraint(r)
    ///   ≡ for each target: target.with_restraint(r.clone())
    /// ```
    ///
    /// Implementation mirrors the equivalence — no separate "global
    /// restraint" storage path in `PackContext`; the restraint is cloned
    /// into each target's `molecule_restraints` when `pack()` is invoked.
    pub fn with_global_restraint(mut self, r: impl Restraint + 'static) -> Self {
        self.global_restraints.push(Arc::new(r));
        self
    }

    /// Convergence precision for `fdist` and `frest` (default `0.01`).
    pub fn with_precision(mut self, p: F) -> Self {
        self.precision = p;
        self
    }

    /// Minimum atom-atom distance tolerance (default `2.0 Å`).
    /// Atom radii are set to `tolerance / 2`.
    pub fn with_tolerance(mut self, t: F) -> Self {
        self.tolerance = t;
        self
    }

    /// GENCAN inner iteration count (default `20`; Packmol `maxit`).
    pub fn with_inner_iterations(mut self, n: usize) -> Self {
        self.inner_iterations = n;
        self
    }

    /// Initialization outer-loop passes (Packmol `nloop0`).
    /// `0` restores the Packmol default of `20 * ntype`.
    pub fn with_init_passes(mut self, n: usize) -> Self {
        self.init_passes = if n == 0 { None } else { Some(n) };
        self
    }

    /// Maximum half-size of the initial placement box (default `1000.0`;
    /// Packmol `sidemax`).
    pub fn with_init_box_half_size(mut self, f: F) -> Self {
        self.init_box_half_size = f;
        self
    }

    /// Declare a global periodic-boundary box (Packmol `pbc`). Every axis
    /// is treated as periodic. When set, the packer's cell grid is built
    /// from `max - min`, bypassing the fallback that derives a box from
    /// post-Phase-1 atom positions — which can be ±`sidemax` wide when
    /// the script has no spatial constraints and drives `ncells` to
    /// 10⁸+ cells.
    ///
    /// If any restraint also declares a `periodic_box()`, the two must
    /// match exactly (bounds + flags) or `pack()` returns
    /// [`PackError::ConflictingPeriodicBoxes`].
    pub fn with_periodic_box(mut self, min: [F; 3], max: [F; 3]) -> Self {
        self.periodic_box = Some((min, max, [true; 3]));
        self
    }

    /// Fraction of molecules re-sampled when packing stalls (default `0.05`;
    /// Packmol `movefrac`).
    pub fn with_perturb_fraction(mut self, f: F) -> Self {
        self.perturb_fraction = f;
        self
    }

    /// Randomize perturbation target selection (default `false`;
    /// Packmol `movebadrandom`).
    pub fn with_random_perturb(mut self, enabled: bool) -> Self {
        self.random_perturb = enabled;
        self
    }

    /// Enable the stall-perturbation heuristic (default `true`;
    /// inverts Packmol's `disable_movebad`). Pass `false` to disable.
    pub fn with_perturb(mut self, enabled: bool) -> Self {
        self.perturb = enabled;
        self
    }

    /// Seed for the internal RNG (default `0` — deterministic).
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    /// Run the pair-kernel reductions on rayon worker threads (default
    /// `false`).
    ///
    /// Parallelism is opt-in because the crossover where rayon pays off
    /// is **workload-shaped, not size-shaped** — the per-call parallel
    /// speedup measured on isolated `compute_fg` doesn't predict
    /// end-to-end `Molpack::pack` wall clock (task-dispatch overhead
    /// accumulates across thousands of calls per pack; the perturbation
    /// pass skips the parallel path entirely; real workloads often
    /// *regress* even when `active_cells.len()` clears any naive
    /// threshold). The user knows their workload shape; the library
    /// doesn't.
    ///
    /// Compile with `--features rayon` for this flag to have any
    /// effect. Without the feature the pair kernel is serial
    /// unconditionally and this setting is silently ignored.
    pub fn with_parallel_eval(mut self, enabled: bool) -> Self {
        self.parallel_eval = enabled;
        self
    }

    /// Run the packing.
    ///
    /// `max_loops` is the outer iteration budget; it is positional
    /// because there is no defensible default (the right value depends
    /// on system size and convergence difficulty). Every other knob
    /// lives on the builder.
    ///
    /// Returns a [`PackResult`] containing the final atom positions and
    /// convergence information (`fdist`, `frest`, `converged`).
    pub fn pack(&mut self, targets: &[Target], max_loops: usize) -> Result<PackResult, PackError> {
        if targets.is_empty() {
            return Err(PackError::NoTargets);
        }

        for (i, t) in targets.iter().enumerate() {
            if t.natoms() == 0 {
                return Err(PackError::EmptyMolecule(i));
            }
        }

        // Scope equivalence (spec §4): broadcast global restraints to each target.
        // If none were added via `Molpack::add_restraint`, this is a no-op.
        let broadcast_targets: Vec<Target>;
        let targets: &[Target] = if self.global_restraints.is_empty() {
            targets
        } else {
            broadcast_targets = targets
                .iter()
                .map(|t| {
                    let mut t = t.clone();
                    for r in &self.global_restraints {
                        t.molecule_restraints.push(Arc::clone(r));
                    }
                    t
                })
                .collect();
            &broadcast_targets
        };
        // Derive the system periodic box from two independent sources:
        //   1. `Molpack::with_periodic_box` — global PBC (script `pbc`).
        //   2. `Restraint::periodic_box` — per-restraint declarations,
        //      e.g. `InsideBoxRestraint` with any axis marked periodic.
        //
        // The two must not disagree; if both are present they must
        // match exactly (same bounds, same per-axis flags). Validate
        // the global PBC extent here since `derive_periodic_box`
        // already does this for restraint-sourced boxes.
        if let Some((min, max, _)) = self.periodic_box {
            let length = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
            if length.iter().any(|&v| v <= 0.0) {
                return Err(PackError::InvalidPBCBox { min, max });
            }
        }
        let pbc = match (self.periodic_box, derive_periodic_box(targets)?) {
            (None, derived) => derived,
            (Some(global), None) => Some(global),
            (Some(global), Some(derived)) if global == derived => Some(global),
            (Some(global), Some(derived)) => {
                return Err(PackError::ConflictingPeriodicBoxes {
                    first: global,
                    second: derived,
                });
            }
        };

        let mut rng = SmallRng::seed_from_u64(self.seed);

        // Split into free and fixed targets
        let free_targets: Vec<&Target> = targets.iter().filter(|t| t.fixed_at.is_none()).collect();
        let fixed_targets: Vec<&Target> = targets.iter().filter(|t| t.fixed_at.is_some()).collect();

        let ntype = free_targets.len();
        let ntype_with_fixed = ntype + fixed_targets.len();

        // Count atoms
        let ntotmol_free: usize = free_targets.iter().map(|t| t.count).sum();
        let ntotat_free: usize = free_targets.iter().map(|t| t.count * t.natoms()).sum();
        let ntotat_fixed: usize = fixed_targets.iter().map(|t| t.natoms()).sum();
        let ntotat = ntotat_free + ntotat_fixed;

        // Variable count: 3N COM + 3N Euler angles (only free molecules)
        let n = 6 * ntotmol_free;

        // Build PackContext
        let mut sys = PackContext::new(ntotat, ntotmol_free, ntype);
        sys.ntype_with_fixed = ntype_with_fixed;
        sys.nfixedat = ntotat_fixed;
        // comptype is initialized with size ntype; resize to include fixed types
        sys.comptype = vec![true; ntype_with_fixed];

        // Fill nmols, natoms, idfirst for free types
        let mut cum_atoms = 0usize;
        let mut coor = Vec::new();
        let mut maxmove_per_type = vec![0usize; ntype];

        sys.nmols = vec![0; ntype_with_fixed];
        sys.natoms = vec![0; ntype_with_fixed];
        sys.idfirst = vec![0; ntype_with_fixed];
        sys.constrain_rot = vec![[false; 3]; ntype];
        sys.rot_bound = vec![[[0.0; 2]; 3]; ntype];

        for (itype, target) in free_targets.iter().enumerate() {
            sys.nmols[itype] = target.count;
            sys.natoms[itype] = target.natoms();
            sys.idfirst[itype] = cum_atoms;
            coor.extend_from_slice(reference_coords(target));
            cum_atoms += target.natoms();

            maxmove_per_type[itype] = target.perturb_budget.unwrap_or(target.count);
            for k in 0..3 {
                if let Some((center, half_width)) = target.rotation_bound[k] {
                    sys.constrain_rot[itype][k] = true;
                    sys.rot_bound[itype][k][0] = center.radians();
                    sys.rot_bound[itype][k][1] = half_width.radians();
                }
            }
        }

        for (fi, target) in fixed_targets.iter().enumerate() {
            let itype = ntype + fi;
            sys.nmols[itype] = 1;
            sys.natoms[itype] = target.natoms();
            sys.idfirst[itype] = cum_atoms;
            coor.extend_from_slice(reference_coords(target));
            cum_atoms += target.natoms();
        }
        sys.coor = coor;

        // Assign radii, element symbols, and per-atom (itype, imol) tags.
        // Packmol uses `radius = tolerance/2` for ALL atoms (packmol.f90 line 283:
        //   `radius(i) = dism/2.d0`), not VdW radii from the PDB file.
        // `ibtype` / `ibmol` are derivable from position in the sequential
        // atom layout, so we set them here once instead of having
        // `insert_atom_in_cell` rewrite the same constants on every eval.
        let atom_radius = self.tolerance / 2.0;
        let mut icart = 0usize;
        for (itype, target) in free_targets.iter().enumerate() {
            for imol in 0..target.count {
                for iatom in 0..target.natoms() {
                    sys.radius[icart] = atom_radius;
                    sys.radius_ini[icart] = atom_radius;
                    sys.ibtype[icart] = itype;
                    sys.ibmol[icart] = imol;
                    sys.elements[icart] = Element::by_symbol(&target.elements[iatom]);
                    icart += 1;
                }
            }
        }
        for (fi, target) in fixed_targets.iter().enumerate() {
            let itype = ntype + fi;
            for iatom in 0..target.natoms() {
                sys.radius[icart] = atom_radius;
                sys.radius_ini[icart] = atom_radius;
                sys.ibtype[icart] = itype;
                sys.ibmol[icart] = 0;
                sys.elements[icart] = Element::by_symbol(&target.elements[iatom]);
                icart += 1;
            }
        }

        // Assign restraints: per-atom
        let mut irest_pool = Vec::new();
        let mut iratom_lists: Vec<Vec<usize>> = vec![Vec::new(); ntotat];
        let mut icart = 0usize;
        for target in free_targets.iter() {
            for _imol in 0..target.count {
                for iatom in 0..target.natoms() {
                    // molecule-level restraints applied to all atoms
                    for r in &target.molecule_restraints {
                        let irest = irest_pool.len();
                        irest_pool.push(std::sync::Arc::clone(r));
                        iratom_lists[icart].push(irest);
                    }
                    // atom-subset restraints
                    for (indices, restraint) in &target.atom_restraints {
                        if indices.contains(&iatom) {
                            let irest = irest_pool.len();
                            irest_pool.push(std::sync::Arc::clone(restraint));
                            iratom_lists[icart].push(irest);
                        }
                    }
                    icart += 1;
                }
            }
        }
        // Fixed atoms: no restraints needed (they are placed directly)
        sys.restraints = irest_pool;
        sys.iratom_offsets.clear();
        sys.iratom_offsets.reserve(ntotat + 1);
        sys.iratom_offsets.push(0);
        for atom_restraints in &iratom_lists {
            let next = sys.iratom_offsets.last().copied().unwrap_or(0) + atom_restraints.len();
            sys.iratom_offsets.push(next);
        }
        sys.iratom_data.clear();
        sys.iratom_data
            .reserve(sys.iratom_offsets.last().copied().unwrap_or(0));
        for atom_restraints in iratom_lists {
            sys.iratom_data.extend(atom_restraints);
        }

        // Handle fixed molecules: place them using eulerfixed
        let free_atoms = ntotat_free;
        let mut fixed_icart = free_atoms;
        for target in fixed_targets.iter() {
            let fp = target.fixed_at.as_ref().unwrap();
            let (v1, v2, v3) = eulerfixed(
                fp.orientation[0].radians(),
                fp.orientation[1].radians(),
                fp.orientation[2].radians(),
            );
            let ref_coords = reference_coords(target);
            for ref_coord in ref_coords.iter().take(target.natoms()) {
                let pos = compcart(&fp.position, ref_coord, &v1, &v2, &v3);
                sys.xcart[fixed_icart] = pos;
                sys.fixedatom[fixed_icart] = true;
                fixed_icart += 1;
            }
        }
        // Populate the AoS `atom_props` mirror from the individual per-atom
        // Vecs now that every hot-loop field is finalized. `sync_atom_props`
        // also refreshes the `any_fixed_atoms` / `any_short_radius`
        // summary flags used by the hot-loop fast paths.
        sys.sync_atom_props();
        // Plumb the builder's opt-in parallel flag through to the
        // objective kernels.
        sys.parallel_pair_eval = self.parallel_eval;

        // Write constant columns (element, mol_id) into the output frame.
        // These don't change during optimization; positions are added at the end.
        crate::frame::init_frame_constants(&mut sys);

        // Initialize x vector
        let mut x = vec![0.0 as F; n];

        // Notify handlers immediately (before any heavy computation)
        for h in self.handlers.iter_mut() {
            h.on_start(ntotat, ntotmol_free);
        }

        // Run initialization
        sys.ntotmol = ntotmol_free;
        let init_passes = self.init_passes.unwrap_or(20 * ntype);
        let movebad_cfg = MoveBadConfig {
            movefrac: self.perturb_fraction,
            maxmove_per_type: &maxmove_per_type,
            movebadrandom: self.random_perturb,
            gencan_maxit: self.inner_iterations,
        };
        initial(
            &mut x,
            &mut sys,
            self.precision,
            self.discale,
            self.init_box_half_size,
            init_passes,
            pbc,
            &movebad_cfg,
            &mut rng,
        );

        // Notify handlers: initialization complete, xcart is valid
        for h in self.handlers.iter_mut() {
            h.on_initialized(&sys);
        }

        // Build relaxer runners from target relaxers (RelaxerRunner carries mutable MC state).
        // Each entry: (type_index, Vec<Box<dyn RelaxerRunner>>).
        let mut relaxer_runners: Vec<(usize, Vec<Box<dyn RelaxerRunner>>)> = free_targets
            .iter()
            .enumerate()
            .filter(|(_, t)| !t.relaxers.is_empty())
            .map(|(i, t)| {
                let base = sys.idfirst[i];
                let na = sys.natoms[i];
                let ref_slice = &sys.coor[base..base + na];
                let runners = t.relaxers.iter().map(|r| r.spawn(ref_slice)).collect();
                (i, runners)
            })
            .collect();

        // max_loops controls the outer loop count, matching Packmol's `nloop` parameter.
        let gencan_params = GencanParams {
            maxit: self.inner_iterations,
            maxfc: self.inner_iterations * 10,
            iprint: 0,
            ..Default::default()
        };

        let mut converged = false;
        let mut gencan_workspace = GencanWorkspace::new();

        // ── Main optimization loop ─────────────────────────────────────────────
        //
        // Matches Packmol's `app/packmol.f90` main loop exactly:
        //   For each type (itype 1..ntype): swaptype(action=1) → pack → restore
        //   Then all types (itype = ntype+1): pack with full x
        //
        // Per-type phases use a compact x (n = nmols[itype]*6) via SwapState,
        // reducing GENCAN problem size by up to 60x vs full n.

        // Save initial full x before phasing (Packmol swaptype action=0 at line 740)
        let mut swap = SwapState::init(&x, &sys);

        let total_phases = ntype + 1;

        for phase in 0..=(ntype) {
            let outcome = run_phase(
                phase,
                ntype,
                ntype_with_fixed,
                total_phases,
                max_loops,
                self.discale,
                self.precision,
                !self.perturb,
                &movebad_cfg,
                &gencan_params,
                &mut sys,
                &mut x,
                &mut swap,
                &mut relaxer_runners,
                &mut self.handlers,
                &mut gencan_workspace,
                &mut rng,
            );
            match outcome {
                PhaseOutcome::Continue => {}
                PhaseOutcome::Converged => {
                    converged = true;
                    break;
                }
            }
        }

        if !converged {
            log::warn!(
                "  Pack did not fully converge (fdist={:.4e}, frest={:.4e})",
                sys.fdist,
                sys.frest
            );
        }

        // Rebuild final xcart from x (all types active)
        for itype in 0..ntype_with_fixed {
            sys.comptype[itype] = true;
        }
        sys.ntotmol = ntotmol_free;
        init_xcart_from_x(&x, &mut sys);

        // Notify handlers of final state
        for h in self.handlers.iter_mut() {
            h.on_finish(&sys);
        }

        // Finalize frame: write positions into frame, move it out (zero-copy)
        let frame = crate::frame::finalize_frame(&mut sys);
        Ok(PackResult {
            frame,
            fdist: sys.fdist,
            frest: sys.frest,
            converged,
        })
    }
}

/// Resolved periodic-box spec: `(min, max, periodic_flags)`. Shared by
/// `derive_periodic_box` and its callers.
type PeriodicSpec = ([F; 3], [F; 3], [bool; 3]);

/// Scan every restraint on every target for a `Restraint::periodic_box`
/// override. Returns `Ok(None)` if no restraint declares one, `Ok(Some(...))`
/// if exactly one unique declaration exists (duplicates with identical
/// bounds + flags are allowed — they come from `with_global_restraint`
/// broadcast and from two targets sharing the same restraint object).
/// Returns `Err(ConflictingPeriodicBoxes)` when two declarations disagree
/// and `Err(InvalidPBCBox)` if the declared box has a non-positive extent
/// on any axis.
fn derive_periodic_box(targets: &[Target]) -> Result<Option<PeriodicSpec>, PackError> {
    let mut found: Option<PeriodicSpec> = None;
    for target in targets {
        let restraints = target
            .molecule_restraints
            .iter()
            .chain(target.atom_restraints.iter().map(|(_, r)| r));
        for r in restraints {
            if let Some(candidate) = r.periodic_box() {
                let (min, max, _periodic) = candidate;
                let length = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
                if length.iter().any(|&v| v <= 0.0) {
                    return Err(PackError::InvalidPBCBox { min, max });
                }
                match found {
                    None => found = Some(candidate),
                    Some(existing) if existing == candidate => {}
                    Some(existing) => {
                        return Err(PackError::ConflictingPeriodicBoxes {
                            first: existing,
                            second: candidate,
                        });
                    }
                }
            }
        }
    }
    Ok(found)
}

fn reference_coords(target: &Target) -> &[[F; 3]] {
    match target.centering {
        CenteringMode::Center => &target.ref_coords,
        CenteringMode::Off => &target.input_coords,
        CenteringMode::Auto => {
            if target.fixed_at.is_some() {
                &target.input_coords
            } else {
                &target.ref_coords
            }
        }
    }
}

/// Evaluate the packing objective once under **unscaled** radii (`radius_ini`),
/// restoring the caller's `radius` values on return.
///
/// On return, `sys.fdist` / `sys.frest` / `sys.fdist_atom` / `sys.frest_atom`
/// reflect the unscaled evaluation (the radius-dependent inner state); only
/// `sys.radius` itself is rolled back to what it was on entry.
///
/// Returns `(f_total, fdist, frest)` from the unscaled evaluation — the exact
/// triple the packer's main loop feeds to `flast` / `fimp` / handler `StepInfo`.
///
/// Pulled out of `pack()` in phase A.4.1 to de-duplicate three inline copies
/// of the same swap-evaluate-restore dance (Packmol `computef` emulation).
pub fn evaluate_unscaled(sys: &mut PackContext, xwork: &[F]) -> (F, F, F) {
    sys.work.radiuswork.copy_from_slice(&sys.radius);
    for i in 0..sys.ntotat {
        sys.set_radius(i, sys.radius_ini[i]);
    }
    let f_total = sys.evaluate(xwork, EvalMode::FOnly, None).f_total;
    let fdist = sys.fdist;
    let frest = sys.frest;
    for i in 0..sys.ntotat {
        sys.set_radius(i, sys.work.radiuswork[i]);
    }
    (f_total, fdist, frest)
}

/// Outcome of one main-loop iteration inside a packing phase.
///
/// Pulled out of `pack()` in phase A.4.3 to isolate the ~140-line per-iteration
/// body that runs movebad → relaxers → pgencan → radii schedule. `Continue`
/// means "run the next iteration"; `Converged` means the convergence predicate
/// fired inside this iteration; `EarlyStop` means a `Handler::should_stop()`
/// returned true.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IterOutcome {
    Continue,
    Converged,
    EarlyStop,
}

/// Run one iteration of a packing phase's main loop.
///
/// Matches Packmol's per-iteration sequence in `app/packmol.f90` lines 815-948:
///
/// 1. `movebad` when `radscale == 1.0` and previous `fimp <= 10%` (unless
///    disabled).
/// 2. Per-target relaxer MC block.
/// 3. `pgencan` on the working coordinate vector.
/// 4. Unscaled-radii statistics (`fdist` / `frest` / `fimp`).
/// 5. Handler `on_step` notification; early stop if any handler opts in.
/// 6. Convergence check (`fdist < precision && frest < precision`).
/// 7. Radii reduction schedule (only when `radscale > 1.0`).
///
/// The function takes each piece of mutable outer-loop state by `&mut` so the
/// caller (the outer phase for-loop in `pack()`) retains ownership across
/// iterations.
#[allow(clippy::too_many_arguments)]
pub fn run_iteration(
    loop_idx: usize,
    max_loops: usize,
    is_all: bool,
    phase: usize,
    phase_info: PhaseInfo,
    precision: F,
    disable_movebad: bool,
    movebad_cfg: &MoveBadConfig,
    gencan_params: &GencanParams,
    sys: &mut PackContext,
    xwork: &mut [F],
    swap: &mut SwapState,
    flast: &mut F,
    fimp_prev: &mut F,
    radscale: &mut F,
    relaxer_runners: &mut Vec<(usize, Vec<Box<dyn RelaxerRunner>>)>,
    handlers: &mut [Box<dyn Handler>],
    gencan_workspace: &mut GencanWorkspace,
    rng: &mut SmallRng,
) -> IterOutcome {
    // movebad: Packmol triggers when radscale==1.0 AND fimp<=10.0
    // (packmol.f90 line 815). fimp here is from the PREVIOUS iteration.
    // After movebad, reset flast to the post-movebad f (Packmol line 821).
    if !disable_movebad && *radscale == 1.0 && *fimp_prev <= 10.0 {
        movebad(xwork, sys, precision, movebad_cfg, rng, gencan_workspace);
        // Reset flast to the post-movebad f value so fimp is measured
        // relative to movebad's starting point.
        *flast = evaluate_unscaled(sys, xwork).0;
    }

    // Relaxer MC block: run per-target relaxers between movebad and pgencan.
    // Each relaxer modifies the reference coords (coor) for its type.
    for (itype, runners) in relaxer_runners.iter_mut() {
        if !is_all && *itype != phase {
            continue;
        }

        let base = sys.idfirst[*itype];
        let na = sys.natoms[*itype];

        for runner in runners.iter_mut() {
            let saved: Vec<[F; 3]> = sys.coor[base..base + na].to_vec();
            let f_before = sys.evaluate(xwork, EvalMode::FOnly, None).f_total;

            let result = runner.on_iter(
                &saved,
                f_before,
                &mut |trial: &[[F; 3]]| {
                    sys.coor[base..base + na].copy_from_slice(trial);
                    let f = sys.evaluate(xwork, EvalMode::FOnly, None).f_total;
                    sys.coor[base..base + na].copy_from_slice(&saved);
                    f
                },
                rng,
            );

            if let Some(new_coords) = result {
                sys.coor[base..base + na].copy_from_slice(&new_coords);
            }
        }
    }

    // GENCAN on working x (compact for per-type, full for all-type)
    sys.reset_eval_counters();
    let res = pgencan(xwork, sys, gencan_params, precision, gencan_workspace);

    // Save compact results back to swap (for restore later)
    if !is_all {
        swap.save_type(phase, xwork, sys);
    }

    // Compute statistics with unscaled radii
    // (Packmol lines 833-841: radiuswork + computef + restore)
    let (fx_unscaled, fdist, frest) = evaluate_unscaled(sys, xwork);

    // fimp: percentage improvement in unscaled f from last iteration
    // Packmol line 846: if(flast>0) fimp = -100*(fx-flast)/flast
    let mut fimp = if *flast > 0.0 {
        -100.0 * (fx_unscaled - *flast) / *flast
    } else if fx_unscaled < objective_small_floor() {
        100.0 // already converged
    } else {
        F::INFINITY
    };
    // Packmol lines 848-849: clamp to [-99.99, 99.99]
    fimp = fimp.clamp(-99.99, 99.99);
    *flast = fx_unscaled;
    *fimp_prev = fimp;

    if !handlers.is_empty() {
        let relaxer_acceptance: Vec<(usize, F)> = relaxer_runners
            .iter()
            .flat_map(|(itype, runners)| runners.iter().map(move |r| (*itype, r.acceptance_rate())))
            .collect();

        let step_info = StepInfo {
            loop_idx,
            max_loops,
            phase: phase_info,
            fdist,
            frest,
            improvement_pct: fimp,
            radscale: *radscale,
            precision,
            relaxer_acceptance,
        };
        for h in handlers.iter_mut() {
            h.on_step(&step_info, sys);
        }

        if handlers.iter().any(|h| h.should_stop()) {
            log::debug!("  Early stop requested at loop {loop_idx}");
            return IterOutcome::EarlyStop;
        }
    }

    log::debug!(
        "    loop={loop_idx} f={:.4e} fdist={:.4e} frest={:.4e} radscale={:.4} fimp={:.2}% ncf={} ncg={} inform={}",
        res.f,
        fdist,
        frest,
        *radscale,
        fimp,
        sys.ncf(),
        sys.ncg(),
        res.inform
    );

    // Check convergence
    if fdist < precision && frest < precision {
        log::debug!("  Converged at phase {phase} loop {loop_idx}");
        return IterOutcome::Converged;
    }

    // Radii reduction schedule (Packmol lines 940-948):
    //   if (fdist<precision && fimp<10%) || fimp<2%: reduce radscale
    if *radscale > 1.0 && (fimp < 2.0 || (fdist < precision && fimp < 10.0)) {
        *radscale = (0.9 * *radscale).max(1.0);
        for i in 0..sys.ntotat {
            let new_r = sys.radius_ini[i].max(0.9 * sys.radius[i]);
            sys.set_radius(i, new_r);
        }
    }

    IterOutcome::Continue
}

/// Outcome of one outer-loop phase in `pack()`.
///
/// Pulled out of `pack()` in phase A.4.2 to isolate the outer per-phase scaffold
/// (handler phase-start notification, comptype reconfiguration, radii reset,
/// swap setup, pre-loop precision short-circuit, inner GENCAN loop, swap
/// restore / xwork-back copy). `Continue` means the outer phase loop should
/// proceed; `Converged` means the all-type phase converged and the outer loop
/// should break.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseOutcome {
    Continue,
    Converged,
}

/// Run one phase of the main packing loop (per-type or all-type).
///
/// Matches the outer `for phase in 0..=ntype` body of Packmol `app/packmol.f90`
/// lines 740-990 (the swaptype / comptype dance bracketing the GENCAN inner
/// loop). For a per-type phase (`phase < ntype`), `xwork` is a compact
/// `nmols[phase] * 6`-element slice produced by `SwapState::set_type`; for the
/// all-type phase (`phase == ntype`), `xwork` is a full `6 * ntotmol_free`
/// clone of `x`.
///
/// The function takes the outer-loop state (`sys`, `x`, `swap`,
/// `relaxer_runners`, `handlers`, `gencan_workspace`, `rng`) by `&mut` so that
/// state persists across phases, exactly as the inlined body did.
///
/// Returns `PhaseOutcome::Converged` **only** when the all-type phase
/// converges (either on its entry precision check or inside the inner loop);
/// every per-type phase returns `Continue` regardless of whether that type
/// converged on its own (Packmol lets the all-type phase decide).
#[allow(clippy::too_many_arguments)]
pub fn run_phase(
    phase: usize,
    ntype: usize,
    ntype_with_fixed: usize,
    total_phases: usize,
    max_loops: usize,
    discale: F,
    precision: F,
    disable_movebad: bool,
    movebad_cfg: &MoveBadConfig,
    gencan_params: &GencanParams,
    sys: &mut PackContext,
    x: &mut Vec<F>,
    swap: &mut SwapState,
    relaxer_runners: &mut Vec<(usize, Vec<Box<dyn RelaxerRunner>>)>,
    handlers: &mut [Box<dyn Handler>],
    gencan_workspace: &mut GencanWorkspace,
    rng: &mut SmallRng,
) -> PhaseOutcome {
    let is_all = phase == ntype;

    let phase_info = PhaseInfo {
        phase,
        total_phases,
        molecule_type: if is_all { None } else { Some(phase) },
    };

    // Reset handler state between phases (e.g. EarlyStopHandler stall counter)
    for h in handlers.iter_mut() {
        h.on_phase_start(&phase_info);
    }

    // Set comptype for this phase
    for itype in 0..ntype_with_fixed {
        sys.comptype[itype] = if is_all {
            true
        } else {
            itype >= ntype || itype == phase
        };
    }

    log::debug!(
        "  Packing phase {phase} ({})",
        if is_all {
            "all".to_string()
        } else {
            format!("type {phase}")
        }
    );

    // Compact x to this type (action=1) or restore full x (all-type phase)
    // Packmol resets radscale = discale at the START of each phase.
    let mut radscale = discale;
    for icart in 0..sys.ntotat {
        sys.set_radius(icart, discale * sys.radius_ini[icart]);
    }

    // Get working x vector (compact for per-type, full for all-type)
    let mut xwork: Vec<F> = if !is_all {
        // Compact: n = nmols[phase] * 6
        // Re-save current x (action=0) then compact (action=1)
        *swap = SwapState::init(x, sys);
        swap.set_type(phase, sys)
    } else {
        // All-type: restore full x (action=3), use it directly
        swap.restore(x, sys);
        x.clone()
    };

    // Packmol checks whether the current approximation is already a solution
    // before entering the GENCAN loop for this phase (packmol.f90 lines 775-782).
    sys.evaluate(&xwork, EvalMode::FOnly, None);
    if sys.fdist < precision && sys.frest < precision {
        if !is_all {
            swap.save_type(phase, &xwork, sys);
            swap.restore(x, sys);
            return PhaseOutcome::Continue;
        } else {
            x.clone_from(&xwork);
            return PhaseOutcome::Converged;
        }
    }

    // Initialize flast = unscaled f before gencanloop
    // (Packmol lines 796-803: compute bestf/flast with unscaled radii)
    let mut flast = evaluate_unscaled(sys, &xwork).0;

    // fimp from previous iteration — used for movebad gate (Packmol packmol.f90 line 798).
    // Initialized to 1e99 so movebad is NOT called on the first iteration.
    let mut fimp_prev = F::INFINITY;
    let mut converged_inner = false;

    for loop_idx in 0..max_loops {
        let outcome = run_iteration(
            loop_idx,
            max_loops,
            is_all,
            phase,
            phase_info,
            precision,
            disable_movebad,
            movebad_cfg,
            gencan_params,
            sys,
            &mut xwork,
            swap,
            &mut flast,
            &mut fimp_prev,
            &mut radscale,
            relaxer_runners,
            handlers,
            gencan_workspace,
            rng,
        );
        match outcome {
            IterOutcome::Continue => {}
            IterOutcome::Converged => {
                converged_inner = true;
                break;
            }
            IterOutcome::EarlyStop => break,
        }
    }

    // After per-type phase: save results + restore full x
    // After all-type phase: copy xwork back to x
    if !is_all {
        // save_type was called inside the loop; restore full x now.
        // Per-type convergence does NOT exit the outer phase loop.
        swap.restore(x, sys);
        PhaseOutcome::Continue
    } else {
        x.clone_from(&xwork);
        if converged_inner {
            PhaseOutcome::Converged
        } else {
            PhaseOutcome::Continue
        }
    }
}
