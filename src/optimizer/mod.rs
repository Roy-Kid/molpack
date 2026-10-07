//! In-loop geometry optimizers driven by [`molrs::optimize::Optimizer`].
//!
//! Callers construct a molrs optimizer (`Lbfgs` over a force-field potential,
//! e.g. `Lbfgs::new(Arc::new(SoftSpec::from_frame(..).potential(None)), ..)`
//! for the soft overlap + 1-2 / 1-3 objective — both live in molrs's `ff`
//! module, which the caller enables on its own molrs dependency — or molpack's
//! [`TorsionMcOptimizer`]) and bind it with
//! [`GencanPack::with_optimizer`](crate::GencanPack::with_optimizer) plus an
//! [`OptimizeSelect`] that names which components to assemble each call.

use molrs::core::Frame;
use molrs::core::Mic;
use molrs::op::F;
use molrs::optimize::{Optimizer, set_free_mask};

use crate::Objective;
use crate::context::PackContext;
use crate::euler::eulerrmat;
use crate::eval::EvalMode;
use crate::target::centered_coords;

mod torsion_mc;
pub use torsion_mc::TorsionMcOptimizer;

/// How selected components are optimized.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptimizeMode {
    /// Each copy of each named target, independently.
    PerCopy,
    /// All copies of the named targets as one joint movable group.
    Joint,
}

/// Which components to assemble into the Frame passed to [`Optimizer::minimize`].
#[derive(Debug, Clone)]
pub struct OptimizeSelect {
    pub names: Vec<String>,
    pub mode: OptimizeMode,
    /// Include frozen neighbour atoms (`atoms.free = false`) within `rcut`.
    pub with_environment: bool,
    pub rcut: F,
}

impl OptimizeSelect {
    pub fn per_copy(names: impl IntoIterator<Item = impl Into<String>>) -> Self {
        Self {
            names: names.into_iter().map(Into::into).collect(),
            mode: OptimizeMode::PerCopy,
            with_environment: false,
            rcut: 5.0,
        }
    }

    pub fn joint(names: impl IntoIterator<Item = impl Into<String>>) -> Self {
        Self {
            names: names.into_iter().map(Into::into).collect(),
            mode: OptimizeMode::Joint,
            with_environment: false,
            rcut: 5.0,
        }
    }

    pub fn with_environment(mut self, rcut: F) -> Self {
        self.with_environment = true;
        self.rcut = rcut;
        self
    }
}

/// One bound optimizer + selection, stored on [`crate::GencanPack`].
pub(crate) struct OptimizerBinding {
    pub select: OptimizeSelect,
    pub optimizer: Box<dyn Optimizer>,
}

/// One [`OptimizerBinding`] **borrowed** for the duration of a single run,
/// with its target names resolved to type indices.
///
/// It appears in the signatures of the GENCAN phase helpers
/// `run_iteration` / `run_phase`.
/// Built by `resolve_bindings` at the start of every run — not constructed
/// by callers.
///
/// It borrows rather than owns because of the stage seam's re-entrancy
/// contract ([`Stage::run`](crate::Stage::run)): the bindings are the
/// stage's own configuration and must still be there on the next run, so a
/// run may resolve them but never take them. Cloning is not the alternative
/// — [`Optimizer`] is a trait object with no `Clone` bound.
pub(crate) struct ResolvedBinding<'a> {
    pub select: &'a OptimizeSelect,
    pub type_indices: Vec<usize>,
    pub optimizer: &'a mut dyn Optimizer,
}

/// Resolve `bindings` against this run's `type_names`, borrowing each one.
///
/// Type indices are recomputed per run rather than cached: it is a name
/// comparison over a handful of targets, and a stage handed a different
/// target set gets the right answer for free.
pub(crate) fn resolve_bindings<'a>(
    bindings: &'a mut [OptimizerBinding],
    type_names: &[Option<String>],
) -> Vec<ResolvedBinding<'a>> {
    let mut out = Vec::with_capacity(bindings.len());
    for binding in bindings.iter_mut() {
        // Split the borrow: `select` is read while `optimizer` is driven.
        let OptimizerBinding { select, optimizer } = binding;
        let mut idxs = Vec::new();
        for name in &select.names {
            let mut found = false;
            for (i, tn) in type_names.iter().enumerate() {
                if tn.as_deref() == Some(name.as_str()) {
                    idxs.push(i);
                    found = true;
                }
            }
            if !found {
                log::warn!("with_optimizer: no target named '{name}', skipping");
            }
        }
        if idxs.is_empty() {
            continue;
        }
        out.push(ResolvedBinding {
            select,
            type_indices: idxs,
            optimizer: &mut **optimizer,
        });
    }
    out
}

/// Run all bound optimizers (all-type phase only for clean COM/Euler indexing).
pub(crate) fn run_optimizer_bindings(
    sys: &mut PackContext,
    xwork: &[F],
    bindings: &mut [ResolvedBinding<'_>],
) {
    if bindings.is_empty() {
        return;
    }
    let xcart_snapshot = sys.xcart.clone();
    let mic = sys.simbox.mic().simplified();

    for binding in bindings.iter_mut() {
        match binding.select.mode {
            OptimizeMode::PerCopy => {
                for &itype in &binding.type_indices {
                    let na = sys.natoms[itype];
                    let mol_offset: usize = sys.nmols[..itype].iter().sum();
                    for imol in 0..sys.nmols[itype] {
                        let start = sys.idfirst[itype] + imol * na;
                        let spans = [CopySpan {
                            copy_start: start,
                            na,
                            group_offset: 0,
                            ilugan: sys.ntotmol * 3 + (mol_offset + imol) * 3,
                        }];
                        let world: Vec<[F; 3]> = xcart_snapshot[start..start + na].to_vec();
                        optimize_group(
                            sys,
                            xwork,
                            &xcart_snapshot,
                            &mic,
                            &spans,
                            &world,
                            binding.select,
                            &mut *binding.optimizer,
                        );
                    }
                }
            }
            OptimizeMode::Joint => {
                let mut spans = Vec::new();
                let mut world = Vec::new();
                for &itype in &binding.type_indices {
                    let na = sys.natoms[itype];
                    let mol_offset: usize = sys.nmols[..itype].iter().sum();
                    for imol in 0..sys.nmols[itype] {
                        let copy_start = sys.idfirst[itype] + imol * na;
                        let group_offset = world.len();
                        for a in 0..na {
                            world.push(xcart_snapshot[copy_start + a]);
                        }
                        spans.push(CopySpan {
                            copy_start,
                            na,
                            group_offset,
                            ilugan: sys.ntotmol * 3 + (mol_offset + imol) * 3,
                        });
                    }
                }
                if world.is_empty() {
                    continue;
                }
                optimize_group(
                    sys,
                    xwork,
                    &xcart_snapshot,
                    &mic,
                    &spans,
                    &world,
                    binding.select,
                    &mut *binding.optimizer,
                );
            }
        }
    }
}

struct CopySpan {
    copy_start: usize,
    na: usize,
    group_offset: usize,
    ilugan: usize,
}

#[allow(clippy::too_many_arguments)]
fn optimize_group(
    sys: &mut PackContext,
    xwork: &[F],
    xcart_snapshot: &[[F; 3]],
    mic: &Mic,
    spans: &[CopySpan],
    world_movable: &[[F; 3]],
    select: &OptimizeSelect,
    optimizer: &mut dyn Optimizer,
) {
    // Build free mask + optional environment.
    let mut world = world_movable.to_vec();
    let mut free = vec![true; world.len()];
    if select.with_environment {
        let mut movable_mask = vec![false; xcart_snapshot.len()];
        for s in spans {
            for a in 0..s.na {
                movable_mask[s.copy_start + a] = true;
            }
        }
        for icart in environment_atoms(
            xcart_snapshot,
            &movable_mask,
            world_movable,
            select.rcut,
            mic,
        ) {
            world.push(xcart_snapshot[icart]);
            free.push(false);
        }
    }

    let mut frame = coords_to_frame(&world);
    if free.iter().any(|&f| !f) {
        let _ = set_free_mask(&mut frame, &free);
    }

    if optimizer.minimize(&mut frame).is_err() {
        return;
    }

    // Read back free (movable) coordinates only.
    let Ok(xyz) = frame.coords() else {
        return;
    };
    let mut world_new = crate::template::coord_rows(&xyz);
    world_new.truncate(world_movable.len());

    // Map world delta → reference coor per copy (Rᵀ + recenter), then non-harm gate.
    //
    // Both evaluations must run on a fresh Cartesian expansion. The cache is
    // keyed on `x`, which this function never touches — it rewrites `coor` —
    // so without invalidation the second call returns the first call's value
    // and the gate below can never fire. `f_before` needs it too: an earlier
    // copy in this same sweep may already have rewritten its own conformer.
    sys.invalidate_geometry_cache();
    let f_before = sys.evaluate(xwork, EvalMode::FOnly, None).f_total;
    let mut saved: Vec<(usize, usize, Vec<[F; 3]>)> = Vec::with_capacity(spans.len());
    for s in spans {
        let (v1, v2, v3) = eulerrmat(xwork[s.ilugan], xwork[s.ilugan + 1], xwork[s.ilugan + 2]);
        let coor_old: Vec<[F; 3]> = sys.coor[s.copy_start..s.copy_start + s.na].to_vec();
        let mut coor_new = vec![[0.0 as F; 3]; s.na];
        for a in 0..s.na {
            let g = s.group_offset + a;
            let d = [
                world_new[g][0] - world_movable[g][0],
                world_new[g][1] - world_movable[g][1],
                world_new[g][2] - world_movable[g][2],
            ];
            let rt = [
                v1[0] * d[0] + v1[1] * d[1] + v1[2] * d[2],
                v2[0] * d[0] + v2[1] * d[1] + v2[2] * d[2],
                v3[0] * d[0] + v3[1] * d[1] + v3[2] * d[2],
            ];
            coor_new[a] = [
                coor_old[a][0] + rt[0],
                coor_old[a][1] + rt[1],
                coor_old[a][2] + rt[2],
            ];
        }
        let coor_new = centered_coords(&coor_new);
        saved.push((s.copy_start, s.na, coor_old));
        sys.coor[s.copy_start..s.copy_start + s.na].copy_from_slice(&coor_new);
    }
    sys.invalidate_geometry_cache();
    let f_after = sys.evaluate(xwork, EvalMode::FOnly, None).f_total;
    if f_after > f_before {
        for (cs, na, coor_old) in &saved {
            sys.coor[*cs..*cs + *na].copy_from_slice(coor_old);
        }
        sys.invalidate_geometry_cache();
    }
}

fn coords_to_frame(world: &[[F; 3]]) -> Frame {
    let mut frame = Frame::new();
    frame
        .set_coords(ndarray::Array2::from(world.to_vec()).view())
        .expect("an N x 3 array always fits a fresh frame");
    frame
}

/// Indices of the non-movable atoms within `rcut` of any movable atom, under
/// the cell's minimum image — the frozen environment an optimizer relaxes against.
fn environment_atoms(
    xcart: &[[F; 3]],
    movable_mask: &[bool],
    world_movable: &[[F; 3]],
    rcut: F,
    mic: &Mic,
) -> Vec<usize> {
    let rcut2 = rcut * rcut;
    (0..xcart.len())
        .filter(|&i| !movable_mask[i])
        .filter(|&i| {
            let p = xcart[i];
            world_movable.iter().any(|w| {
                let d = mic.apply([p[0] - w[0], p[1] - w[1], p[2] - w[2]]);
                d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < rcut2
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use molrs::core::SimBox;
    use ndarray::array;

    use super::*;

    /// A slab (`z` not periodic): an atom across the `z` face is far, one
    /// across the `x` face is near. Wrapping every axis once any axis is
    /// periodic pulled the `z` atom in as environment.
    #[test]
    fn environment_wraps_only_periodic_axes() {
        let bx = SimBox::ortho(
            array![10.0, 10.0, 10.0],
            array![0.0, 0.0, 0.0],
            [true, true, false],
        )
        .unwrap();
        let mic = bx.mic().simplified();
        let xcart = [[0.5, 5.0, 0.5], [9.5, 5.0, 0.5], [0.5, 5.0, 9.5]];
        let movable = [true, false, false];
        let near = environment_atoms(&xcart, &movable, &xcart[..1], 2.0, &mic);
        assert_eq!(near, vec![1]);
    }

    /// No periodic axis: nothing wraps, so neither face neighbour is near.
    #[test]
    fn environment_without_pbc_does_not_wrap() {
        let xcart = [[0.5, 5.0, 0.5], [9.5, 5.0, 0.5], [0.5, 5.0, 9.5]];
        let movable = [true, false, false];
        let near = environment_atoms(&xcart, &movable, &xcart[..1], 2.0, &Mic::Free);
        assert!(near.is_empty());
    }

    /// Records the atom count of every frame it is handed and leaves the
    /// frame alone, so the non-harm gate never fires and the sizes are the
    /// one thing a test observes.
    struct FrameSizes(std::sync::Arc<std::sync::Mutex<Vec<usize>>>);

    impl Optimizer for FrameSizes {
        fn minimize(
            &mut self,
            frame: &mut Frame,
        ) -> Result<molrs::optimize::OptimizationReport, String> {
            let n = frame.coords().map_err(|e| e.to_string())?.nrows();
            self.0.lock().unwrap().push(n);
            Ok(molrs::optimize::OptimizationReport {
                converged: true,
                n_steps: 0,
                final_energy: 0.0,
                final_fmax: 0.0,
                final_grad_rms: 0.0,
            })
        }
    }

    /// Pack 8 monomers `a` and 8 dimers `b` into a 4 Å cube — unsatisfiable,
    /// so the all-type phase iterates and reaches the optimizer block — with
    /// one optimizer bound to `select`. Returns the atom count of every frame
    /// the optimizer saw.
    fn frame_sizes(select: OptimizeSelect) -> Vec<usize> {
        use crate::{GencanPack, PackEngine, Target};

        let cube = || crate::test_fixtures::inside_box([0.0; 3], [4.0; 3]);
        let a = Target::from_coords(&[[0.0; 3]], &[1.0], 8)
            .with_name("a")
            .with_restraint(cube());
        let b = Target::from_coords(&[[0.0; 3], [1.5, 0.0, 0.0]], &[1.0, 1.0], 8)
            .with_name("b")
            .with_restraint(cube());
        let sizes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        GencanPack::new()
            .with_tolerance(2.0)
            .with_seed(7)
            .with_optimizer(select, FrameSizes(std::sync::Arc::clone(&sizes)))
            .run(&[a, b], 4)
            .expect("an unsatisfiable pack still returns a state");
        sizes.lock().unwrap().clone()
    }

    /// `PerCopy` hands the optimizer one copy at a time.
    #[test]
    fn per_copy_optimizes_one_copy_per_call() {
        let sizes = frame_sizes(OptimizeSelect::per_copy(["b"]));
        assert!(!sizes.is_empty(), "fixture guard: the optimizer never ran");
        assert!(sizes.iter().all(|&n| n == 2), "{sizes:?}");
    }

    /// `Joint` hands every copy of every named target over as one movable
    /// group: 8 monomers + 8 dimers = 24 atoms per call.
    #[test]
    fn joint_optimizes_every_named_copy_together() {
        let sizes = frame_sizes(OptimizeSelect::joint(["a", "b"]));
        assert!(!sizes.is_empty(), "fixture guard: the optimizer never ran");
        assert!(sizes.iter().all(|&n| n == 24), "{sizes:?}");
    }

    /// With the environment on, the frame also carries the frozen neighbours
    /// of the joint group — here the copies of `a`, which is not selected.
    #[test]
    fn joint_environment_adds_frozen_neighbours() {
        let sizes = frame_sizes(OptimizeSelect::joint(["b"]).with_environment(8.0));
        assert!(!sizes.is_empty(), "fixture guard: the optimizer never ran");
        // Every monomer sits in the 4 Å cube, within 8 Å of the dimers.
        assert!(sizes.iter().all(|&n| n == 16 + 8), "{sizes:?}");
    }

    /// molrs's soft-overlap L-BFGS bound jointly over two species, in a box
    /// crowded enough that the all-type phase iterates: it is called, the
    /// non-harm gate keeps it from worsening the pack, and the pack converges.
    #[test]
    fn joint_soft_lbfgs_over_two_species_converges() {
        use molrs::ff::potential::soft::SoftSpec;
        use molrs::optimize::{Lbfgs, LbfgsSettings};

        use crate::{GencanPack, PackEngine, Target};

        /// Counts the calls it forwards to the wrapped optimizer.
        struct Counted(Lbfgs, std::sync::Arc<std::sync::atomic::AtomicUsize>);
        impl Optimizer for Counted {
            fn minimize(
                &mut self,
                frame: &mut Frame,
            ) -> Result<molrs::optimize::OptimizationReport, String> {
                self.1.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                self.0.minimize(frame)
            }
        }

        let cube = || crate::test_fixtures::inside_box([0.0; 3], [10.0; 3]);
        let ion = Target::from_coords(&[[0.0; 3]], &[1.0], 30)
            .with_name("ion")
            .with_restraint(cube());
        let water = Target::from_coords(&[[0.0; 3], [2.0, 0.0, 0.0]], &[1.0, 1.0], 30)
            .with_name("water")
            .with_restraint(cube());
        // A coordinates-only group has no bonds: a pure soft-overlap push.
        let soft = Lbfgs::new(
            std::sync::Arc::new(SoftSpec::from_frame(&Frame::new()).potential(None)),
            LbfgsSettings {
                max_steps: 200,
                ..LbfgsSettings::DEFAULT
            },
        );
        let calls = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let state = GencanPack::new()
            .with_tolerance(2.0)
            .with_precision(1e-2)
            .with_seed(7)
            .with_optimizer(
                OptimizeSelect::joint(["ion", "water"]),
                Counted(soft, std::sync::Arc::clone(&calls)),
            )
            .run(&[ion, water], 400)
            .expect("pack runs");
        let calls = calls.load(std::sync::atomic::Ordering::Relaxed);
        assert!(calls > 0, "fixture guard: the optimizer never ran");
        assert!(
            state.converged,
            "calls={calls} fdist={} frest={}",
            state.fdist, state.frest
        );
    }

    /// A name that matches no target is skipped, not an error: the pack runs
    /// and the optimizer is never called.
    #[test]
    fn an_unmatched_name_binds_nothing() {
        let sizes = frame_sizes(OptimizeSelect::joint(["no_such_target"]));
        assert!(sizes.is_empty(), "{sizes:?}");
    }
}
