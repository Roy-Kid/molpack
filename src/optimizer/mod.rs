//! In-loop geometry optimizers driven by [`molrs::optimize::Optimizer`].
//!
//! Callers construct a molrs optimizer (e.g. `LBFGS`, `SoftSpec::into_optimizer`,
//! or molpack's [`TorsionMcOptimizer`]) and bind it with
//! [`Molpack::with_optimizer`](crate::Molpack::with_optimizer) plus an
//! [`OptimizeSelect`] that names which components to assemble each call.

#![cfg(feature = "ff")]

use molrs::optimize::{Optimizer, set_free_mask};
use molrs::store::block::Block;
use molrs::store::frame::Frame;
use molrs::types::F;
use ndarray::Array1;

use crate::context::PackContext;
use crate::constraints::EvalMode;
use crate::Objective;
use crate::euler::eulerrmat;


pub mod torsion_mc;
pub use torsion_mc::TorsionMcOptimizer;

/// How selected components are optimized.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptimizeMode {
    /// Each copy of each named target, independently.
    PerCopy,
    /// All copies of the named targets as one joint movable group.
    Joint,
}

/// Which components to assemble into the Frame passed to [`Optimizer::run`].
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

    pub fn without_environment(mut self) -> Self {
        self.with_environment = false;
        self
    }
}

/// One bound optimizer + selection, stored on [`crate::Molpack`].
pub struct OptimizerBinding {
    pub select: OptimizeSelect,
    pub optimizer: Box<dyn Optimizer>,
}

/// Resolved type indices for a binding (matched by target name at pack start).
pub(crate) struct ResolvedBinding {
    pub select: OptimizeSelect,
    pub type_indices: Vec<usize>,
    pub optimizer: Box<dyn Optimizer>,
}

pub(crate) fn resolve_bindings(
    bindings: Vec<OptimizerBinding>,
    type_names: &[Option<String>],
) -> Vec<ResolvedBinding> {
    let mut out = Vec::with_capacity(bindings.len());
    for b in bindings {
        let mut idxs = Vec::new();
        for name in &b.select.names {
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
            select: b.select,
            type_indices: idxs,
            optimizer: b.optimizer,
        });
    }
    out
}

/// Run all bound optimizers (all-type phase only for clean COM/Euler indexing).
pub(crate) fn run_optimizer_bindings(
    sys: &mut PackContext,
    xwork: &[F],
    bindings: &mut [ResolvedBinding],
) {
    if bindings.is_empty() {
        return;
    }
    let xcart_snapshot = sys.xcart.clone();
    let pbc = optimizer_pbc(sys);

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
                            pbc,
                            &spans,
                            &world,
                            &binding.select,
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
                    pbc,
                    &spans,
                    &world,
                    &binding.select,
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
    pbc: Option<[F; 3]>,
    spans: &[CopySpan],
    world_movable: &[[F; 3]],
    select: &OptimizeSelect,
    optimizer: &mut dyn Optimizer,
) {
    // Build free mask + optional environment.
    let mut world = world_movable.to_vec();
    let mut free = vec![true; world.len()];
    if select.with_environment {
        let rcut2 = select.rcut * select.rcut;
        let mut movable_mask = vec![false; xcart_snapshot.len()];
        for s in spans {
            for a in 0..s.na {
                movable_mask[s.copy_start + a] = true;
            }
        }
        for (icart, p) in xcart_snapshot.iter().enumerate() {
            if movable_mask[icart] {
                continue;
            }
            let near = world_movable.iter().any(|w| {
                let d = min_image([p[0] - w[0], p[1] - w[1], p[2] - w[2]], pbc);
                d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < rcut2
            });
            if near {
                world.push(*p);
                free.push(false);
            }
        }
    }

    let mut frame = coords_to_frame(&world);
    if free.iter().any(|&f| !f) {
        let _ = set_free_mask(&mut frame, &free);
    }

    if optimizer.run(&mut frame).is_err() {
        return;
    }

    // Read back free (movable) coordinates only.
    let new_flat = match molrs::ff::potential::extract_coords(&frame) {
        Ok(c) => c,
        Err(_) => return,
    };
    let n_movable = world_movable.len();
    let mut world_new: Vec<[F; 3]> = (0..n_movable)
        .map(|i| {
            [
                new_flat[3 * i],
                new_flat[3 * i + 1],
                new_flat[3 * i + 2],
            ]
        })
        .collect();

    // Map world delta → reference coor per copy (Rᵀ + recenter), then non-harm gate.
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
        recenter(&mut coor_new);
        saved.push((s.copy_start, s.na, coor_old));
        sys.coor[s.copy_start..s.copy_start + s.na].copy_from_slice(&coor_new);
    }
    let f_after = sys.evaluate(xwork, EvalMode::FOnly, None).f_total;
    if f_after > f_before {
        for (cs, na, coor_old) in &saved {
            sys.coor[*cs..*cs + *na].copy_from_slice(coor_old);
        }
    }
    let _ = world_new; // silence if unused after refactor
}

fn coords_to_frame(world: &[[F; 3]]) -> Frame {
    let n = world.len();
    let mut atoms = Block::new();
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    let mut z = Vec::with_capacity(n);
    for p in world {
        x.push(p[0]);
        y.push(p[1]);
        z.push(p[2]);
    }
    let _ = atoms.insert("x", Array1::from_vec(x).into_dyn());
    let _ = atoms.insert("y", Array1::from_vec(y).into_dyn());
    let _ = atoms.insert("z", Array1::from_vec(z).into_dyn());
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    frame
}

fn recenter(coords: &mut [[F; 3]]) {
    if coords.is_empty() {
        return;
    }
    let n = coords.len() as F;
    let mut c = [0.0 as F; 3];
    for p in coords.iter() {
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    c[0] /= n;
    c[1] /= n;
    c[2] /= n;
    for p in coords.iter_mut() {
        p[0] -= c[0];
        p[1] -= c[1];
        p[2] -= c[2];
    }
}

fn optimizer_pbc(sys: &PackContext) -> Option<[F; 3]> {
    let pbc = sys.pbc_periodic();
    if pbc.iter().any(|&p| p) {
        let l = sys.simbox.lengths();
        Some([l[0], l[1], l[2]])
    } else {
        None
    }
}

#[inline]
fn min_image(d: [F; 3], pbc: Option<[F; 3]>) -> [F; 3] {
    match pbc {
        Some(l) => std::array::from_fn(|k| {
            if l[k] > 0.0 {
                d[k] - (d[k] / l[k]).round() * l[k]
            } else {
                d[k]
            }
        }),
        None => d,
    }
}
