//! Single-molecule constraint solve.
//! Port of `restmol.f90`.
//!
//! Placement and the bad-move heuristic both call this. It calls neither
//! of them.

use molrs::op::types::F;

use crate::Objective;
use crate::context::PackContext;
use crate::eval::EvalMode;
use crate::pack::gencan::{GencanParams, GencanWorkspace, pgencan};

/// Scoped state override for `restmol`; restores context on drop.
struct RestmolScope<'a> {
    sys: &'a mut PackContext,
    itype: usize,
    ntotmol: usize,
    nmols_itype: usize,
    comptype: Vec<bool>,
    init1: bool,
}

impl<'a> RestmolScope<'a> {
    fn enter(sys: &'a mut PackContext, itype: usize) -> Self {
        let saved = Self {
            ntotmol: sys.ntotmol,
            nmols_itype: sys.nmols[itype],
            comptype: sys.comptype.clone(),
            init1: sys.init1,
            itype,
            sys,
        };

        saved.sys.ntotmol = 1;
        // Only reduce the active type to 1 molecule.
        // Other types keep their original nmols so compute_f's icart counter advances
        // correctly past them — preserving the constraint array index alignment.
        // (Packmol restmol.f90 line 34: only nmols(itype) = 1, others unchanged.)
        saved.sys.nmols[itype] = 1;
        for i in 0..saved.sys.ntype_with_fixed {
            saved.sys.comptype[i] = i == itype;
        }
        saved.sys.init1 = true; // constraint-only, no cell list

        saved
    }

    fn ctx_mut(&mut self) -> &mut PackContext {
        self.sys
    }
}

impl Drop for RestmolScope<'_> {
    fn drop(&mut self) {
        self.sys.ntotmol = self.ntotmol;
        self.sys.nmols[self.itype] = self.nmols_itype;
        self.sys.comptype.clone_from(&self.comptype);
        self.sys.init1 = self.init1;
    }
}

/// Run a single-molecule GENCAN solve (restmol).
/// Port of `restmol.f90`.
///
/// `ilubar` is the offset in `x` for the COM of this molecule.
/// Euler angles are at `x[ilubar + ntotmol*3 ..]`.
///
/// - `solve = false`: evaluate constraint function only (no optimization).
/// - `solve = true`: run GENCAN to minimize constraint violations.
///
/// On return, `sys.frest` holds the constraint violation for this molecule.
#[allow(clippy::too_many_arguments)]
pub fn restmol(
    itype: usize,
    ilubar: usize,
    x: &mut [F],
    sys: &mut PackContext,
    precision: F,
    gencan_maxit: usize,
    solve: bool,
    workspace: &mut GencanWorkspace,
) {
    let ilugan_offset = sys.ntotmol * 3;
    let mut xmol = vec![0.0 as F; 6];
    xmol[0] = x[ilubar];
    xmol[1] = x[ilubar + 1];
    xmol[2] = x[ilubar + 2];
    xmol[3] = x[ilubar + ilugan_offset];
    xmol[4] = x[ilubar + ilugan_offset + 1];
    xmol[5] = x[ilubar + ilugan_offset + 2];

    {
        let mut scope = RestmolScope::enter(sys, itype);
        let sys = scope.ctx_mut();
        if !solve {
            sys.evaluate(&xmol, EvalMode::FOnly, None);
        } else {
            let params = GencanParams {
                maxit: gencan_maxit,
                maxfc: gencan_maxit * 10,
                ..Default::default()
            };
            pgencan(&mut xmol, sys, &params, precision, workspace);
        }
    }

    x[ilubar] = xmol[0];
    x[ilubar + 1] = xmol[1];
    x[ilubar + 2] = xmol[2];
    x[ilubar + ilugan_offset] = xmol[3];
    x[ilubar + ilugan_offset + 1] = xmol[4];
    x[ilubar + ilugan_offset + 2] = xmol[5];
    // sys.frest retains the value from the restmol compute_f
}
