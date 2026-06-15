//! Python wrappers for molpack's in-loop relaxers (relaxation-assisted packing).
//!
//! Two relaxers attach to a `Target` via :meth:`Target.with_relaxer` and run
//! *during* the pack loop, between GENCAN iterations:
//!
//! * [`PyTorsionMcRelaxer`] — Monte-Carlo torsion-angle sampling
//!   ([`molpack::TorsionMcRelaxer`]). Pure molpack, **no external engine and no
//!   force field**: it folds one flexible chain by proposing random rotations
//!   about its rotatable bonds and accepting them against the packer objective
//!   (Metropolis). Built straight from the molecule's **frame topology** — the
//!   bond graph crosses the boundary in the same zero-copy frame the `Target`
//!   uses, promoted to a molrs [`Atomistic`] so rotatable bonds can be detected.
//! * [`PyLBFGSRelaxer`] — force-field L-BFGS geometry relaxer
//!   ([`molpack::LBFGSRelaxer`], `ff` feature). Constructed from a
//!   `molrs.ForceField` / `molpy.ForceField` (zero-copy FFI capsule, see
//!   [`crate::interop`]); the potential is compiled lazily on first use.

use pyo3::prelude::*;
use pyo3::types::PyAny;

// ── TorsionMcRelaxer (core; needs no `ff` feature) ──────────────────────────

use molpack::TorsionMcRelaxer;

/// Monte-Carlo torsion-angle relaxer, attachable to a `Target` via
/// :meth:`Target.with_relaxer`. Engine-free and force-field-free.
#[pyclass(name = "TorsionMcRelaxer", from_py_object)]
#[derive(Clone)]
pub struct PyTorsionMcRelaxer {
    pub(crate) inner: TorsionMcRelaxer,
}

#[pymethods]
impl PyTorsionMcRelaxer {
    /// Build a torsion-MC relaxer from a molecule frame.
    ///
    /// Parameters
    /// ----------
    /// frame : molrs.Frame | molpy.Frame
    ///     The chain's frame, carrying an ``"atoms"`` block and its **bond**
    ///     topology. Rotatable bonds (single, acyclic, non-terminal) are
    ///     detected automatically from the graph; their downstream atom sets are
    ///     computed by BFS. Resolved zero-copy through the frame's FFI capsule.
    #[new]
    fn new(frame: &Bound<'_, PyAny>) -> PyResult<Self> {
        let rust_frame = crate::interop::owned_frame_from_py(frame)?;
        let graph = molrs::system::atomistic::Atomistic::from_frame(&rust_frame).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "could not build a molecular graph from the frame (need a bond \
                 topology for torsion MC): {e}"
            ))
        })?;
        Ok(Self {
            inner: TorsionMcRelaxer::new(&graph),
        })
    }

    /// Return a copy with the Metropolis temperature set (reduced units of the
    /// packing objective). Default 1.0.
    fn with_temperature(&self, t: f64) -> Self {
        Self {
            inner: self.inner.clone().with_temperature(t),
        }
    }

    /// Return a copy with the number of MC steps proposed per packing iteration.
    /// Default 10.
    fn with_steps(&self, n: usize) -> Self {
        Self {
            inner: self.inner.clone().with_steps(n),
        }
    }

    /// Return a copy with the maximum per-step rotation (radians). Default π/6.
    fn with_max_delta(&self, rad: f64) -> Self {
        Self {
            inner: self.inner.clone().with_max_delta(rad),
        }
    }

    /// Return a copy with intra-molecular self-avoidance enabled: atom pairs
    /// (excluding 1-2/1-3/1-4 bonded neighbours) closer than ``2 * radius`` get a
    /// quadratic overlap penalty, so the chain folds into physically valid coils
    /// instead of through itself. ``0.0`` disables (default).
    fn with_self_avoidance(&self, radius: f64) -> Self {
        Self {
            inner: self.inner.clone().with_self_avoidance(radius),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "TorsionMcRelaxer(rotatable_bonds={}, steps={}, T={}, max_delta={:.4}, self_avoid={})",
            self.inner.bonds.len(),
            self.inner.steps,
            self.inner.temperature,
            self.inner.max_delta,
            self.inner.self_avoidance_radius,
        )
    }
}

// ── LBFGSRelaxer (force-field geometry relaxer; `ff` feature) ────────────────

#[cfg(feature = "ff")]
use molpack::LBFGSRelaxer;

/// L-BFGS force-field geometry relaxer, attachable to a `Target` via
/// :meth:`Target.with_relaxer`.
#[cfg(feature = "ff")]
#[pyclass(name = "LBFGSRelaxer", from_py_object)]
#[derive(Clone)]
pub struct PyLBFGSRelaxer {
    pub(crate) inner: LBFGSRelaxer,
}

#[cfg(feature = "ff")]
#[pymethods]
impl PyLBFGSRelaxer {
    /// Build a relaxer from a force field.
    ///
    /// Parameters
    /// ----------
    /// forcefield : molrs.ForceField | molpy.ForceField
    ///     The molecule's force field, e.g. from
    ///     ``molpy.io.forcefield.read_lammps("melt.ff")``. Resolved zero-copy
    ///     through its FFI capsule; the potential is compiled lazily against the
    ///     molecule's frame when packing starts.
    #[new]
    fn new(forcefield: &Bound<'_, PyAny>) -> PyResult<Self> {
        // Resolve the molrs/molpy ForceField over the interop bridge and take an
        // owned copy to hand to the (lazy-compiling) core relaxer.
        let ff = crate::interop::forcefield_from_py(forcefield)?.with_forcefield(|ff| ff.clone());
        Ok(Self {
            inner: LBFGSRelaxer::from_forcefield(ff),
        })
    }

    /// Return a copy whose convergence stops when the max per-atom force drops
    /// below `fmax` (kcal/mol/Å). Default 0.05.
    fn with_fmax(&self, fmax: f64) -> Self {
        Self {
            inner: self.inner.clone().with_fmax(fmax),
        }
    }

    /// Return a copy whose L-BFGS is capped at `max_steps` iterations per
    /// relaxation call. Default 500.
    fn with_max_steps(&self, max_steps: usize) -> Self {
        Self {
            inner: self.inner.clone().with_max_steps(max_steps),
        }
    }

    fn __repr__(&self) -> String {
        "LBFGSRelaxer(<force-field L-BFGS geometry relaxer>)".to_string()
    }
}
