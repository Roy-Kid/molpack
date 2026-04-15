//! Python wrapper for the Molpack molecular packer.
//!
//! Provides [`PyPacker`] (a Packmol-grade molecular packer) and
//! [`PyPackResult`] (the output container).

use crate::helpers::{NpF, pack_error_to_pyerr};
use crate::target::PyTarget;
use molpack::F;
use molpack::handler::ProgressHandler;
use molpack::packer::{Molpack, PackResult};
use numpy::IntoPyArray;
use numpy::PyArray2;
use pyo3::prelude::*;

#[pyclass(name = "PackResult", from_py_object)]
#[derive(Clone)]
pub struct PyPackResult {
    inner: PackResult,
}

#[pymethods]
impl PyPackResult {
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<NpF>> {
        let pos = self.inner.positions();
        let n = pos.len();
        let flat: Vec<F> = pos.iter().flat_map(|p| [p[0], p[1], p[2]]).collect();
        let arr = ndarray::Array2::from_shape_vec((n, 3), flat).expect("positions shape");
        arr.into_pyarray(py)
    }

    #[getter]
    fn elements(&self) -> Vec<String> {
        let atoms = self.inner.frame.get("atoms").expect("no atoms block");
        atoms
            .get_string("element")
            .expect("no element column")
            .iter()
            .cloned()
            .collect()
    }

    #[getter]
    fn converged(&self) -> bool {
        self.inner.converged
    }

    #[getter]
    fn fdist(&self) -> F {
        self.inner.fdist
    }

    #[getter]
    fn frest(&self) -> F {
        self.inner.frest
    }

    #[getter]
    fn natoms(&self) -> usize {
        self.inner.natoms()
    }

    fn __repr__(&self) -> String {
        format!(
            "PackResult(converged={}, fdist={:.4}, frest={:.4}, natoms={})",
            self.inner.converged,
            self.inner.fdist,
            self.inner.frest,
            self.inner.natoms()
        )
    }
}

#[pyclass(name = "Packer")]
pub struct PyPacker {
    tolerance: F,
    precision: F,
    maxit: usize,
    nloop0: usize,
    sidemax: F,
    movefrac: F,
    movebadrandom: bool,
    disable_movebad: bool,
    pbc: Option<([F; 3], [F; 3])>,
    progress: bool,
}

#[pymethods]
impl PyPacker {
    #[new]
    #[pyo3(signature = (tolerance=2.0, precision=0.01))]
    fn new(tolerance: NpF, precision: NpF) -> Self {
        PyPacker {
            tolerance,
            precision,
            maxit: 20,
            nloop0: 0,
            sidemax: 1000.0,
            movefrac: 0.05,
            movebadrandom: false,
            disable_movebad: false,
            pbc: None,
            progress: true,
        }
    }

    fn with_tolerance(&self, tolerance: NpF) -> Self {
        PyPacker {
            tolerance,
            ..self.clone_fields()
        }
    }

    fn with_precision(&self, precision: NpF) -> Self {
        PyPacker {
            precision,
            ..self.clone_fields()
        }
    }

    fn with_maxit(&self, maxit: usize) -> Self {
        PyPacker {
            maxit,
            ..self.clone_fields()
        }
    }

    fn with_nloop0(&self, nloop0: usize) -> Self {
        PyPacker {
            nloop0,
            ..self.clone_fields()
        }
    }

    fn with_sidemax(&self, sidemax: NpF) -> Self {
        PyPacker {
            sidemax,
            ..self.clone_fields()
        }
    }

    fn with_movefrac(&self, movefrac: NpF) -> Self {
        PyPacker {
            movefrac,
            ..self.clone_fields()
        }
    }

    fn with_movebadrandom(&self, enabled: bool) -> Self {
        PyPacker {
            movebadrandom: enabled,
            ..self.clone_fields()
        }
    }

    fn with_disable_movebad(&self, disabled: bool) -> Self {
        PyPacker {
            disable_movebad: disabled,
            ..self.clone_fields()
        }
    }

    fn with_pbc(&self, min: [NpF; 3], max: [NpF; 3]) -> Self {
        PyPacker {
            pbc: Some((min, max)),
            ..self.clone_fields()
        }
    }

    fn with_pbc_box(&self, lengths: [NpF; 3]) -> Self {
        self.with_pbc([0.0, 0.0, 0.0], lengths)
    }

    fn with_progress(&self, enabled: bool) -> Self {
        PyPacker {
            progress: enabled,
            ..self.clone_fields()
        }
    }

    #[pyo3(signature = (targets, max_loops=200, seed=None))]
    fn pack(
        &self,
        targets: Vec<PyTarget>,
        max_loops: usize,
        seed: Option<u64>,
    ) -> PyResult<PyPackResult> {
        let rust_targets: Vec<_> = targets.into_iter().map(|t| t.inner).collect();

        let mut packer = Molpack::new()
            .tolerance(self.tolerance)
            .precision(self.precision)
            .maxit(self.maxit)
            .nloop0(self.nloop0)
            .sidemax(self.sidemax)
            .movefrac(self.movefrac)
            .movebadrandom(self.movebadrandom)
            .disable_movebad(self.disable_movebad);

        if let Some((min, max)) = self.pbc {
            packer = packer.pbc(min, max);
        }

        if self.progress {
            packer = packer.add_handler(ProgressHandler::new());
        }

        let result = packer
            .pack(&rust_targets, max_loops, seed)
            .map_err(pack_error_to_pyerr)?;

        Ok(PyPackResult { inner: result })
    }

    fn __repr__(&self) -> String {
        format!(
            "Packer(tolerance={:.2}, precision={:.4})",
            self.tolerance, self.precision
        )
    }
}

impl PyPacker {
    fn clone_fields(&self) -> PyPacker {
        PyPacker {
            tolerance: self.tolerance,
            precision: self.precision,
            maxit: self.maxit,
            nloop0: self.nloop0,
            sidemax: self.sidemax,
            movefrac: self.movefrac,
            movebadrandom: self.movebadrandom,
            disable_movebad: self.disable_movebad,
            pbc: self.pbc,
            progress: self.progress,
        }
    }
}
