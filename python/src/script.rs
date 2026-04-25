//! Python binding for the molpack script loader.
//!
//! Exposes a single function :func:`load_script` that parses an `.inp`
//! script and returns a ready-to-run :class:`Molpack` plus target list.
//! Everything downstream — attaching handlers, running ``pack()``,
//! writing output — stays in Python hands.

use std::path::PathBuf;

use pyo3::exceptions::PyOSError;
use pyo3::prelude::*;

use molpack::script;

use crate::helpers::script_error_to_pyerr;
use crate::packer::PyPacker;
use crate::target::PyTarget;

/// Output of [`load_script`] — four fields bundled as a PyClass so
/// Python callers access them by attribute *and* iterate it for
/// tuple-unpacking (``packer, targets, output, nloop = load_script(...)``).
#[pyclass(name = "ScriptJob", module = "molpack", sequence)]
pub struct PyScriptJob {
    /// Packer pre-configured with ``tolerance`` and ``seed`` from the script.
    #[pyo3(get)]
    pub packer: Py<PyPacker>,
    /// Targets ready to be packed.
    #[pyo3(get)]
    pub targets: Vec<PyTarget>,
    /// Resolved output file path (relative paths are resolved against the
    /// script's parent directory).
    #[pyo3(get)]
    pub output: PathBuf,
    /// Outer-loop iteration cap (``nloop`` keyword; default 400).
    #[pyo3(get)]
    pub nloop: usize,
}

#[pymethods]
impl PyScriptJob {
    fn __repr__(&self) -> String {
        format!(
            "ScriptJob(targets={}, output={:?}, nloop={})",
            self.targets.len(),
            self.output,
            self.nloop
        )
    }

    fn __len__(&self) -> usize {
        4
    }

    fn __getitem__<'py>(&self, py: Python<'py>, idx: isize) -> PyResult<Bound<'py, PyAny>> {
        let i = if idx < 0 { idx + 4 } else { idx };
        match i {
            0 => Ok(self.packer.clone_ref(py).into_bound(py).into_any()),
            1 => Ok(self.targets.clone().into_pyobject(py)?.into_any()),
            2 => Ok(self.output.clone().into_pyobject(py)?.into_any()),
            3 => Ok(self.nloop.into_pyobject(py)?.into_any()),
            _ => Err(pyo3::exceptions::PyIndexError::new_err(
                "ScriptJob index out of range (0..4)",
            )),
        }
    }
}

/// Parse and lower a molpack `.inp` script.
///
/// Parameters
/// ----------
/// path : str | pathlib.Path
///     Path to a ``.inp`` script. Relative file paths inside the script
///     (structures, output) are resolved against the script's parent
///     directory.
///
/// Returns
/// -------
/// ScriptJob
///     Bundle of ``(packer, targets, output, nloop)`` — supports both
///     attribute access and tuple unpacking.
#[pyfunction]
pub fn load_script(py: Python<'_>, path: PathBuf) -> PyResult<PyScriptJob> {
    let src = std::fs::read_to_string(&path)
        .map_err(|e| PyOSError::new_err(format!("reading {}: {e}", path.display())))?;

    let script_ast = script::parse(&src).map_err(script_error_to_pyerr)?;

    let base_dir = path
        .canonicalize()
        .unwrap_or_else(|_| path.clone())
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));

    let built = script_ast.build(&base_dir).map_err(script_error_to_pyerr)?;

    let mut packer = PyPacker::default();
    packer.tolerance = Some(script_ast.tolerance);
    packer.seed = script_ast.seed;
    if let Some(pbc) = script_ast.pbc {
        packer.periodic_box = Some((pbc.min, pbc.max));
    }

    let targets: Vec<PyTarget> = built
        .targets
        .into_iter()
        .map(|inner| PyTarget { inner })
        .collect();

    Ok(PyScriptJob {
        packer: Py::new(py, packer)?,
        targets,
        output: built.output,
        nloop: built.nloop,
    })
}
