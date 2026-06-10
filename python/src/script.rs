//! Python binding for the molpack script loader.
//!
//! Exposes a single function :func:`load_script` that parses an `.inp`
//! script and returns a ready-to-run :class:`Molpack` plus target list.
//! Everything downstream — attaching handlers, running ``pack()``,
//! writing output — stays in Python hands.
//!
//! The loader does **not** touch molecule files in Rust. Each
//! ``structure``'s template is read on the Python side, defaulting to
//! :mod:`molrs` (``molrs.read_pdb`` / ``read_xyz``) but pluggable via the
//! ``read_frame`` argument. This keeps the PyO3 wheel free of
//! ``molrs-io`` and lets users plug in their own loader (mdtraj, ASE,
//! plain dicts, etc.).

use std::path::PathBuf;

use pyo3::exceptions::{PyImportError, PyOSError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use molpack::script::{self, ScriptPlan, StructurePlan};

use crate::helpers::script_error_to_pyerr;
use crate::packer::PyPacker;
use crate::target::{PyTarget, target_from_frame};

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
/// read_frame : callable, optional
///     Callable ``(path: str, filetype: str | None) -> Frame`` used to
///     load each ``structure`` template. The returned object only needs
///     a ``frame["atoms"]`` block exposing ``x`` / ``y`` / ``z`` and an
///     ``element`` (or ``symbol``) column — :class:`molrs.Frame` and
///     plain dicts both work. Defaults to :mod:`molrs`'s ``read_pdb`` /
///     ``read_xyz`` dispatched by file extension.
///
/// Returns
/// -------
/// ScriptJob
///     Bundle of ``(packer, targets, output, nloop)`` — supports both
///     attribute access and tuple unpacking.
#[pyfunction]
#[pyo3(signature = (path, *, read_frame = None))]
pub fn load_script(
    py: Python<'_>,
    path: PathBuf,
    read_frame: Option<Py<PyAny>>,
) -> PyResult<PyScriptJob> {
    let src = std::fs::read_to_string(&path)
        .map_err(|e| PyOSError::new_err(format!("reading {}: {e}", path.display())))?;

    let script_ast = script::parse(&src).map_err(script_error_to_pyerr)?;

    let base_dir = path
        .canonicalize()
        .unwrap_or_else(|_| path.clone())
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));

    let plan: ScriptPlan = script_ast.lower(&base_dir).map_err(script_error_to_pyerr)?;

    let loader = match read_frame {
        Some(cb) => cb,
        None => default_molrs_loader(py)?,
    };

    let targets: Vec<PyTarget> = plan
        .structures
        .iter()
        .map(|sp| build_target(py, sp, plan.filetype.as_deref(), &loader))
        .collect::<PyResult<_>>()?;

    let mut packer = PyPacker::default();
    packer.tolerance = Some(script_ast.tolerance);
    packer.seed = script_ast.seed;
    if let Some(pbc) = script_ast.pbc {
        packer.periodic_box = Some((pbc.min, pbc.max));
    }

    Ok(PyScriptJob {
        packer: Py::new(py, packer)?,
        targets,
        output: plan.output,
        nloop: plan.nloop,
    })
}

/// Call the user-supplied loader with `(path, filetype)`, then build a
/// [`PyTarget`] from the returned frame and stamp on the structure's
/// restraints / centering / fixed pose.
fn build_target(
    py: Python<'_>,
    sp: &StructurePlan,
    filetype: Option<&str>,
    loader: &Py<PyAny>,
) -> PyResult<PyTarget> {
    let path_str = sp.filepath.to_string_lossy().into_owned();
    let frame_obj = loader.bind(py).call1((path_str, filetype))?;

    let target = target_from_frame(&frame_obj, sp.number)?;
    Ok(PyTarget {
        inner: sp.apply(target),
    })
}

/// Build a default loader: import :mod:`molrs` and dispatch on the
/// file extension or the script's ``filetype`` keyword.
///
/// Returns a Python callable; failures (e.g. ``molrs`` not installed)
/// surface as :class:`ImportError` from the script-loading site.
fn default_molrs_loader(py: Python<'_>) -> PyResult<Py<PyAny>> {
    let molrs = py.import("molrs").map_err(|e| {
        PyImportError::new_err(format!(
            "loading template files needs `molcrafts-molrs` (or pass read_frame=...): {e}"
        ))
    })?;

    let globals = PyDict::new(py);
    globals.set_item("molrs", molrs)?;

    let code = pyo3::ffi::c_str!(
        r#"
def _loader(path, filetype):
    fmt = (filetype or '').lower() or path.rsplit('.', 1)[-1].lower()
    if fmt == 'pdb':
        return molrs.read_pdb(path)
    if fmt == 'xyz':
        return molrs.read_xyz(path)
    raise ValueError(
        f"default loader handles .pdb / .xyz only - pass read_frame=... for {fmt!r}"
    )
"#
    );
    py.run(code, Some(&globals), None)
        .map_err(|e| PyValueError::new_err(format!("failed to build default frame loader: {e}")))?;

    let loader = globals.get_item("_loader")?.ok_or_else(|| {
        PyValueError::new_err("internal error: _loader missing from default-loader globals")
    })?;
    Ok(loader.unbind())
}
