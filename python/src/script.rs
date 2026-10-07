//! Python binding for the molpack script loader.
//!
//! Exposes a single function :func:`load_script` that parses an `.inp`
//! script and returns a ready-to-run :class:`GencanPack` plus target list.
//! Everything downstream — attaching callbacks, running ``run()``,
//! writing output — stays in Python hands.
//!
//! The loader does **not** touch molecule files in Rust. Each
//! ``structure``'s template is read on the Python side, defaulting to the
//! ``molrs.io`` reader of its [`StructureFormat`] (the script's ``filetype``,
//! else the file name: ``molrs.io.read_pdb``, ``read_xyz``, …) but pluggable
//! via the ``read_frame`` argument. This keeps the PyO3 wheel free of
//! ``molrs-io`` and lets users plug in their own loader (mdtraj, ASE, …) as
//! long as it returns a ``molrs.core.Frame``.

use std::path::PathBuf;

use pyo3::exceptions::{PyImportError, PyOSError};
use pyo3::prelude::*;
use pyo3::types::PyModule;

use molpack::script::{self, ScriptPlan, StructureFormat, StructurePlan};

use crate::errors::script_error_to_pyerr;
use crate::packing_methods::PyGencanPack;
use crate::target::{PyTarget, target_from_frame};

/// Output of [`load_script`] — four fields bundled as a PyClass so
/// Python callers access them by attribute *and* iterate it for
/// tuple-unpacking (``packer, targets, output, nloop = load_script(...)``).
#[pyclass(name = "ScriptJob", module = "molpack", sequence)]
pub struct PyScriptJob {
    /// `GencanPack` pre-configured with ``tolerance`` / ``seed`` /
    /// periodic box from the script.
    #[pyo3(get)]
    pub packer: Py<PyGencanPack>,
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
///     ``element`` column. Must be a :class:`molrs.core.Frame` (``molpy.Frame``
///     is the same class). Defaults to the ``molrs.io`` reader of the
///     format ``filetype`` or the file name names (``read_pdb``,
///     ``read_xyz``, …).
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
        Some(callable) => TemplateLoader::Callable(callable),
        None => TemplateLoader::Molrs(import_molrs_io(py)?),
    };

    let targets: Vec<PyTarget> = plan
        .structures
        .iter()
        .map(|sp| build_target(py, sp, plan.filetype.as_deref(), &loader))
        .collect::<PyResult<_>>()?;

    let packer = PyGencanPack::from_script(
        Some(script_ast.tolerance),
        script_ast.seed,
        script_ast.pbc.map(|pbc| (pbc.min, pbc.max)),
    );

    Ok(PyScriptJob {
        packer: Py::new(py, packer)?,
        targets,
        output: plan.output,
        nloop: plan.nloop,
    })
}

/// How a template file becomes a frame: the caller's callable, or the
/// `molrs.io` reader of the file's format.
enum TemplateLoader<'py> {
    Callable(Py<PyAny>),
    Molrs(Bound<'py, PyModule>),
}

/// Read the structure's template through `loader`, then build a
/// [`PyTarget`] from the returned frame and stamp on the structure's
/// restraints / centering / fixed pose.
fn build_target<'py>(
    py: Python<'py>,
    sp: &StructurePlan,
    filetype: Option<&str>,
    loader: &TemplateLoader<'py>,
) -> PyResult<PyTarget> {
    let path_str = sp.filepath.to_string_lossy().into_owned();
    let frame_obj = match loader {
        TemplateLoader::Callable(callable) => callable.bind(py).call1((path_str, filetype))?,
        TemplateLoader::Molrs(io) => {
            let format =
                StructureFormat::resolve(&sp.filepath, filetype).map_err(script_error_to_pyerr)?;
            read_with_molrs(io, format, &path_str)?
        }
    };

    let target = target_from_frame(&frame_obj, sp.number)?;
    Ok(PyTarget {
        inner: sp.apply(target),
    })
}

/// ``molrs.io``, the default loader's home.
///
/// Failures (e.g. ``molrs`` not installed) surface as :class:`ImportError`
/// from the script-loading site.
fn import_molrs_io(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    py.import("molrs.io").map_err(|e| {
        PyImportError::new_err(format!(
            "loading template files needs `molcrafts-molrs` (or pass read_frame=...): {e}"
        ))
    })
}

/// Read the first structure of `path` through the ``molrs.io`` reader of
/// `format` — the Python twin of `StructureFormat::read`.
fn read_with_molrs<'py>(
    io: &Bound<'py, PyModule>,
    format: StructureFormat,
    path: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let reader = match format {
        StructureFormat::Pdb => "read_pdb",
        StructureFormat::Xyz => "read_xyz",
        StructureFormat::Sdf => "read_sdf",
        StructureFormat::Mol2 => "read_mol2",
        StructureFormat::Gro => "read_gro",
        StructureFormat::Cif => "read_cif",
        StructureFormat::VaspPoscar => "read_vasp_poscar",
        StructureFormat::Xsf => "read_xsf",
        StructureFormat::Cube => "read_cube",
        StructureFormat::AmberInpcrd => "read_amber_inpcrd",
        StructureFormat::LammpsData => "read_lammps_data",
        StructureFormat::LammpsDump => {
            // A dump is a trajectory: its first snapshot is the template.
            let trajectory = io.getattr("read_lammps_trajectory")?.call1((path,))?;
            return trajectory.call_method1("read_frame", (0,));
        }
    };
    io.getattr(reader)?.call1((path,))
}
