//! Python wrapper for packing `Target`.
//!
//! [`PyTarget`] describes one type of molecule to pack: its template
//! geometry, topology, and the number of copies.
//!
//! The constructor accepts a real ``molrs.Frame`` (``molpy.Frame`` is the same
//! class) carrying an ``"atoms"`` block. The frame crosses the
//! language boundary **zero-copy** through its stable-FFI capsule (see
//! [`crate::interop`]) — no dict marshalling, no consumer-side data type. The
//! full frame, with topology, is handed to the core [`Target`], which owns the
//! assembly.

use crate::constraint::{extract_collective_restraint, extract_restraint, try_region};
use crate::helpers::NpF;
use crate::types::{PyAngle, PyAxis, PyCenteringMode};
use molpack::F;
use molpack::target::Target;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

/// Build a [`Target`] from any frame-like Python object plus a copy count.
///
/// Shared by [`PyTarget::new`] and the script loader. The frame is converted to
/// a Rust [`molrs::Frame`] so the core retains its full topology.
pub(crate) fn target_from_frame(frame: &Bound<'_, PyAny>, count: usize) -> PyResult<Target> {
    let rust_frame = crate::interop::owned_frame_from_py(frame)?;
    // `Target::new` panics on a frame without float coordinates; answer that
    // here as a Python error instead of a panic across the boundary.
    rust_frame
        .coords()
        .map_err(|e| PyValueError::new_err(format!("frame has no coordinates: {e}")))?;
    Ok(Target::new(rust_frame, count))
}

#[pyclass(name = "Target", from_py_object)]
#[derive(Clone)]
pub struct PyTarget {
    pub(crate) inner: Target,
}

#[pymethods]
impl PyTarget {
    /// Create a packing target from a molecule frame.
    ///
    /// Parameters
    /// ----------
    /// frame : molrs.Frame
    ///     A frame with an ``"atoms"`` block (``x`` / ``y`` / ``z``
    ///     columns). Resolved zero-copy via its FFI capsule — a plain ``dict``
    ///     is no longer accepted; build a ``molrs.Frame`` first.
    /// count : int
    ///     Number of copies to pack.
    ///
    /// Add a display name via :meth:`with_name`.
    #[new]
    #[pyo3(signature = (frame, count))]
    fn new(frame: &Bound<'_, PyAny>, count: usize) -> PyResult<Self> {
        Ok(PyTarget {
            inner: target_from_frame(frame, count)?,
        })
    }

    fn with_name(&self, name: &str) -> Self {
        PyTarget {
            inner: self.inner.clone().with_name(name),
        }
    }

    /// Override the per-copy total mass (amu) used by
    /// ``with_density`` when element symbols cannot provide one.
    fn with_mass(&self, amu: crate::helpers::NpF) -> Self {
        PyTarget {
            inner: self.inner.clone().with_mass(amu),
        }
    }

    /// Attach a restraint to this target — the single unified extension point.
    ///
    /// Accepts:
    ///
    /// * a molrs **region** (``molrs.Sphere``, ``Cuboid``, ``Parallelepiped``,
    ///   ``HalfSpace``, ``Cylinder``, ``Ellipsoid``, ``Polyhedron``,
    ///   ``SphereUnion``, or a ``&`` / ``|`` / ``~`` composition) — lifted
    ///   through ``RegionRestraint``, so every atom must stay inside it;
    /// * a built-in **distribution** restraint (:class:`GaussianPlane`,
    ///   :class:`GaussianPoint`, :class:`ExponentialPlane`,
    ///   :class:`ExponentialPoint`, :class:`TabulatedPlane`,
    ///   :class:`TabulatedPoint`);
    /// * any object with callable ``f`` / ``fg`` — the duck-typed extension
    ///   point.
    ///
    /// A region is evaluated per atom. For the
    /// distribution and custom (duck-typed) restraints,
    /// ``f(coords, scale, scale2)`` / ``fg(coords, scale, scale2)`` see **every
    /// copy's** ``(x, y, z)`` (``coords`` is the full list) and ``fg`` returns
    /// ``(energy, [(gx, gy, gz), ...])`` — one gradient triple per copy.
    fn with_restraint(&self, restraint: &Bound<'_, pyo3::types::PyAny>) -> PyResult<Self> {
        if let Some(atom_r) = try_region(restraint)? {
            Ok(PyTarget {
                inner: self.inner.clone().with_restraint(atom_r),
            })
        } else {
            let group_r = extract_collective_restraint(restraint)?;
            Ok(PyTarget {
                inner: self.inner.clone().with_collective_restraint(group_r),
            })
        }
    }

    /// Attach a restraint to selected atoms of every copy.
    ///
    /// ``indices`` are **0-based** (Rust/Python native). Porting from a
    /// Packmol ``.inp`` file? Subtract 1 at the call site.
    fn with_atom_restraint(
        &self,
        indices: Vec<usize>,
        restraint: &Bound<'_, pyo3::types::PyAny>,
    ) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        let r = extract_restraint(restraint)?;
        Ok(PyTarget {
            inner: self.inner.clone().with_atom_restraint(&indices, r),
        })
    }

    /// Set the packing radius for **every atom** of this target.
    ///
    /// Packmol's structure-level ``radius``. The packer separates two atoms by
    /// the sum of their radii; without this every atom uses the global
    /// ``tolerance / 2``. Van der Waals radii from the source file are not used
    /// as packing radii.
    ///
    /// Raises ``ValueError`` if ``radius`` is not positive.
    fn with_radius(&self, radius: F) -> PyResult<Self> {
        check_positive(radius, "packing radius")?;
        Ok(PyTarget {
            inner: self.inner.clone().with_radius(radius),
        })
    }

    /// Set the packing radius for selected atoms of every copy.
    ///
    /// Packmol's ``radius`` inside an ``atoms ... end atoms`` block.
    /// ``indices`` are **0-based** (Rust/Python native); a Packmol ``.inp``
    /// uses 1-based indices, so subtract 1 when porting.
    ///
    /// Raises ``ValueError`` if ``radius`` is not positive or an index is out
    /// of range.
    fn with_atom_radius(&self, indices: Vec<usize>, radius: F) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        check_positive(radius, "packing radius")?;
        Ok(PyTarget {
            inner: self.inner.clone().with_atom_radius(&indices, radius),
        })
    }

    /// Weight this target's atoms in the overlap penalty (Packmol ``fscale``).
    ///
    /// The pair term is multiplied by ``fscale_i * fscale_j``, so a value below
    /// 1 makes a species *softer* without changing the distance it is asked to
    /// keep. Default ``1.0``.
    ///
    /// Raises ``ValueError`` if ``fscale`` is not positive.
    fn with_fscale(&self, fscale: F) -> PyResult<Self> {
        check_positive(fscale, "fscale")?;
        Ok(PyTarget {
            inner: self.inner.clone().with_fscale(fscale),
        })
    }

    /// Weight selected atoms in the overlap penalty. ``indices`` are **0-based**.
    fn with_atom_fscale(&self, indices: Vec<usize>, fscale: F) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        check_positive(fscale, "fscale")?;
        Ok(PyTarget {
            inner: self.inner.clone().with_atom_fscale(&indices, fscale),
        })
    }

    /// Give this target's atoms a second, shorter penalty radius
    /// (Packmol ``short_radius``). Must be smaller than the packing radius.
    fn with_short_radius(&self, short_radius: F) -> PyResult<Self> {
        check_positive(short_radius, "short radius")?;
        Ok(PyTarget {
            inner: self.inner.clone().with_short_radius(short_radius),
        })
    }

    /// Per-atom counterpart of :meth:`with_short_radius`; ``indices`` are
    /// **0-based**.
    fn with_atom_short_radius(&self, indices: Vec<usize>, short_radius: F) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        check_positive(short_radius, "short radius")?;
        Ok(PyTarget {
            inner: self
                .inner
                .clone()
                .with_atom_short_radius(&indices, short_radius),
        })
    }

    /// Weight the short-radius penalty (Packmol ``short_radius_scale``).
    fn with_short_radius_scale(&self, scale: F) -> PyResult<Self> {
        check_positive(scale, "short radius scale")?;
        Ok(PyTarget {
            inner: self.inner.clone().with_short_radius_scale(scale),
        })
    }

    /// Per-atom counterpart of :meth:`with_short_radius_scale`; ``indices`` are
    /// **0-based**.
    fn with_atom_short_radius_scale(&self, indices: Vec<usize>, scale: F) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        check_positive(scale, "short radius scale")?;
        Ok(PyTarget {
            inner: self
                .inner
                .clone()
                .with_atom_short_radius_scale(&indices, scale),
        })
    }

    /// Set this target's intramolecular skip table.
    ///
    /// Slot 0 is the 1-2 weight; the last slot is the 1-N tail. Default is
    /// ``[0, 0, 0, 1]`` (depth 3: 1-2/1-3/1-4 exempt). This is not a
    /// force-field ``special_bonds`` triple.
    ///
    /// Fractional weights (Amber 1-4 ``0.5``) are stored here and refused
    /// later when growth compiles the skip set. All-atom explicit hydrogen
    /// keeps the default table and shrinks hydrogen via
    /// :meth:`with_atom_radius`.
    ///
    /// Raises ``ValueError`` if the table is empty or a weight is outside
    /// ``[0, 1]`` or not finite.
    fn with_special_bonds(&self, table: Vec<NpF>) -> PyResult<Self> {
        let table = validate_special_bonds(table)?;
        Ok(PyTarget {
            inner: self.inner.clone().with_special_bonds(table),
        })
    }

    /// Name the atoms growth treats as hydrogens (**0-based**), replacing the
    /// default rule (element symbol ``H``). Lattice growth places hydrogens
    /// off their backbone neighbour, never on a lattice site; a
    /// coarse-grained model with no hydrogens passes ``[]``.
    ///
    /// Raises ``ValueError`` if an index is out of range.
    fn with_hydrogens(&self, indices: Vec<usize>) -> PyResult<Self> {
        validate_atom_indices(&indices, self.inner.natoms())?;
        Ok(PyTarget {
            inner: self.inner.clone().with_hydrogens(&indices),
        })
    }

    fn with_perturb_budget(&self, budget: usize) -> Self {
        PyTarget {
            inner: self.inner.clone().with_perturb_budget(budget),
        }
    }

    /// One fixed obstacle target holding a previous pack's entire output,
    /// coordinates kept verbatim — the named chaining primitive: grow first,
    /// then pack the next stage around the frozen matrix.
    #[staticmethod]
    fn fixed_from(result: &crate::result::PyState) -> Self {
        Self {
            inner: molpack::Target::fixed_from(&result.inner),
        }
    }

    fn with_centering(&self, mode: PyCenteringMode) -> Self {
        PyTarget {
            inner: self.inner.clone().with_centering(mode.into()),
        }
    }

    fn with_rotation_bound(&self, axis: PyAxis, center: PyAngle, half_width: PyAngle) -> Self {
        PyTarget {
            inner: self.inner.clone().with_rotation_bound(
                axis.into(),
                center.inner,
                half_width.inner,
            ),
        }
    }

    fn fixed_at(&self, position: [NpF; 3]) -> Self {
        PyTarget {
            inner: self.inner.clone().fixed_at(position),
        }
    }

    /// Apply an Euler orientation to a previously-fixed target.
    ///
    /// Takes a 3-tuple of :class:`Angle` — e.g.
    /// ``(Angle.from_degrees(45), Angle.ZERO, Angle.from_degrees(90))``.
    /// Must be called after :meth:`fixed_at`.
    fn with_orientation(&self, orientation: (PyAngle, PyAngle, PyAngle)) -> Self {
        PyTarget {
            inner: self.inner.clone().with_orientation([
                orientation.0.inner,
                orientation.1.inner,
                orientation.2.inner,
            ]),
        }
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }

    #[getter]
    fn natoms(&self) -> usize {
        self.inner.natoms()
    }

    #[getter]
    fn count(&self) -> usize {
        self.inner.count
    }

    #[getter]
    fn elements(&self) -> Vec<String> {
        self.inner.elements.clone()
    }

    #[getter]
    fn radii(&self) -> Vec<F> {
        self.inner.radii.clone()
    }

    #[getter]
    fn special_bonds(&self) -> Vec<F> {
        self.inner.special_bonds.as_slice().to_vec()
    }

    #[getter]
    fn is_fixed(&self) -> bool {
        self.inner.fixed_at.is_some()
    }

    fn __repr__(&self) -> String {
        format!(
            "Target(natoms={}, count={}, name={:?})",
            self.inner.natoms(),
            self.inner.count,
            self.inner.name
        )
    }
}

fn validate_atom_indices(indices: &[usize], natoms: usize) -> PyResult<()> {
    for &index in indices {
        if index >= natoms {
            return Err(PyValueError::new_err(format!(
                "atom indices are 0-based and must be in 0..{natoms}, got {index}",
            )));
        }
    }
    Ok(())
}

/// Reject a non-positive per-atom property value with a Python `ValueError`.
fn check_positive(value: F, what: &str) -> PyResult<()> {
    if value <= 0.0 || value.is_nan() {
        return Err(PyValueError::new_err(format!(
            "{what} must be positive, got {value}"
        )));
    }
    Ok(())
}

/// Marshal a Python weight list into [`molpack::BondDistanceWeights`].
///
/// Rejects empty, non-finite, or out-of-range entries with ``ValueError``.
/// Fractional weights are legal here; growth refuses them later.
fn validate_special_bonds(weights: Vec<NpF>) -> PyResult<molpack::BondDistanceWeights> {
    molpack::BondDistanceWeights::new(weights).map_err(|e| PyValueError::new_err(e.to_string()))
}
