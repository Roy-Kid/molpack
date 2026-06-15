//! Python wrappers for molecular packing restraints.
//!
//! Each built-in geometric restraint is a `#[pyclass]` named `*Restraint` to
//! mirror the Rust type names. Custom Python-defined restraints are supported via
//! **duck typing**: any object exposing callable ``f(coords, scale, scale2)`` and
//! ``fg(coords, scale, scale2)`` attributes may be passed to
//! ``Target.with_restraint``; see [`PyCallableRestraint`] for the group contract.

use std::sync::Arc;

use crate::helpers::{NpF, stash_err};
use molpack::F;
use molpack::restraint::{
    AbovePlaneRestraint, AtomRestraint, BelowPlaneRestraint, ExponentialPlane, ExponentialPoint,
    GaussianPlane, GaussianPoint, InsideBoxRestraint, InsideSphereRestraint,
    OutsideSphereRestraint, Restraint, TabulatedPlane, TabulatedPoint,
};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

// Pass-through wrapper: an owned `Arc<dyn AtomRestraint>` that itself
// implements `AtomRestraint`, so it can be fed into
// `Target::with_restraint(impl AtomRestraint)`. Adding a new restraint type
// only requires a new arm in `extract_restraint`.

#[derive(Clone)]
pub(crate) struct SharedAtomRestraint(pub Arc<dyn AtomRestraint>);

impl std::fmt::Debug for SharedAtomRestraint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SharedAtomRestraint({})", self.0.name())
    }
}

impl AtomRestraint for SharedAtomRestraint {
    #[inline]
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F {
        self.0.f(x, scale, scale2)
    }
    #[inline]
    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        self.0.fg(x, scale, scale2, g)
    }
    #[inline]
    fn is_parallel_safe(&self) -> bool {
        self.0.is_parallel_safe()
    }
    #[inline]
    fn name(&self) -> &'static str {
        self.0.name()
    }
    #[inline]
    fn periodic_box(&self) -> Option<([F; 3], [F; 3], [bool; 3])> {
        self.0.periodic_box()
    }
}

// ============================================================================
// Extractor: try each built-in `#[pyclass]`, else duck-type on `f`/`fg`.
// ============================================================================

pub(crate) fn extract_restraint(
    obj: &Bound<'_, pyo3::types::PyAny>,
) -> PyResult<SharedAtomRestraint> {
    if let Ok(c) = obj.extract::<PyInsideBoxRestraint>() {
        return Ok(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyInsideSphereRestraint>() {
        return Ok(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyOutsideSphereRestraint>() {
        return Ok(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyAbovePlaneRestraint>() {
        return Ok(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyBelowPlaneRestraint>() {
        return Ok(SharedAtomRestraint(Arc::new(c.inner)));
    }

    // Duck-typed Python restraint: object with callable `f` and `fg`
    // methods. Bound methods are resolved once here so the hot path
    // skips per-call attribute lookups.
    if let (Ok(f_method), Ok(fg_method)) = (obj.getattr("f"), obj.getattr("fg")) {
        return Ok(SharedAtomRestraint(Arc::new(PyCallableAtomRestraint {
            f_method: f_method.unbind(),
            fg_method: fg_method.unbind(),
        })));
    }

    Err(PyTypeError::new_err(
        "expected a restraint: one of InsideBoxRestraint / InsideSphereRestraint / \
         OutsideSphereRestraint / AbovePlaneRestraint / BelowPlaneRestraint, or an \
         object with callable `f(x, scale, scale2)` and `fg(x, scale, scale2)` methods",
    ))
}

/// Try ONLY the built-in geometric per-atom pyclasses (no duck-typing).
///
/// Returns `Some(..)` if `obj` is one of the five geometric built-ins
/// (`InsideBoxRestraint`, `InsideSphereRestraint`, `OutsideSphereRestraint`,
/// `AbovePlaneRestraint`, `BelowPlaneRestraint`); otherwise `None`. The unified
/// [`crate::target::PyTarget::with_restraint`] entry point uses this to route a
/// geometric built-in to the per-atom path and everything else (built-in
/// distribution restraints + duck-typed objects) to the group path.
pub(crate) fn try_atom_builtin(obj: &Bound<'_, pyo3::types::PyAny>) -> Option<SharedAtomRestraint> {
    if let Ok(c) = obj.extract::<PyInsideBoxRestraint>() {
        return Some(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyInsideSphereRestraint>() {
        return Some(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyOutsideSphereRestraint>() {
        return Some(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyAbovePlaneRestraint>() {
        return Some(SharedAtomRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyBelowPlaneRestraint>() {
        return Some(SharedAtomRestraint(Arc::new(c.inner)));
    }
    None
}

// PyCallableAtomRestraint — bridge from the Rust `AtomRestraint` trait to a
// Python object. Stores the **bound methods** directly (resolved at
// attach time) instead of the host object, because each restraint
// evaluation goes through `fg` inside the GENCAN inner loop — a per-atom
// string lookup there is measurable on larger systems.
//
// Python contract:
//   obj.f(x, scale, scale2)  -> float
//   obj.fg(x, scale, scale2) -> (float, (gx, gy, gz))
//
// `x` is passed as a 3-tuple; `fg`'s returned gradient is a flat tuple
// that Rust accumulates into `g` with `+=` on the caller's behalf.
//
// `is_parallel_safe() -> false` — the GIL serializes callbacks, and
// pretending otherwise would deadlock rayon reductions.

pub(crate) struct PyCallableAtomRestraint {
    f_method: Py<pyo3::types::PyAny>,
    fg_method: Py<pyo3::types::PyAny>,
}

impl std::fmt::Debug for PyCallableAtomRestraint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PyCallableAtomRestraint")
    }
}

impl AtomRestraint for PyCallableAtomRestraint {
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F {
        Python::attach(|py| {
            let args = ((x[0], x[1], x[2]), scale, scale2);
            match self.f_method.bind(py).call1(args) {
                Ok(res) => match res.extract::<F>() {
                    Ok(v) => v,
                    Err(e) => {
                        stash_err(e);
                        0.0
                    }
                },
                Err(e) => {
                    stash_err(e);
                    0.0
                }
            }
        })
    }

    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        Python::attach(|py| {
            let args = ((x[0], x[1], x[2]), scale, scale2);
            match self.fg_method.bind(py).call1(args) {
                Ok(res) => match res.extract::<(F, (F, F, F))>() {
                    Ok((v, (gx, gy, gz))) => {
                        g[0] += gx;
                        g[1] += gy;
                        g[2] += gz;
                        v
                    }
                    Err(e) => {
                        stash_err(PyTypeError::new_err(format!(
                            "callable restraint `fg` must return (float, (gx, gy, gz)); {e}",
                        )));
                        0.0
                    }
                },
                Err(e) => {
                    stash_err(e);
                    0.0
                }
            }
        })
    }

    fn is_parallel_safe(&self) -> bool {
        false
    }

    fn name(&self) -> &'static str {
        "PyCallableAtomRestraint"
    }
}

// ============================================================================
// Built-in geometric restraint `#[pyclass]` wrappers. Names match the Rust
// types (`*Restraint` suffix); parameter order mirrors the Rust constructors.
// ============================================================================

#[pyclass(name = "InsideBoxRestraint", from_py_object)]
#[derive(Clone)]
pub struct PyInsideBoxRestraint {
    pub(crate) inner: InsideBoxRestraint,
}

#[pymethods]
impl PyInsideBoxRestraint {
    #[new]
    #[pyo3(signature = (min, max, periodic=(false, false, false)))]
    fn new(min: [NpF; 3], max: [NpF; 3], periodic: (bool, bool, bool)) -> Self {
        Self {
            inner: InsideBoxRestraint::new(min, max, [periodic.0, periodic.1, periodic.2]),
        }
    }

    fn __repr__(&self) -> String {
        let p = self.inner.periodic;
        format!(
            "InsideBoxRestraint(min={:?}, max={:?}, periodic=({}, {}, {}))",
            self.inner.min, self.inner.max, p[0], p[1], p[2],
        )
    }
}

#[pyclass(name = "InsideSphereRestraint", from_py_object)]
#[derive(Clone)]
pub struct PyInsideSphereRestraint {
    pub(crate) inner: InsideSphereRestraint,
}

#[pymethods]
impl PyInsideSphereRestraint {
    #[new]
    #[pyo3(signature = (center, radius))]
    fn new(center: [NpF; 3], radius: NpF) -> Self {
        Self {
            inner: InsideSphereRestraint::new(center, radius),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "InsideSphereRestraint(center={:?}, radius={})",
            self.inner.center, self.inner.radius,
        )
    }
}

#[pyclass(name = "OutsideSphereRestraint", from_py_object)]
#[derive(Clone)]
pub struct PyOutsideSphereRestraint {
    pub(crate) inner: OutsideSphereRestraint,
}

#[pymethods]
impl PyOutsideSphereRestraint {
    #[new]
    #[pyo3(signature = (center, radius))]
    fn new(center: [NpF; 3], radius: NpF) -> Self {
        Self {
            inner: OutsideSphereRestraint::new(center, radius),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "OutsideSphereRestraint(center={:?}, radius={})",
            self.inner.center, self.inner.radius,
        )
    }
}

#[pyclass(name = "AbovePlaneRestraint", from_py_object)]
#[derive(Clone)]
pub struct PyAbovePlaneRestraint {
    pub(crate) inner: AbovePlaneRestraint,
}

#[pymethods]
impl PyAbovePlaneRestraint {
    #[new]
    #[pyo3(signature = (normal, distance))]
    fn new(normal: [NpF; 3], distance: NpF) -> Self {
        Self {
            inner: AbovePlaneRestraint::new(normal, distance),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AbovePlaneRestraint(normal={:?}, distance={})",
            self.inner.normal, self.inner.distance,
        )
    }
}

#[pyclass(name = "BelowPlaneRestraint", from_py_object)]
#[derive(Clone)]
pub struct PyBelowPlaneRestraint {
    pub(crate) inner: BelowPlaneRestraint,
}

#[pymethods]
impl PyBelowPlaneRestraint {
    #[new]
    #[pyo3(signature = (normal, distance))]
    fn new(normal: [NpF; 3], distance: NpF) -> Self {
        Self {
            inner: BelowPlaneRestraint::new(normal, distance),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "BelowPlaneRestraint(normal={:?}, distance={})",
            self.inner.normal, self.inner.distance,
        )
    }
}

// ============================================================================
// Collective (group-level) restraints — the `with_collective_restraint` path.
//
// Mirror of the per-atom machinery above, one level up: where a [`AtomRestraint`]
// sees one atom, a [`Restraint`] sees every copy of a species at once.
// `extract_collective_restraint` tries the built-in distribution-matching
// pyclasses ([`PyGaussianPlane`], [`PyGaussianPoint`]), else duck-types
// on group-level `f`/`fg`.
// ============================================================================

// Pass-through wrapper so an owned `Arc<dyn Restraint>` itself
// implements `Restraint` and can feed
// `Target::with_collective_restraint(impl Restraint)`.
#[derive(Clone)]
pub(crate) struct SharedRestraint(pub Arc<dyn Restraint>);

impl std::fmt::Debug for SharedRestraint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SharedRestraint({})", self.0.name())
    }
}

impl Restraint for SharedRestraint {
    #[inline]
    fn f(&self, coords: &[[F; 3]], scale: F, scale2: F) -> F {
        self.0.f(coords, scale, scale2)
    }
    #[inline]
    fn fg(&self, coords: &[[F; 3]], scale: F, scale2: F, grads: &mut [[F; 3]]) -> F {
        self.0.fg(coords, scale, scale2, grads)
    }
    #[inline]
    fn is_parallel_safe(&self) -> bool {
        self.0.is_parallel_safe()
    }
    #[inline]
    fn name(&self) -> &'static str {
        self.0.name()
    }
}

pub(crate) fn extract_collective_restraint(
    obj: &Bound<'_, pyo3::types::PyAny>,
) -> PyResult<SharedRestraint> {
    if let Ok(c) = obj.extract::<PyGaussianPlane>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyGaussianPoint>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyExponentialPlane>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyExponentialPoint>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyTabulatedPlane>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }
    if let Ok(c) = obj.extract::<PyTabulatedPoint>() {
        return Ok(SharedRestraint(Arc::new(c.inner)));
    }

    // Duck-typed Python collective restraint: callable `f`/`fg` taking the
    // whole group. Bound methods resolved once, like the per-atom path.
    if let (Ok(f_method), Ok(fg_method)) = (obj.getattr("f"), obj.getattr("fg")) {
        return Ok(SharedRestraint(Arc::new(PyCallableRestraint {
            f_method: f_method.unbind(),
            fg_method: fg_method.unbind(),
        })));
    }

    Err(PyTypeError::new_err(
        "expected a restraint: a {Gaussian,Exponential,Tabulated}{Plane,Point} \
         distribution restraint, or an object with callable `f(coords, scale, scale2)` \
         and `fg(coords, scale, scale2)` methods, where `coords` is every copy's (x, y, z)",
    ))
}

// Bridge from the Rust `Restraint` trait to a Python object.
//
// Python contract:
//   obj.f(coords, scale, scale2)  -> float
//   obj.fg(coords, scale, scale2) -> (float, [(gx, gy, gz), ...])
//
// `coords` is a list of N `(x, y, z)` tuples (all copies of one species);
// `fg`'s returned gradient list must have the same length N and is accumulated
// into the caller's gradient with `+=`.
pub(crate) struct PyCallableRestraint {
    f_method: Py<pyo3::types::PyAny>,
    fg_method: Py<pyo3::types::PyAny>,
}

impl std::fmt::Debug for PyCallableRestraint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PyCallableRestraint")
    }
}

impl Restraint for PyCallableRestraint {
    fn f(&self, coords: &[[F; 3]], scale: F, scale2: F) -> F {
        Python::attach(|py| {
            let pts: Vec<(F, F, F)> = coords.iter().map(|p| (p[0], p[1], p[2])).collect();
            match self.f_method.bind(py).call1((pts, scale, scale2)) {
                Ok(res) => res.extract::<F>().unwrap_or_else(|e| {
                    stash_err(e);
                    0.0
                }),
                Err(e) => {
                    stash_err(e);
                    0.0
                }
            }
        })
    }

    fn fg(&self, coords: &[[F; 3]], scale: F, scale2: F, grads: &mut [[F; 3]]) -> F {
        Python::attach(|py| {
            let pts: Vec<(F, F, F)> = coords.iter().map(|p| (p[0], p[1], p[2])).collect();
            match self.fg_method.bind(py).call1((pts, scale, scale2)) {
                Ok(res) => match res.extract::<(F, Vec<[F; 3]>)>() {
                    Ok((v, g)) => {
                        if g.len() == grads.len() {
                            for (acc, gi) in grads.iter_mut().zip(g.iter()) {
                                acc[0] += gi[0];
                                acc[1] += gi[1];
                                acc[2] += gi[2];
                            }
                        } else {
                            stash_err(PyTypeError::new_err(format!(
                                "collective `fg` returned {} gradients for {} atoms",
                                g.len(),
                                grads.len(),
                            )));
                        }
                        v
                    }
                    Err(e) => {
                        stash_err(PyTypeError::new_err(format!(
                            "collective `fg` must return (float, [(gx, gy, gz), ...]); {e}",
                        )));
                        0.0
                    }
                },
                Err(e) => {
                    stash_err(e);
                    0.0
                }
            }
        })
    }

    fn is_parallel_safe(&self) -> bool {
        false
    }

    fn name(&self) -> &'static str {
        "PyCallableRestraint"
    }
}

/// Distribution-matching restraint for a **slab**: drives a species' signed
/// distance to a plane (`xi = x . n_hat - offset`) to a Gaussian `N(mu, sigma)`
/// via the squared 1-D Wasserstein (sorted-CDF) penalty. The compiled Rust
/// [`GaussianPlane`].
#[pyclass(name = "GaussianPlane", from_py_object)]
#[derive(Clone)]
pub struct PyGaussianPlane {
    pub(crate) inner: GaussianPlane,
}

#[pymethods]
impl PyGaussianPlane {
    /// Parameters
    /// ----------
    /// normal : (float, float, float)
    ///     Plane normal; the reaction coordinate is ``xi = x . n_hat - offset``.
    /// offset : float
    ///     Plane offset along the (normalised) normal.
    /// strength : float
    ///     Overall penalty multiplier ``lambda``.
    /// mu : float
    ///     Target Gaussian mean (Å, in ``xi``).
    /// sigma : float
    ///     Target Gaussian standard deviation (Å); must be > 0.
    #[new]
    #[pyo3(signature = (normal, offset, strength, mu, sigma))]
    fn new(normal: [NpF; 3], offset: NpF, strength: NpF, mu: NpF, sigma: NpF) -> PyResult<Self> {
        if sigma <= 0.0 {
            return Err(PyValueError::new_err("GaussianPlane sigma must be > 0"));
        }
        let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if norm <= 0.0 {
            return Err(PyValueError::new_err(
                "GaussianPlane normal must be non-zero",
            ));
        }
        Ok(Self {
            inner: GaussianPlane::new(normal, offset, strength, mu, sigma),
        })
    }

    fn __repr__(&self) -> String {
        "GaussianPlane(...)".to_string()
    }
}

/// Distribution-matching restraint for a **spherical shell**: drives a species'
/// distance to a centre (`xi = ||x - center||`) to a Gaussian `N(mu, sigma)`,
/// i.e. a shell of radius ``mu`` and thickness ``sigma`` (use ``mu`` >~ 3 ``sigma``
/// so the shell stays at positive radius). The compiled Rust [`GaussianPoint`].
#[pyclass(name = "GaussianPoint", from_py_object)]
#[derive(Clone)]
pub struct PyGaussianPoint {
    pub(crate) inner: GaussianPoint,
}

#[pymethods]
impl PyGaussianPoint {
    /// Parameters
    /// ----------
    /// center : (float, float, float)
    ///     Point the shell is centred on.
    /// strength : float
    ///     Overall penalty multiplier ``lambda``.
    /// mu : float
    ///     Target shell radius (Å).
    /// sigma : float
    ///     Target shell thickness (Å); must be > 0.
    #[new]
    #[pyo3(signature = (center, strength, mu, sigma))]
    fn new(center: [NpF; 3], strength: NpF, mu: NpF, sigma: NpF) -> PyResult<Self> {
        if sigma <= 0.0 {
            return Err(PyValueError::new_err("GaussianPoint sigma must be > 0"));
        }
        Ok(Self {
            inner: GaussianPoint::new(center, strength, mu, sigma),
        })
    }

    fn __repr__(&self) -> String {
        "GaussianPoint(...)".to_string()
    }
}

/// Distribution-matching restraint for a **diffuse layer**: drives a species'
/// signed distance to a plane (`xi = x . n_hat - offset`) to an exponential
/// distribution (density proportional to ``exp(-xi/lambda)``, ``xi >= 0``) — a
/// layer densest at the plane and decaying with length ``lambda``. The compiled
/// Rust [`ExponentialPlane`].
#[pyclass(name = "ExponentialPlane", from_py_object)]
#[derive(Clone)]
pub struct PyExponentialPlane {
    pub(crate) inner: ExponentialPlane,
}

#[pymethods]
impl PyExponentialPlane {
    /// Parameters
    /// ----------
    /// normal : (float, float, float)
    ///     Plane normal; the reaction coordinate is ``xi = x . n_hat - offset``.
    /// offset : float
    ///     Plane offset (the wall sits at ``xi = 0``).
    /// strength : float
    ///     Overall penalty multiplier.
    /// lambda_ : float
    ///     Exponential decay length (Å); must be > 0.
    #[new]
    #[pyo3(signature = (normal, offset, strength, lambda_))]
    fn new(normal: [NpF; 3], offset: NpF, strength: NpF, lambda_: NpF) -> PyResult<Self> {
        if lambda_ <= 0.0 {
            return Err(PyValueError::new_err("ExponentialPlane lambda must be > 0"));
        }
        let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if norm <= 0.0 {
            return Err(PyValueError::new_err(
                "ExponentialPlane normal must be non-zero",
            ));
        }
        Ok(Self {
            inner: ExponentialPlane::new(normal, offset, strength, lambda_),
        })
    }

    fn __repr__(&self) -> String {
        "ExponentialPlane(...)".to_string()
    }
}

/// Distribution-matching restraint for a **radial atmosphere**: drives a species'
/// distance to a centre (`xi = ||x - center||`) to an exponential distribution
/// (density proportional to ``exp(-xi/lambda)``) — densest at the centre and
/// decaying radially with length ``lambda``. The compiled Rust [`ExponentialPoint`].
#[pyclass(name = "ExponentialPoint", from_py_object)]
#[derive(Clone)]
pub struct PyExponentialPoint {
    pub(crate) inner: ExponentialPoint,
}

#[pymethods]
impl PyExponentialPoint {
    /// Parameters
    /// ----------
    /// center : (float, float, float)
    ///     Point the decay is measured from.
    /// strength : float
    ///     Overall penalty multiplier.
    /// lambda_ : float
    ///     Radial decay length (Å); must be > 0.
    #[new]
    #[pyo3(signature = (center, strength, lambda_))]
    fn new(center: [NpF; 3], strength: NpF, lambda_: NpF) -> PyResult<Self> {
        if lambda_ <= 0.0 {
            return Err(PyValueError::new_err("ExponentialPoint lambda must be > 0"));
        }
        Ok(Self {
            inner: ExponentialPoint::new(center, strength, lambda_),
        })
    }

    fn __repr__(&self) -> String {
        "ExponentialPoint(...)".to_string()
    }
}

/// Distribution-matching restraint for an **arbitrary prior** along a plane:
/// drives a species' signed distance to a plane (`xi = x . n_hat - offset`) to
/// any target density supplied as a grid ``(xs, rho)`` (e.g. a Gouy–Chapman
/// counter-ion profile sampled on a grid). The compiled Rust [`TabulatedPlane`].
#[pyclass(name = "TabulatedPlane", from_py_object)]
#[derive(Clone)]
pub struct PyTabulatedPlane {
    pub(crate) inner: TabulatedPlane,
}

#[pymethods]
impl PyTabulatedPlane {
    /// Parameters
    /// ----------
    /// normal : (float, float, float)
    ///     Plane normal; the reaction coordinate is ``xi = x . n_hat - offset``.
    /// offset : float
    ///     Plane offset.
    /// strength : float
    ///     Overall penalty multiplier.
    /// xs : list[float]
    ///     Strictly-ascending grid of ``xi`` values.
    /// rho : list[float]
    ///     Target density at each ``xs`` (>= 0, positive total mass).
    #[new]
    #[pyo3(signature = (normal, offset, strength, xs, rho))]
    fn new(
        normal: [NpF; 3],
        offset: NpF,
        strength: NpF,
        xs: Vec<NpF>,
        rho: Vec<NpF>,
    ) -> PyResult<Self> {
        validate_grid(&xs, &rho)?;
        let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if norm <= 0.0 {
            return Err(PyValueError::new_err(
                "TabulatedPlane normal must be non-zero",
            ));
        }
        Ok(Self {
            inner: TabulatedPlane::new(normal, offset, strength, &xs, &rho),
        })
    }

    fn __repr__(&self) -> String {
        "TabulatedPlane(...)".to_string()
    }
}

/// Distribution-matching restraint for an **arbitrary radial prior**: drives a
/// species' distance to a centre (`xi = ||x - center||`) to any target radial
/// density supplied as a grid ``(xs, rho)``. The compiled Rust [`TabulatedPoint`].
#[pyclass(name = "TabulatedPoint", from_py_object)]
#[derive(Clone)]
pub struct PyTabulatedPoint {
    pub(crate) inner: TabulatedPoint,
}

#[pymethods]
impl PyTabulatedPoint {
    /// Parameters
    /// ----------
    /// center : (float, float, float)
    ///     Point distances are measured from.
    /// strength : float
    ///     Overall penalty multiplier.
    /// xs : list[float]
    ///     Strictly-ascending grid of radii.
    /// rho : list[float]
    ///     Target radial density at each ``xs`` (>= 0, positive total mass).
    #[new]
    #[pyo3(signature = (center, strength, xs, rho))]
    fn new(center: [NpF; 3], strength: NpF, xs: Vec<NpF>, rho: Vec<NpF>) -> PyResult<Self> {
        validate_grid(&xs, &rho)?;
        Ok(Self {
            inner: TabulatedPoint::new(center, strength, &xs, &rho),
        })
    }

    fn __repr__(&self) -> String {
        "TabulatedPoint(...)".to_string()
    }
}

/// Validate a tabulated target grid before handing it to the Rust constructor
/// (which would otherwise panic). Mirrors `Quantile::from_grid`'s contract.
fn validate_grid(xs: &[NpF], rho: &[NpF]) -> PyResult<()> {
    if xs.len() < 2 {
        return Err(PyValueError::new_err(
            "tabulated grid needs at least 2 points",
        ));
    }
    if xs.len() != rho.len() {
        return Err(PyValueError::new_err("xs and rho must have equal length"));
    }
    if !xs.windows(2).all(|w| w[1] > w[0]) {
        return Err(PyValueError::new_err("xs must be strictly ascending"));
    }
    if rho.iter().any(|&r| r < 0.0) {
        return Err(PyValueError::new_err("rho must be >= 0"));
    }
    if rho.iter().sum::<NpF>() <= 0.0 {
        return Err(PyValueError::new_err("rho must have positive total mass"));
    }
    Ok(())
}
