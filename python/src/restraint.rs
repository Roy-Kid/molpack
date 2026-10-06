//! Python wrappers for molecular packing restraints.
//!
//! Geometry is a molrs region object (``molrs.spatial.Sphere``, ``Cuboid``,
//! ``Parallelepiped``, ``HalfSpace``, ``Cylinder``, ``Ellipsoid``,
//! ``Polyhedron``, ``SphereUnion``, or any ``&`` / ``|`` / ``~`` composition);
//! it reaches this wheel as a ``molrs.RegionRef/<line>`` capsule and is lifted
//! through [`RegionRestraint`]. molpack defines no geometric restraint class of
//! its own. Custom Python-defined restraints are supported via **duck typing**:
//! any object exposing callable ``f(coords, scale, scale2)`` and
//! ``fg(coords, scale, scale2)`` attributes may be passed to
//! ``Target.with_restraint``; see [`PyCallableRestraint`] for the group contract.

use std::sync::Arc;

use crate::errors::stash_err;
use crate::interop::region_from_py;
use molpack::RegionRestraint;
use molpack::{
    AtomRestraint, ExponentialPlane, ExponentialPoint, GaussianPlane, GaussianPoint, GroupCtx,
    Restraint, SelfSeparation, TabulatedPlane, TabulatedPoint,
};
use molrs::op::types::F;
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
    fn holds_along(&self, shift: [F; 3]) -> bool {
        self.0.holds_along(shift)
    }
    #[inline]
    fn declared_cell(&self) -> Option<molrs::spatial::SimBox> {
        self.0.declared_cell()
    }
}

// ============================================================================
// Extractor: a molrs region (by capsule), else duck-type on `f`/`fg`.
// ============================================================================

/// Whether `obj` is a molrs region object: it exports the region capsule.
fn is_region(obj: &Bound<'_, pyo3::types::PyAny>) -> PyResult<bool> {
    obj.hasattr("_ffi_regionref_capsule")
}

pub(crate) fn extract_restraint(
    obj: &Bound<'_, pyo3::types::PyAny>,
) -> PyResult<SharedAtomRestraint> {
    if let Some(region) = try_region(obj)? {
        return Ok(region);
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
        "expected a restraint: a molrs region (any object exposing \
         `_ffi_regionref_capsule()` — molrs.spatial.Sphere / Cuboid / Parallelepiped / HalfSpace / \
         Cylinder / Ellipsoid / Polyhedron / SphereUnion or a `&` / `|` / `~` composition), \
         or an object with callable `f(x, scale, scale2)` and `fg(x, scale, scale2)` methods",
    ))
}

/// The per-atom region lift of `obj` when it is a molrs region, `None`
/// otherwise (no duck typing).
///
/// The unified [`crate::target::PyTarget::with_restraint`] entry point uses
/// this to route a region to the per-atom path and everything else (built-in
/// distribution restraints + duck-typed objects) to the group path. A region
/// from another molrs minor line is an error here, never a fall-through.
pub(crate) fn try_region(
    obj: &Bound<'_, pyo3::types::PyAny>,
) -> PyResult<Option<SharedAtomRestraint>> {
    if !is_region(obj)? {
        return Ok(None);
    }
    let region = region_from_py(obj)?;
    Ok(Some(SharedAtomRestraint(Arc::new(RegionRestraint(region)))))
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
    fn f(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>) -> F {
        self.0.f(coords, ctx)
    }
    #[inline]
    fn fg(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>, grads: &mut [[F; 3]]) -> F {
        self.0.fg(coords, ctx, grads)
    }
    #[inline]
    fn is_bound(&self) -> bool {
        self.0.is_bound()
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
    if let Ok(c) = obj.extract::<PySelfSeparation>() {
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
        "expected a restraint: a molrs region (any object exposing \
         `_ffi_regionref_capsule()`), a {Gaussian,Exponential,Tabulated}{Plane,Point} \
         distribution restraint, a SelfSeparation restraint, or an object with callable \
         `f(coords, scale, scale2)` and `fg(coords, scale, scale2)` methods, where \
         `coords` is every copy's (x, y, z)",
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
    fn f(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>) -> F {
        Python::attach(|py| {
            let pts: Vec<(F, F, F)> = coords.iter().map(|p| (p[0], p[1], p[2])).collect();
            match self.f_method.bind(py).call1((pts, ctx.scale, ctx.scale2)) {
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

    fn fg(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>, grads: &mut [[F; 3]]) -> F {
        Python::attach(|py| {
            let pts: Vec<(F, F, F)> = coords.iter().map(|p| (p[0], p[1], p[2])).collect();
            match self.fg_method.bind(py).call1((pts, ctx.scale, ctx.scale2)) {
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
    fn new(normal: [F; 3], offset: F, strength: F, mu: F, sigma: F) -> PyResult<Self> {
        if sigma <= 0.0 {
            return Err(PyValueError::new_err("GaussianPlane sigma must be > 0"));
        }
        if molrs::op::vec3::normalize(normal).is_none() {
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
    fn new(center: [F; 3], strength: F, mu: F, sigma: F) -> PyResult<Self> {
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
    fn new(normal: [F; 3], offset: F, strength: F, lambda_: F) -> PyResult<Self> {
        if lambda_ <= 0.0 {
            return Err(PyValueError::new_err("ExponentialPlane lambda must be > 0"));
        }
        if molrs::op::vec3::normalize(normal).is_none() {
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
    fn new(center: [F; 3], strength: F, lambda_: F) -> PyResult<Self> {
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
    fn new(normal: [F; 3], offset: F, strength: F, xs: Vec<F>, rho: Vec<F>) -> PyResult<Self> {
        validate_grid(&xs, &rho)?;
        if molrs::op::vec3::normalize(normal).is_none() {
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
    fn new(center: [F; 3], strength: F, xs: Vec<F>, rho: Vec<F>) -> PyResult<Self> {
        validate_grid(&xs, &rho)?;
        Ok(Self {
            inner: TabulatedPoint::new(center, strength, &xs, &rho),
        })
    }

    fn __repr__(&self) -> String {
        "TabulatedPoint(...)".to_string()
    }
}

/// Keep every pair of copies of one species at least ``d_min`` apart, measured
/// **centre to centre** — the anti-clustering restraint. The compiled Rust
/// [`SelfSeparation`].
///
/// The packer's own pair term only stops molecules overlapping; nothing in it
/// distinguishes two copies of one species from a copy of each of two, so a
/// species may pile its copies into one corner. This states the missing bound.
///
/// The penalty is silent above ``d_min`` and grows as ``(d_min - D)^2`` below
/// it, where ``D`` is the minimum-image centre-to-centre distance. Quadratic in
/// the length by which the bound is missed, like a geometric restraint's
/// penalty, so it is commensurate with the rest of ``frest``.
#[pyclass(name = "SelfSeparation", from_py_object)]
#[derive(Clone)]
pub struct PySelfSeparation {
    pub(crate) inner: SelfSeparation,
}

#[pymethods]
impl PySelfSeparation {
    /// Parameters
    /// ----------
    /// d_min : float
    ///     Minimum centre-to-centre distance between two copies (Å); must be > 0.
    /// strength : float, default 1.0
    ///     Overall penalty multiplier ``lambda``; ``1.0`` weights a shortfall
    ///     like a geometric restraint weights an equal overshoot, below ``1.0``
    ///     makes the bound softer. Must be > 0.
    #[new]
    #[pyo3(signature = (d_min, strength = 1.0))]
    fn new(d_min: F, strength: F) -> PyResult<Self> {
        if d_min <= 0.0 {
            return Err(PyValueError::new_err("SelfSeparation d_min must be > 0"));
        }
        if strength <= 0.0 {
            return Err(PyValueError::new_err("SelfSeparation strength must be > 0"));
        }
        Ok(Self {
            inner: SelfSeparation::new(d_min, strength),
        })
    }

    /// The minimum centre-to-centre distance this restraint asks for.
    #[getter]
    fn d_min(&self) -> F {
        self.inner.d_min()
    }

    fn __repr__(&self) -> String {
        format!("SelfSeparation(d_min={})", self.inner.d_min())
    }
}

/// Validate a tabulated target grid before handing it to the Rust constructor
/// (which would otherwise panic). Mirrors `Quantile::from_grid`'s contract.
fn validate_grid(xs: &[F], rho: &[F]) -> PyResult<()> {
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
    if rho.iter().sum::<F>() <= 0.0 {
        return Err(PyValueError::new_err("rho must have positive total mass"));
    }
    Ok(())
}
