//! Thread-pool introspection and control for the parallel evaluator.
//!
//! molpack's parallel path runs on rayon's process-global thread pool. That
//! pool is built **once** per process and cannot be resized afterwards — a
//! scaling study that loops over thread counts in a single interpreter only
//! ever sees the first value take effect. These helpers make that contract
//! explicit instead of silent:
//!
//! - [`rayon_enabled`] — was the wheel compiled with the `rayon` feature?
//! - [`num_threads`]   — worker count the global pool will actually use.
//! - [`init_thread_pool`] — pin the global pool size, fail-fast if already set.
//!
//! For a proper scaling sweep, set the thread count **once** per process
//! (via [`init_thread_pool`] before the first pack, or the `RAYON_NUM_THREADS`
//! environment variable) and launch one process per data point.

use pyo3::prelude::*;

/// True when the wheel was built with the `rayon` feature, i.e. the parallel
/// evaluator is actually compiled in. Used by `with_parallel_eval` to
/// fail-fast and by scaling scripts to assert before measuring.
#[inline]
pub(crate) fn rayon_compiled() -> bool {
    cfg!(feature = "rayon")
}

/// `True` if this wheel was built with parallel evaluation compiled in.
#[pyfunction]
pub(crate) fn rayon_enabled() -> bool {
    rayon_compiled()
}

/// Number of worker threads the parallel evaluator will use.
///
/// Returns the rayon global pool size when built with `rayon` (defaulting to
/// the CPU count, or `RAYON_NUM_THREADS` if set), otherwise `1` — the serial
/// build only ever runs on the calling thread.
#[pyfunction]
pub(crate) fn num_threads() -> usize {
    #[cfg(feature = "rayon")]
    {
        rayon::current_num_threads().max(1)
    }
    #[cfg(not(feature = "rayon"))]
    {
        1
    }
}

/// Pin the rayon global thread pool to `n` workers.
///
/// Must be called **before** the first parallel pack — rayon's global pool is
/// immutable once built. Raises ``RuntimeError`` if the wheel lacks the
/// `rayon` feature, if `n == 0`, or if the pool was already initialized
/// (e.g. a prior pack already triggered it). This is deliberately fail-fast:
/// silently ignoring a second resize is exactly what flattens scaling curves.
#[pyfunction]
pub(crate) fn init_thread_pool(n: usize) -> PyResult<()> {
    if !rayon_compiled() {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(
            "init_thread_pool requires a wheel built with the `rayon` feature; \
             rebuild with `maturin develop --release`",
        ));
    }
    if n == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "init_thread_pool: n must be >= 1",
        ));
    }
    #[cfg(feature = "rayon")]
    {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!(
                    "rayon global thread pool already initialized ({e}); the thread \
                     count is fixed for the life of the process — set it once before \
                     the first pack, and launch one process per scaling data point"
                ))
            })?;
    }
    Ok(())
}
