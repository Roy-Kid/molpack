//! Generic marshalling between a Python frame and a Rust [`molrs::Frame`].
//!
//! molpack links `molrs-core` (Rust) but not the `molrs` Python extension — a
//! separate cdylib whose Rust types are not ABI-compatible with ours. So a
//! Python `molrs.Frame` cannot be borrowed directly; it crosses the boundary
//! through the `to_dict` / `from_dict` shape every frame flavour shares.
//!
//! This module carries **no packing logic** — it only converts data. All
//! topology assembly lives in [`molpack::assemble`] so every binding gets an
//! identical result.

use molrs::block::{Block, Column, DType};
use molrs::types::{F, I, U};
use ndarray::Array1;
use numpy::IntoPyArray;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyModule};

/// Convert any frame-like Python object — `molrs.Frame`, `molpy.Frame`, or a
/// plain `dict` (with or without a `"blocks"` envelope) — into a Rust frame.
pub fn pyframe_to_rust(frame: &Bound<'_, PyAny>) -> PyResult<molrs::Frame> {
    let blocks = blocks_dict(frame)?;
    let mut out = molrs::Frame::new();
    for (name, table) in blocks.iter() {
        let name: String = name.extract()?;
        let table = table.cast::<PyDict>()?;
        let mut block = Block::new();
        for (key, column) in table.iter() {
            let key: String = key.extract()?;
            block
                .insert_column(&key, pycolumn_to_rust(&column)?)
                .map_err(|e| PyValueError::new_err(format!("column '{key}': {e}")))?;
        }
        out.insert(name, block);
    }
    Ok(out)
}

/// Build a Python `molrs.Frame` from a Rust frame, stamping an orthorhombic
/// periodic box from `box_bounds` (`(min, max)` corners) when present.
pub fn rust_frame_to_py<'py>(
    py: Python<'py>,
    frame: &molrs::Frame,
    box_bounds: Option<([F; 3], [F; 3])>,
) -> PyResult<Bound<'py, PyAny>> {
    let molrs = PyModule::import(py, "molrs")?;

    let blocks = PyDict::new(py);
    for (name, block) in frame.iter() {
        let columns = PyDict::new(py);
        for (key, column) in block.iter() {
            columns.set_item(key, rust_column_to_py(py, column)?)?;
        }
        blocks.set_item(name, columns)?;
    }
    let data = PyDict::new(py);
    data.set_item("blocks", blocks)?;
    let frame = molrs.getattr("Frame")?.call_method1("from_dict", (data,))?;

    if let Some((lo, hi)) = box_bounds {
        let lengths = vec![hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]].into_pyarray(py);
        let kwargs = PyDict::new(py);
        kwargs.set_item("origin", lo.to_vec().into_pyarray(py))?;
        kwargs.set_item("pbc", vec![true, true, true].into_pyarray(py))?;
        let simbox = molrs
            .getattr("Box")?
            .call_method("ortho", (lengths,), Some(&kwargs))?;
        frame.setattr("box", simbox)?;
    }
    Ok(frame)
}

/// Normalize any frame flavour to a `{block: {column: array}}` mapping:
///
/// - rich `molrs.Frame` / `molpy.Frame` — via `to_dict()["blocks"]`;
/// - plain `dict` — used directly (unwrapping a `"blocks"` envelope);
/// - native `molrs` frame (e.g. from `read_pdb`) — enumerated through its
///   `keys()` / `[block]` / `block.keys()` / `block.view(col)` surface.
fn blocks_dict<'py>(frame: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyDict>> {
    if frame.hasattr("to_dict")? {
        let data = frame.call_method0("to_dict")?;
        let blocks = data
            .get_item("blocks")
            .map_err(|_| PyValueError::new_err("frame.to_dict() has no 'blocks' key"))?;
        return blocks.cast_into::<PyDict>().map_err(Into::into);
    }
    if let Ok(dict) = frame.cast::<PyDict>() {
        if dict.contains("blocks")? && !dict.contains("atoms")? {
            let blocks = dict.get_item("blocks")?.expect("checked above");
            return blocks.cast_into::<PyDict>().map_err(Into::into);
        }
        return Ok(dict.clone());
    }
    enumerate_native_frame(frame)
}

/// Read a native frame block-by-block, column-by-column into a fresh mapping.
fn enumerate_native_frame<'py>(frame: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyDict>> {
    let py = frame.py();
    let blocks = PyDict::new(py);
    for block_name in frame.call_method0("keys")?.try_iter()? {
        let block_name = block_name?;
        let block = frame.get_item(&block_name)?;
        let columns = PyDict::new(py);
        for column_name in block.call_method0("keys")?.try_iter()? {
            let column_name = column_name?;
            let array = block.call_method1("view", (&column_name,))?;
            columns.set_item(&column_name, array)?;
        }
        blocks.set_item(&block_name, columns)?;
    }
    Ok(blocks)
}

/// Convert one column (numpy array or list) to a Rust [`Column`], dispatching
/// on the numpy dtype kind; anything non-numeric becomes a string column.
fn pycolumn_to_rust(array: &Bound<'_, PyAny>) -> PyResult<Column> {
    // Plain-dict columns may be Python lists; coerce so a dtype is inferable.
    let owned;
    let array = if array.hasattr("dtype")? {
        array
    } else {
        let np = PyModule::import(array.py(), "numpy")?;
        owned = np.call_method1("asarray", (array,))?;
        &owned
    };
    let kind = array
        .getattr("dtype")
        .and_then(|d| d.getattr("kind"))
        .and_then(|k| k.extract::<String>())
        .unwrap_or_else(|_| "O".to_string());
    let column = match kind.as_str() {
        "f" => Column::from_float(Array1::from_vec(array.extract::<Vec<F>>()?).into_dyn()),
        "i" | "u" => Column::from_int(Array1::from_vec(array.extract::<Vec<I>>()?).into_dyn()),
        "b" => Column::from_bool(Array1::from_vec(array.extract::<Vec<bool>>()?).into_dyn()),
        _ => Column::from_string(Array1::from_vec(array.extract::<Vec<String>>()?).into_dyn()),
    };
    Ok(column)
}

/// Convert one Rust [`Column`] to a numpy array (string columns via `np.asarray`).
fn rust_column_to_py<'py>(py: Python<'py>, column: &Column) -> PyResult<Bound<'py, PyAny>> {
    let array = match column.dtype() {
        DType::Float => floats(column).into_pyarray(py).into_any(),
        DType::Int => ints(column).into_pyarray(py).into_any(),
        DType::Bool => bools(column).into_pyarray(py).into_any(),
        DType::UInt => uints(column).into_pyarray(py).into_any(),
        DType::U8 => u8s(column).into_pyarray(py).into_any(),
        DType::String => {
            let values = column.as_string().expect("string column").iter().cloned();
            let list = PyList::new(py, values)?;
            let np = PyModule::import(py, "numpy")?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("dtype", "str")?;
            np.call_method("asarray", (list,), Some(&kwargs))?
        }
    };
    Ok(array)
}

fn floats(c: &Column) -> Vec<F> {
    c.as_float()
        .expect("float column")
        .iter()
        .copied()
        .collect()
}
fn ints(c: &Column) -> Vec<I> {
    c.as_int().expect("int column").iter().copied().collect()
}
fn bools(c: &Column) -> Vec<bool> {
    c.as_bool().expect("bool column").iter().copied().collect()
}
fn uints(c: &Column) -> Vec<U> {
    c.as_uint().expect("uint column").iter().copied().collect()
}
fn u8s(c: &Column) -> Vec<u8> {
    c.as_u8().expect("u8 column").iter().copied().collect()
}
