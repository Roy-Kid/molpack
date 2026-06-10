//! Build a topology-complete `molrs.Frame` from packed coordinates + templates.
//!
//! molpack's numeric core packs *coordinates only* — bonds/angles/dihedrals/
//! impropers and per-atom metadata (`id`, `type`, `charge`, `mass`, …) are
//! dropped at the Rust boundary because the optimizer never needs them. To let
//! [`crate::packer::PyPackResult::frame`] hand back a ready-to-write frame, the
//! wheel keeps a reference to each input template and replays it here: every
//! template is replicated `count` times onto the packed positions, topology
//! indices are offset per copy, and `id` / `mol_id` are regenerated.
//!
//! The result is a genuine `molrs.Frame`: molpack links `molrs-core` (Rust) but
//! not the `molrs-python` bindings (a cdylib), so the only way to produce the
//! Python type is to drive the user's installed `molrs` module at runtime —
//! the same `py.import("molrs")` pattern the `.inp` script loader uses. Column
//! math goes through `numpy`; molpack just orchestrates. Force fields are out
//! of scope — callers merge those separately.

use molpack::F;
use numpy::IntoPyArray;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyModule};

/// A block's columns as `(name, array)` pairs, in declaration order.
type NamedColumns<'py> = Vec<(String, Bound<'py, PyAny>)>;

const COORD_COLS: [&str; 3] = ["x", "y", "z"];
/// Regenerated from scratch, never carried from the template.
const REGENERATED: [&str; 3] = ["id", "mol_id", "mol"];
const TOPO_BLOCKS: [&str; 4] = ["bonds", "angles", "dihedrals", "impropers"];

/// 0-based atom-index columns per topology block (offset per copy); every other
/// column is carried verbatim.
fn topo_index_cols(block: &str) -> &'static [&'static str] {
    match block {
        "bonds" => &["atomi", "atomj"],
        "angles" => &["atomi", "atomj", "atomk"],
        "dihedrals" | "impropers" => &["atomi", "atomj", "atomk", "atoml"],
        _ => &[],
    }
}

/// Read a block's column as a numpy-ready object, duck-typed over the frame
/// flavours molpack accepts: `molpy.Block` / `dict` (via `[]`, which also
/// exposes object columns like `element`) then `molrs.Block` (via `.view`).
fn get_column<'py>(block: &Bound<'py, PyAny>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    if let Ok(col) = block.get_item(name) {
        return Ok(col);
    }
    block.call_method1("view", (name,))
}

/// Enumerate `(name, column)` for every column of a block, in declaration order.
fn read_block_columns<'py>(block: &Bound<'py, PyAny>) -> PyResult<NamedColumns<'py>> {
    let keys = block.call_method0("keys")?;
    let mut cols = Vec::new();
    for key in keys.try_iter()? {
        let name: String = key?.extract()?;
        let col = get_column(block, &name)?;
        cols.push((name, col));
    }
    Ok(cols)
}

/// Append `value` to the per-column accumulator list `dict[key]`, creating it on
/// first use so column order follows first insertion.
fn push<'py>(
    py: Python<'py>,
    dict: &Bound<'py, PyDict>,
    key: &str,
    value: Bound<'py, PyAny>,
) -> PyResult<()> {
    match dict.get_item(key)? {
        Some(existing) => existing.cast_into::<PyList>()?.append(value)?,
        None => {
            let list = PyList::empty(py);
            list.append(value)?;
            dict.set_item(key, list)?;
        }
    }
    Ok(())
}

/// Concatenate an accumulator list of arrays and insert it into `block`,
/// coercing object-dtype columns (e.g. molpy's `element`) to unicode so
/// `molrs.Block.insert` accepts them.
fn insert_concat<'py>(
    np: &Bound<'py, PyModule>,
    block: &Bound<'py, PyAny>,
    name: &str,
    parts: &Bound<'py, PyAny>,
) -> PyResult<usize> {
    let mut arr = np.call_method1("concatenate", (parts,))?;
    let kind: String = arr.getattr("dtype")?.getattr("kind")?.extract()?;
    if kind == "O" {
        arr = arr.call_method1("astype", ("str",))?;
    }
    let n: usize = arr.len()?;
    block.call_method1("insert", (name, arr))?;
    Ok(n)
}

/// Stamp an orthorhombic periodic box on `frame` from `(min, max)` corners.
fn set_box<'py>(
    py: Python<'py>,
    molrs: &Bound<'py, PyModule>,
    frame: &Bound<'py, PyAny>,
    bounds: Option<([F; 3], [F; 3])>,
) -> PyResult<()> {
    let Some((lo, hi)) = bounds else {
        return Ok(());
    };
    let lengths = vec![hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]].into_pyarray(py);
    let origin = lo.to_vec().into_pyarray(py);
    let kwargs = PyDict::new(py);
    kwargs.set_item("origin", origin)?;
    kwargs.set_item("pbc", vec![true, true, true].into_pyarray(py))?;
    let simbox = molrs
        .getattr("Box")?
        .call_method("ortho", (lengths,), Some(&kwargs))?;
    frame.setattr("box", simbox)?;
    Ok(())
}

/// Build a coordinates-only `molrs.Frame` (single `atoms` block with `x`, `y`,
/// `z`, `element`) — the fallback when no source templates are retained
/// (`.inp` script packing).
pub fn coords_only_frame<'py>(
    py: Python<'py>,
    xs: Vec<F>,
    ys: Vec<F>,
    zs: Vec<F>,
    elements: Vec<String>,
    bounds: Option<([F; 3], [F; 3])>,
) -> PyResult<Bound<'py, PyAny>> {
    let molrs = PyModule::import(py, "molrs")?;
    let frame = molrs.getattr("Frame")?.call0()?;
    let block = molrs.getattr("Block")?.call0()?;
    block.call_method1("insert", ("x", xs.into_pyarray(py)))?;
    block.call_method1("insert", ("y", ys.into_pyarray(py)))?;
    block.call_method1("insert", ("z", zs.into_pyarray(py)))?;
    block.call_method1("insert", ("element", PyList::new(py, elements)?))?;
    frame.set_item("atoms", block)?;
    set_box(py, &molrs, &frame, bounds)?;
    Ok(frame)
}

/// Replay each template `count` times onto the packed `positions` and return a
/// boxed `molrs.Frame` with full, re-indexed topology.
pub fn topology_frame<'py>(
    py: Python<'py>,
    positions: &[[F; 3]],
    templates: &[(Py<PyAny>, usize)],
    bounds: Option<([F; 3], [F; 3])>,
) -> PyResult<Bound<'py, PyAny>> {
    let molrs = PyModule::import(py, "molrs")?;
    let np = PyModule::import(py, "numpy")?;

    // Per-column accumulators: atoms is name -> [arrays]; topo is block ->
    // (name -> [arrays]). Python dicts preserve first-insertion order.
    let atom_cols = PyDict::new(py);
    let topo = PyDict::new(py);
    let mut atom_base: usize = 0;
    let mut mol_id: usize = 0;
    let mut cursor: usize = 0;

    for (template, count) in templates {
        let template = template.bind(py);
        let atoms = template.get_item("atoms")?;
        let cols = read_block_columns(&atoms)?;
        let n: usize = cols.first().map(|(_, c)| c.len()).transpose()?.unwrap_or(0);

        // Columns replayed verbatim each copy: all except coords (overwritten)
        // and id/mol_id (regenerated).
        let carry: Vec<&(String, Bound<'py, PyAny>)> = cols
            .iter()
            .filter(|(k, _)| {
                !COORD_COLS.contains(&k.as_str()) && !REGENERATED.contains(&k.as_str())
            })
            .collect();

        // Topology blocks present on this template, pre-read once.
        let mut topo_blocks: Vec<(&str, NamedColumns<'py>)> = Vec::new();
        for &block in &TOPO_BLOCKS {
            if let Ok(b) = template.get_item(block) {
                if !b.is_none() {
                    topo_blocks.push((block, read_block_columns(&b)?));
                }
            }
        }

        for _ in 0..*count {
            mol_id += 1;
            let end = cursor + n;
            if end > positions.len() {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "packed positions hold {} atoms but templates need at least {end}; \
                     template/count order must match pack()",
                    positions.len()
                )));
            }

            let ids: Vec<i64> = ((atom_base as i64 + 1)..=(atom_base as i64 + n as i64)).collect();
            push(py, &atom_cols, "id", ids.into_pyarray(py).into_any())?;
            let mols: Vec<i64> = vec![mol_id as i64; n];
            push(py, &atom_cols, "mol_id", mols.into_pyarray(py).into_any())?;
            for (name, col) in &carry {
                push(py, &atom_cols, name, col.clone())?;
            }
            let xs: Vec<F> = positions[cursor..end].iter().map(|p| p[0]).collect();
            let ys: Vec<F> = positions[cursor..end].iter().map(|p| p[1]).collect();
            let zs: Vec<F> = positions[cursor..end].iter().map(|p| p[2]).collect();
            cursor = end;
            push(py, &atom_cols, "x", xs.into_pyarray(py).into_any())?;
            push(py, &atom_cols, "y", ys.into_pyarray(py).into_any())?;
            push(py, &atom_cols, "z", zs.into_pyarray(py).into_any())?;

            for (block, bcols) in &topo_blocks {
                let idx_cols = topo_index_cols(block);
                let dst = match topo.get_item(block)? {
                    Some(d) => d.cast_into::<PyDict>()?,
                    None => {
                        let d = PyDict::new(py);
                        topo.set_item(block, &d)?;
                        d
                    }
                };
                for (name, val) in bcols {
                    if name == "id" {
                        continue;
                    }
                    let out = if idx_cols.contains(&name.as_str()) {
                        val.call_method1("__add__", (atom_base as i64,))?
                    } else {
                        val.clone()
                    };
                    push(py, &dst, name, out)?;
                }
            }
            atom_base += n;
        }
    }

    if cursor != positions.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "packed positions hold {} atoms but templates account for {cursor}; \
             template/count order must match pack()",
            positions.len()
        )));
    }

    let frame = molrs.getattr("Frame")?.call0()?;
    let atoms_block = molrs.getattr("Block")?.call0()?;
    for (name, parts) in atom_cols.iter() {
        let name: String = name.extract()?;
        insert_concat(&np, &atoms_block, &name, &parts)?;
    }
    frame.set_item("atoms", atoms_block)?;

    for (block, cols) in topo.iter() {
        let block: String = block.extract()?;
        let cols = cols.cast_into::<PyDict>()?;
        let topo_block = molrs.getattr("Block")?.call0()?;
        let mut nrows = 0usize;
        for (name, parts) in cols.iter() {
            let name: String = name.extract()?;
            nrows = insert_concat(&np, &topo_block, &name, &parts)?;
        }
        let ids: Vec<i64> = (1..=nrows as i64).collect();
        topo_block.call_method1("insert", ("id", ids.into_pyarray(py)))?;
        frame.set_item(&block, topo_block)?;
    }

    set_box(py, &molrs, &frame, bounds)?;
    Ok(frame)
}
