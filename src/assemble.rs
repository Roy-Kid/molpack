//! Build a topology-complete [`molrs::Frame`] from packed coordinates.
//!
//! The numeric core packs *coordinates only*; topology (bonds/angles/dihedrals/
//! impropers) and per-atom metadata ride along on each [`Target`]'s source
//! frame ([`Target::template`]). This module replays every template `count`
//! times onto the packed positions — atom columns tiled per copy, topology
//! atom-index columns offset per copy, and `id` / `mol_id` regenerated — so the
//! result is a single ready-to-use frame.
//!
//! Keeping this in the Rust core (rather than in each language binding) means
//! every binding gets identical output: bindings only marshal the frame across
//! the language boundary, never re-derive it.

use molrs::block::{Block, Column, DType};
use molrs::types::{F, I, U};
use ndarray::{Array1, ArrayD, Axis, concatenate};

use crate::target::Target;

/// Topology blocks, in a fixed order, paired with the columns that hold
/// atom indices (offset per copy). Every other column is carried verbatim.
const TOPOLOGY: [(&str, &[&str]); 4] = [
    ("bonds", &["atomi", "atomj"]),
    ("angles", &["atomi", "atomj", "atomk"]),
    ("dihedrals", &["atomi", "atomj", "atomk", "atoml"]),
    ("impropers", &["atomi", "atomj", "atomk", "atoml"]),
];
/// Atom columns overwritten from packed positions, never carried.
const COORD_COLUMNS: [&str; 3] = ["x", "y", "z"];
/// Atom columns regenerated per system, never carried.
const REGENERATED: [&str; 3] = ["id", "mol_id", "mol"];

/// Assemble a frame from each target's template replayed onto `positions`.
///
/// `positions` is the packed layout — atoms ordered target-by-target,
/// copy-by-copy, atom-by-atom (the order the packer emits). When every target
/// carries a [`Target::template`], the result has full topology; otherwise it
/// falls back to a coordinates-only `atoms` block built from element symbols.
pub fn assemble_frame(targets: &[Target], positions: &[[F; 3]]) -> molrs::Frame {
    if targets.iter().all(|t| t.template.is_some()) {
        topology_frame(targets, positions)
    } else {
        coords_only_frame(targets, positions)
    }
}

fn topology_frame(targets: &[Target], positions: &[[F; 3]]) -> molrs::Frame {
    // Union of carried atom columns across all templates (first-seen order),
    // each paired with its dtype so a template missing one can be default-filled.
    let mut carried: Vec<(String, DType)> = Vec::new();
    for target in targets {
        let atoms = template_atoms(target);
        for key in atoms.keys() {
            let carried_here = !COORD_COLUMNS.contains(&key) && !REGENERATED.contains(&key);
            if carried_here && !carried.iter().any(|(k, _)| k == key) {
                carried.push((key.to_string(), atoms.dtype(key).expect("dtype")));
            }
        }
    }

    // Per-column accumulators (one part per target, concatenated at the end).
    let mut atom_parts: Vec<Vec<Column>> = vec![Vec::new(); carried.len()];
    let mut topo_parts: Vec<Vec<(String, Vec<Column>)>> =
        TOPOLOGY.iter().map(|_| Vec::new()).collect();
    let (mut xs, mut ys, mut zs) = (Vec::new(), Vec::new(), Vec::new());
    let mut ids: Vec<I> = Vec::new();
    let mut mol_ids: Vec<I> = Vec::new();

    let mut atom_base: usize = 0;
    let mut mol_base: usize = 0;
    let mut cursor: usize = 0;

    for target in targets {
        let atoms = template_atoms(target);
        let n = atoms.nrows().unwrap_or(0);
        let count = target.count;
        let span = n * count;

        for p in &positions[cursor..cursor + span] {
            xs.push(p[0]);
            ys.push(p[1]);
            zs.push(p[2]);
        }
        cursor += span;
        ids.extend((atom_base + 1..=atom_base + span).map(|i| i as I));
        for copy in 0..count {
            mol_ids.extend(std::iter::repeat_n((mol_base + copy + 1) as I, n));
        }

        for (slot, (key, dtype)) in carried.iter().enumerate() {
            let part = match atoms.get(key) {
                Some(col) => tile_column(col, count),
                None => default_column(*dtype, span),
            };
            atom_parts[slot].push(part);
        }

        for (block_idx, (block, index_columns)) in TOPOLOGY.iter().enumerate() {
            let Some(table) = target_template(target).get(block) else {
                continue;
            };
            let rows = table.nrows().unwrap_or(0);
            for key in table.keys() {
                if key == "id" {
                    continue;
                }
                let col = table.get(key).expect("column");
                let tiled = tile_column(col, count);
                let out = if index_columns.contains(&key) {
                    offset_index_column(tiled, atom_base, n, rows)
                } else {
                    tiled
                };
                push_topo_part(&mut topo_parts[block_idx], key, out);
            }
        }

        atom_base += span;
        mol_base += count;
    }

    let mut atoms = Block::new();
    insert_int(&mut atoms, "id", ids);
    insert_int(&mut atoms, "mol_id", mol_ids);
    insert_float(&mut atoms, "x", xs);
    insert_float(&mut atoms, "y", ys);
    insert_float(&mut atoms, "z", zs);
    for ((key, _), parts) in carried.iter().zip(atom_parts) {
        atoms
            .insert_column(key, concat_columns(parts))
            .expect("atom column insert");
    }

    let mut frame = molrs::Frame::new();
    let natoms = atoms.nrows().unwrap_or(0);
    frame.insert("atoms", atoms);

    for ((block, _), columns) in TOPOLOGY.iter().zip(topo_parts) {
        if columns.is_empty() {
            continue;
        }
        let mut table = Block::new();
        for (key, parts) in columns {
            table
                .insert_column(&key, concat_columns(parts))
                .expect("topology column insert");
        }
        let nrows = table.nrows().unwrap_or(0);
        insert_int(&mut table, "id", (1..=nrows as I).collect());
        frame.insert(*block, table);
    }

    frame.meta.insert("natoms".to_string(), natoms.to_string());
    frame
}

fn coords_only_frame(targets: &[Target], positions: &[[F; 3]]) -> molrs::Frame {
    let n = positions.len();
    let mut elements: Vec<String> = Vec::with_capacity(n);
    let mut mol_ids: Vec<I> = Vec::with_capacity(n);
    let mut mol = 0usize;
    for target in targets {
        for _ in 0..target.count {
            mol += 1;
            elements.extend(target.elements.iter().cloned());
            mol_ids.extend(std::iter::repeat_n(mol as I, target.elements.len()));
        }
    }

    let mut atoms = Block::new();
    insert_int(&mut atoms, "id", (1..=n as I).collect());
    insert_float(&mut atoms, "x", positions.iter().map(|p| p[0]).collect());
    insert_float(&mut atoms, "y", positions.iter().map(|p| p[1]).collect());
    insert_float(&mut atoms, "z", positions.iter().map(|p| p[2]).collect());
    insert_int(&mut atoms, "mol_id", mol_ids);
    atoms
        .insert("element", Array1::from_vec(elements).into_dyn())
        .expect("element insert");

    let mut frame = molrs::Frame::new();
    frame.insert("atoms", atoms);
    frame.meta.insert("natoms".to_string(), n.to_string());
    frame
}

fn target_template(target: &Target) -> &molrs::Frame {
    target.template.as_ref().expect("target has a template")
}

fn template_atoms(target: &Target) -> &Block {
    target_template(target)
        .get("atoms")
        .expect("template has an 'atoms' block")
}

fn push_topo_part(parts: &mut Vec<(String, Vec<Column>)>, key: &str, value: Column) {
    if let Some((_, slot)) = parts.iter_mut().find(|(k, _)| k == key) {
        slot.push(value);
    } else {
        parts.push((key.to_string(), vec![value]));
    }
}

fn insert_float(block: &mut Block, key: &str, values: Vec<F>) {
    block
        .insert(key, Array1::from_vec(values).into_dyn())
        .expect("float column insert");
}

fn insert_int(block: &mut Block, key: &str, values: Vec<I>) {
    block
        .insert(key, Array1::from_vec(values).into_dyn())
        .expect("int column insert");
}

/// Repeat `arr` `count` times along axis 0 (numpy `tile`).
fn tile_arr<T: Clone>(arr: &ArrayD<T>, count: usize) -> ArrayD<T> {
    let views: Vec<_> = (0..count).map(|_| arr.view()).collect();
    concatenate(Axis(0), &views).expect("tile concatenate")
}

fn tile_column(col: &Column, count: usize) -> Column {
    match col.dtype() {
        DType::Float => Column::from_float(tile_arr(col.as_float().unwrap(), count)),
        DType::Int => Column::from_int(tile_arr(col.as_int().unwrap(), count)),
        DType::Bool => Column::from_bool(tile_arr(col.as_bool().unwrap(), count)),
        DType::UInt => Column::from_uint(tile_arr(col.as_uint().unwrap(), count)),
        DType::U8 => Column::from_u8(tile_arr(col.as_u8().unwrap(), count)),
        DType::String => Column::from_string(tile_arr(col.as_string().unwrap(), count)),
    }
}

/// A `len`-long default column matching `dtype`: numeric `0`, text `""`.
fn default_column(dtype: DType, len: usize) -> Column {
    match dtype {
        DType::Float => Column::from_float(ArrayD::zeros(vec![len])),
        DType::Int => Column::from_int(ArrayD::zeros(vec![len])),
        DType::Bool => Column::from_bool(ArrayD::from_elem(vec![len], false)),
        DType::UInt => Column::from_uint(ArrayD::zeros(vec![len])),
        DType::U8 => Column::from_u8(ArrayD::zeros(vec![len])),
        DType::String => Column::from_string(ArrayD::from_elem(vec![len], String::new())),
    }
}

fn concat_columns(parts: Vec<Column>) -> Column {
    let dtype = parts.first().expect("at least one part").dtype();
    match dtype {
        DType::Float => Column::from_float(concat_arr(parts.iter().map(|c| c.as_float().unwrap()))),
        DType::Int => Column::from_int(concat_arr(parts.iter().map(|c| c.as_int().unwrap()))),
        DType::Bool => Column::from_bool(concat_arr(parts.iter().map(|c| c.as_bool().unwrap()))),
        DType::UInt => Column::from_uint(concat_arr(parts.iter().map(|c| c.as_uint().unwrap()))),
        DType::U8 => Column::from_u8(concat_arr(parts.iter().map(|c| c.as_u8().unwrap()))),
        DType::String => {
            Column::from_string(concat_arr(parts.iter().map(|c| c.as_string().unwrap())))
        }
    }
}

fn concat_arr<'a, T: Clone + 'a>(arrs: impl Iterator<Item = &'a ArrayD<T>>) -> ArrayD<T> {
    let owned: Vec<&ArrayD<T>> = arrs.collect();
    let views: Vec<_> = owned.iter().map(|a| a.view()).collect();
    concatenate(Axis(0), &views).expect("concatenate columns")
}

/// Add `atom_base + copy * n` to each copy's segment of a tiled index column,
/// preserving its integer dtype (readers emit `int` or `uint` atom indices).
/// A non-integer index column is left unchanged.
fn offset_index_column(tiled: Column, atom_base: usize, n: usize, rows: usize) -> Column {
    // Element `i` of the tiled `[copy0, copy1, …]` layout belongs to copy `i / rows`.
    let offset_at = |i: usize| atom_base + (i / rows) * n;
    match tiled.dtype() {
        DType::Int => {
            let mut arr = tiled.as_int().expect("int column").clone();
            for (i, v) in arr.iter_mut().enumerate() {
                *v += offset_at(i) as I;
            }
            Column::from_int(arr)
        }
        DType::UInt => {
            let mut arr = tiled.as_uint().expect("uint column").clone();
            for (i, v) in arr.iter_mut().enumerate() {
                *v += offset_at(i) as U;
            }
            Column::from_uint(arr)
        }
        DType::U8 => {
            let mut arr = tiled.as_u8().expect("u8 column").clone();
            for (i, v) in arr.iter_mut().enumerate() {
                *v += offset_at(i) as u8;
            }
            Column::from_u8(arr)
        }
        _ => tiled,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn col_int(frame: &molrs::Frame, block: &str, key: &str) -> Vec<I> {
        frame
            .get(block)
            .unwrap()
            .get_int(key)
            .unwrap()
            .iter()
            .copied()
            .collect()
    }

    fn col_str(frame: &molrs::Frame, block: &str, key: &str) -> Vec<String> {
        frame
            .get(block)
            .unwrap()
            .get_string(key)
            .unwrap()
            .iter()
            .cloned()
            .collect()
    }

    fn diatomic() -> molrs::Frame {
        let mut atoms = Block::new();
        atoms
            .insert(
                "element",
                Array1::from_vec(vec!["C".to_string(), "H".to_string()]).into_dyn(),
            )
            .unwrap();
        atoms
            .insert(
                "type",
                Array1::from_vec(vec!["A".to_string(), "B".to_string()]).into_dyn(),
            )
            .unwrap();
        for c in COORD_COLUMNS {
            atoms
                .insert(c, Array1::from_vec(vec![0.0 as F, 1.0 as F]).into_dyn())
                .unwrap();
        }
        let mut bonds = Block::new();
        bonds
            .insert("atomi", Array1::from_vec(vec![0 as I]).into_dyn())
            .unwrap();
        bonds
            .insert("atomj", Array1::from_vec(vec![1 as I]).into_dyn())
            .unwrap();
        let mut frame = molrs::Frame::new();
        frame.insert("atoms", atoms);
        frame.insert("bonds", bonds);
        frame
    }

    fn argon() -> molrs::Frame {
        let mut atoms = Block::new();
        atoms
            .insert(
                "element",
                Array1::from_vec(vec!["Ar".to_string()]).into_dyn(),
            )
            .unwrap();
        for c in COORD_COLUMNS {
            atoms
                .insert(c, Array1::from_vec(vec![0.0 as F]).into_dyn())
                .unwrap();
        }
        let mut frame = molrs::Frame::new();
        frame.insert("atoms", atoms);
        frame
    }

    fn positions(n: usize) -> Vec<[F; 3]> {
        (0..n).map(|i| [i as F, 0.0, (2 * i) as F]).collect()
    }

    #[test]
    fn topology_indices_offset_per_copy() {
        let target = Target::new(diatomic(), 3);
        let frame = assemble_frame(&[target], &positions(6));

        assert_eq!(col_int(&frame, "bonds", "atomi"), [0, 2, 4]);
        assert_eq!(col_int(&frame, "bonds", "atomj"), [1, 3, 5]);
        assert_eq!(col_int(&frame, "bonds", "id"), [1, 2, 3]);
        assert_eq!(col_int(&frame, "atoms", "id"), [1, 2, 3, 4, 5, 6]);
        assert_eq!(col_int(&frame, "atoms", "mol_id"), [1, 1, 2, 2, 3, 3]);
        assert_eq!(
            col_str(&frame, "atoms", "type"),
            ["A", "B", "A", "B", "A", "B"]
        );
    }

    #[test]
    fn coords_come_from_positions() {
        let target = Target::new(diatomic(), 2);
        let frame = assemble_frame(&[target], &positions(4));
        let xs: Vec<F> = frame
            .get("atoms")
            .unwrap()
            .get_float("x")
            .unwrap()
            .iter()
            .copied()
            .collect();
        assert_eq!(xs, [0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn heterogeneous_schema_is_union_filled() {
        // diatomic has a `type` column; argon does not — argon's atoms get "".
        let targets = [Target::new(diatomic(), 2), Target::new(argon(), 2)];
        let frame = assemble_frame(&targets, &positions(6));

        assert_eq!(
            col_str(&frame, "atoms", "type"),
            ["A", "B", "A", "B", "", ""]
        );
        assert_eq!(col_int(&frame, "atoms", "mol_id"), [1, 1, 2, 2, 3, 4]);
        // Bonds belong only to the two diatomics.
        assert_eq!(col_int(&frame, "bonds", "atomi"), [0, 2]);
    }

    #[test]
    fn no_template_falls_back_to_coords_only() {
        let target = Target::from_coords(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &[1.5, 1.5], 2);
        let frame = assemble_frame(&[target], &positions(4));

        assert_eq!(col_int(&frame, "atoms", "id"), [1, 2, 3, 4]);
        assert!(frame.get("bonds").is_none());
    }
}
