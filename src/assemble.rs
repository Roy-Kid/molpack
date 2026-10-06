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

use molrs::store::block::{Block, Column};
use molrs::store::keys;
use molrs::store::schema::block_names::ATOMS;
use molrs::store::schema::{self, RowKind, RowReference, relation_endpoints};
use molrs::types::{F, Idx};
use ndarray::{Array1, Array2};

use crate::error::PackError;
use crate::target::Target;

/// Atom columns never carried from a template: coordinates and image flags
/// are replaced by the packed positions, and `id` / `mol_id` / `mol` are
/// regenerated per system.
const NOT_CARRIED: [&str; 9] = [
    keys::X,
    keys::Y,
    keys::Z,
    keys::IX,
    keys::IY,
    keys::IZ,
    keys::ID,
    keys::MOL_ID,
    "mol",
];

/// Assemble a frame from each target's template replayed onto `positions`.
///
/// `positions` is the packed layout — atoms ordered target-by-target,
/// copy-by-copy, atom-by-atom (the order the packer emits). When every target
/// carries a [`Target::template`], the result has full topology; otherwise it
/// falls back to a coordinates-only `atoms` block built from element symbols.
///
/// Topology is every canonical relation block of the molrs schema (bonds,
/// angles, dihedrals, impropers, pairs, exclusions) a template carries, in
/// the order the templates first name them. Other blocks are not replayed.
///
/// # Errors
/// [`PackError::TemplateColumns`] when two templates carry one column under
/// different dtypes, so the copies cannot share one column, or when a
/// template's relation block lacks a 1-D `UInt` endpoint column, so its
/// copies could not be offset.
pub fn assemble_frame(targets: &[Target], positions: &[[F; 3]]) -> Result<molrs::Frame, PackError> {
    if targets.iter().all(|t| t.template.is_some()) {
        let counts: Vec<usize> = targets.iter().map(|t| t.count).collect();
        topology_frame(targets, &counts, positions)
    } else {
        Ok(coords_only_frame(targets, positions))
    }
}

/// Check up front that the templates can be assembled into one frame, so a
/// dtype clash is named before a run rather than after it. One copy of each
/// target is enough to meet every column.
pub(crate) fn check_templates(targets: &[Target]) -> Result<(), PackError> {
    if targets.iter().all(|t| t.template.is_some()) {
        let ones = vec![1; targets.len()];
        let origin: usize = targets.iter().map(Target::natoms).sum();
        topology_frame(targets, &ones, &vec![[0.0; 3]; origin])?;
    }
    Ok(())
}

/// The topology-complete frame, with `counts[i]` copies of `targets[i]`.
///
/// Each target's replayed template is copied with [`molrs::Frame::replicate`]
/// (endpoints offset per copy), its endpoints shifted past the atoms of the
/// targets before it, and the per-target parts joined with [`Block::stack`] —
/// a column one template lacks is null on the other templates' rows.
fn topology_frame(
    targets: &[Target],
    counts: &[usize],
    positions: &[[F; 3]],
) -> Result<molrs::Frame, PackError> {
    let mut atom_parts: Vec<Block> = Vec::with_capacity(targets.len());
    let mut relation_parts: Vec<(String, Vec<Block>)> = Vec::new();
    let mut ids: Vec<Idx> = Vec::new();
    let mut mol_ids: Vec<Idx> = Vec::new();

    let mut atom_base: usize = 0;
    let mut mol_base: usize = 0;
    for (target, &count) in targets.iter().zip(counts) {
        let one = replayed(target_template(target))?;
        let n = one.get(ATOMS).and_then(Block::nrows).unwrap_or(0);
        let span = n * count;

        ids.extend((atom_base + 1..=atom_base + span).map(|i| i as Idx));
        for copy in 0..count {
            mol_ids.extend(std::iter::repeat_n((mol_base + copy + 1) as Idx, n));
        }

        let (copies, _, _) = one.replicate(count).map_err(column_error)?.into_inner();
        let bases = row_bases(&relation_parts);
        for (name, mut block) in copies {
            if name == ATOMS {
                atom_parts.push(block);
                continue;
            }
            // `replicate` offset each copy within this target; shift every
            // local row reference past the rows earlier targets put in the
            // block it indexes.
            let declared: Vec<(&str, &str)> = block.targets().collect();
            let refs = relation_endpoints(&name, |k| block.contains_key(k), &declared);
            for r in refs.into_iter().filter(RowReference::is_local) {
                let base = if r.target == ATOMS {
                    atom_base
                } else {
                    bases
                        .iter()
                        .find(|(k, _)| *k == r.target)
                        .map_or(0, |&(_, rows)| rows)
                };
                if let Some(index) = block.get_mut(&r.column).and_then(Column::as_uint_mut) {
                    *index += base as Idx;
                }
            }
            match relation_parts.iter_mut().find(|(k, _)| *k == name) {
                Some((_, parts)) => parts.push(block),
                None => relation_parts.push((name, vec![block])),
            }
        }

        atom_base += span;
        mol_base += count;
    }

    let mut atoms = Block::stack(&atom_parts).map_err(column_error)?;
    insert_front(&mut atoms, keys::ID, 0, Array1::from_vec(ids))?;
    insert_front(&mut atoms, keys::MOL_ID, 1, Array1::from_vec(mol_ids))?;
    atoms
        .set_coords(xyz(positions).view())
        .map_err(column_error)?;
    for (slot, key) in keys::COORDS.into_iter().enumerate() {
        atoms.move_column(key, 2 + slot).map_err(column_error)?;
    }

    let mut frame = molrs::Frame::new();
    frame.insert(ATOMS, atoms);
    for (name, parts) in relation_parts {
        let mut table = Block::stack(&parts).map_err(column_error)?;
        let nrows = table.nrows().unwrap_or(0);
        insert_front(&mut table, keys::ID, 0, (1..=nrows as Idx).collect())?;
        frame.insert(name, table);
    }
    Ok(frame)
}

/// One copy of what a template contributes: its carried atom columns (the row
/// count kept even when nothing is carried) and every canonical relation block
/// of the molrs schema it has, without the per-row `id` that is regenerated.
fn replayed(template: &molrs::Frame) -> Result<molrs::Frame, PackError> {
    let atoms = template.get(ATOMS).expect("template has an 'atoms' block");
    let carried: Vec<&str> = atoms.keys().filter(|k| !NOT_CARRIED.contains(k)).collect();
    let carried = if carried.is_empty() {
        let mut rows = Block::new();
        rows.resize(atoms.nrows().unwrap_or(0))
            .map_err(column_error)?;
        rows
    } else {
        atoms.select_columns(&carried).map_err(column_error)?
    };
    let mut one = molrs::Frame::new();
    one.insert(ATOMS, carried);
    for (name, table) in template.iter() {
        let is_relation = schema::block(name)
            .is_some_and(|spec| matches!(spec.row_kind, RowKind::Relation { .. }));
        if is_relation {
            let mut table = table.clone();
            table.remove(keys::ID);
            one.insert(name, table);
        }
    }
    Ok(one)
}

fn coords_only_frame(targets: &[Target], positions: &[[F; 3]]) -> molrs::Frame {
    let n = positions.len();
    let mut elements: Vec<String> = Vec::with_capacity(n);
    let mut mol_ids: Vec<Idx> = Vec::with_capacity(n);
    let mut mol = 0usize;
    for target in targets {
        for _ in 0..target.count {
            mol += 1;
            elements.extend(target.elements.iter().cloned());
            mol_ids.extend(std::iter::repeat_n(mol as Idx, target.elements.len()));
        }
    }

    let mut atoms = Block::new();
    let inserted = atoms
        .insert(keys::ID, (1..=n as Idx).collect::<Array1<Idx>>().into_dyn())
        .and_then(|()| atoms.set_coords(xyz(positions).view()))
        .and_then(|()| atoms.insert(keys::MOL_ID, Array1::from_vec(mol_ids).into_dyn()))
        .and_then(|()| atoms.insert(keys::ELEMENT, Array1::from_vec(elements).into_dyn()));
    inserted.expect("canonical columns of one length always insert");

    let mut frame = molrs::Frame::new();
    frame.insert(ATOMS, atoms);
    frame
}

/// Rows each relation block holds so far: the base a later target's
/// references into that block are shifted by.
fn row_bases(parts: &[(String, Vec<Block>)]) -> Vec<(String, usize)> {
    parts
        .iter()
        .map(|(name, blocks)| {
            let rows = blocks.iter().map(|b| b.nrows().unwrap_or(0)).sum();
            (name.clone(), rows)
        })
        .collect()
}

/// `positions` as the `N × 3` array [`Block::set_coords`] takes.
fn xyz(positions: &[[F; 3]]) -> Array2<F> {
    Array2::from(positions.to_vec())
}

fn target_template(target: &Target) -> &molrs::Frame {
    target.template.as_ref().expect("target has a template")
}

/// Insert a generated column and move it to `index`.
fn insert_front<T: molrs::store::block::BlockDtype>(
    block: &mut Block,
    key: &str,
    index: usize,
    values: Array1<T>,
) -> Result<(), PackError> {
    block.insert(key, values.into_dyn()).map_err(column_error)?;
    block.move_column(key, index).map_err(column_error)
}

fn column_error(err: impl std::fmt::Display) -> PackError {
    PackError::TemplateColumns {
        detail: err.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molrs::store::block::Column;
    use ndarray::ArrayD;

    fn col_uint(frame: &molrs::Frame, block: &str, key: &str) -> Vec<Idx> {
        frame
            .get(block)
            .unwrap()
            .get(key)
            .and_then(Column::as_uint)
            .unwrap()
            .iter()
            .copied()
            .collect()
    }

    fn col_str(frame: &molrs::Frame, block: &str, key: &str) -> Vec<String> {
        frame
            .get(block)
            .unwrap()
            .get(key)
            .and_then(Column::as_string)
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
        for c in keys::COORDS {
            atoms
                .insert(c, Array1::from_vec(vec![0.0 as F, 1.0 as F]).into_dyn())
                .unwrap();
        }
        let mut bonds = Block::new();
        bonds
            .insert("atomi", Array1::from_vec(vec![0 as Idx]).into_dyn())
            .unwrap();
        bonds
            .insert("atomj", Array1::from_vec(vec![1 as Idx]).into_dyn())
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
        for c in keys::COORDS {
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
        let frame = assemble_frame(&[target], &positions(6)).unwrap();

        assert_eq!(col_uint(&frame, "bonds", "atomi"), [0, 2, 4]);
        assert_eq!(col_uint(&frame, "bonds", "atomj"), [1, 3, 5]);
        assert_eq!(col_uint(&frame, "bonds", "id"), [1, 2, 3]);
        assert_eq!(col_uint(&frame, "atoms", "id"), [1, 2, 3, 4, 5, 6]);
        assert_eq!(col_uint(&frame, "atoms", "mol_id"), [1, 1, 2, 2, 3, 3]);
        assert_eq!(
            col_str(&frame, "atoms", "type"),
            ["A", "B", "A", "B", "A", "B"]
        );
    }

    #[test]
    fn coords_come_from_positions() {
        let target = Target::new(diatomic(), 2);
        let frame = assemble_frame(&[target], &positions(4)).unwrap();
        let xs: Vec<F> = frame
            .get("atoms")
            .unwrap()
            .get("x")
            .and_then(Column::as_float)
            .unwrap()
            .iter()
            .copied()
            .collect();
        assert_eq!(xs, [0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn heterogeneous_schema_is_union_filled() {
        // diatomic has a `type` column; argon does not — argon's atoms get a
        // null "" there.
        let targets = [Target::new(diatomic(), 2), Target::new(argon(), 2)];
        let frame = assemble_frame(&targets, &positions(6)).unwrap();

        assert_eq!(
            col_str(&frame, "atoms", "type"),
            ["A", "B", "A", "B", "", ""]
        );
        assert_eq!(
            frame.get("atoms").unwrap().validity("type"),
            Some(&[true, true, true, true, false, false][..])
        );
        assert_eq!(frame.get("atoms").unwrap().validity("element"), None);
        assert_eq!(col_uint(&frame, "atoms", "mol_id"), [1, 1, 2, 2, 3, 4]);
        // Bonds belong only to the two diatomics.
        assert_eq!(col_uint(&frame, "bonds", "atomi"), [0, 2]);
    }

    /// A later target's relation endpoints are shifted past every atom of the
    /// targets before it, on top of the per-copy offset.
    #[test]
    fn relations_of_a_later_target_are_offset_past_earlier_atoms() {
        let targets = [Target::new(argon(), 1), Target::new(diatomic(), 2)];
        let frame = assemble_frame(&targets, &positions(5)).unwrap();

        assert_eq!(col_uint(&frame, "bonds", "atomi"), [1, 3]);
        assert_eq!(col_uint(&frame, "bonds", "atomj"), [2, 4]);
        assert_eq!(col_uint(&frame, "bonds", "id"), [1, 2]);
        assert_eq!(col_uint(&frame, "atoms", "id"), [1, 2, 3, 4, 5]);
        assert_eq!(col_uint(&frame, "atoms", "mol_id"), [1, 2, 2, 3, 3]);
        let keys: Vec<&str> = frame.get("atoms").unwrap().keys().collect();
        assert_eq!(keys, ["id", "mol_id", "x", "y", "z", "element", "type"]);
    }

    #[test]
    fn no_template_falls_back_to_coords_only() {
        let target = Target::from_coords(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &[1.5, 1.5], 2);
        let frame = assemble_frame(&[target], &positions(4)).unwrap();

        assert_eq!(col_uint(&frame, "atoms", "id"), [1, 2, 3, 4]);
        assert!(frame.get("bonds").is_none());
    }

    fn col_i16(frame: &molrs::Frame, block: &str, key: &str) -> Vec<i16> {
        let column = frame.get(block).unwrap().get(key).unwrap();
        match column {
            Column::Int16(h) => h.array().iter().copied().collect(),
            other => panic!("expected i16, got {}", other.dtype()),
        }
    }

    /// A column in a narrow width tiles like any other; the hand-rolled
    /// tiler panicked on everything but the default widths.
    #[test]
    fn narrow_dtype_columns_are_tiled() {
        let mut frame = diatomic();
        let flag = ArrayD::from_shape_vec(vec![2], vec![-1i16, 1]).unwrap();
        frame
            .get_mut("atoms")
            .unwrap()
            .insert_column("flag", Column::from_i16(flag))
            .unwrap();
        let out = assemble_frame(&[Target::new(frame, 2)], &positions(4)).unwrap();
        assert_eq!(col_i16(&out, "atoms", "flag"), [-1, 1, -1, 1]);
    }

    /// Every canonical relation block is replayed, `pairs` included, with its
    /// endpoints offset per copy.
    #[test]
    fn pairs_are_replayed_with_offsets() {
        let mut frame = diatomic();
        let mut pairs = Block::new();
        pairs
            .insert(keys::ATOMI, Array1::from_vec(vec![0 as Idx]).into_dyn())
            .unwrap();
        pairs
            .insert(keys::ATOMJ, Array1::from_vec(vec![1 as Idx]).into_dyn())
            .unwrap();
        frame.insert("pairs", pairs);
        let out = assemble_frame(&[Target::new(frame, 2)], &positions(4)).unwrap();
        assert_eq!(col_uint(&out, "pairs", "atomi"), [0, 2]);
        assert_eq!(col_uint(&out, "pairs", "atomj"), [1, 3]);
    }

    /// A template whose atoms carry nothing beyond coordinates still gets its
    /// full row count when another template's column is filled in for it.
    #[test]
    fn a_coordinates_only_template_is_filled_to_its_row_count() {
        let mut bare = argon();
        bare.get_mut("atoms").unwrap().remove(keys::ELEMENT);
        let targets = [Target::new(bare, 3), Target::new(diatomic(), 1)];
        let out = assemble_frame(&targets, &positions(5)).unwrap();
        assert_eq!(col_str(&out, "atoms", "type"), ["", "", "", "A", "B"]);
    }

    /// Image flags describe the template's own coordinates, which packing
    /// replaces; they are not carried onto the packed copies.
    #[test]
    fn image_flags_are_not_carried() {
        let mut frame = argon();
        frame
            .get_mut("atoms")
            .unwrap()
            .insert(keys::IX, Array1::from_vec(vec![1i32]).into_dyn())
            .unwrap();
        let out = assemble_frame(&[Target::new(frame, 2)], &positions(2)).unwrap();
        assert!(!out.get("atoms").unwrap().contains_key(keys::IX));
    }

    /// Two templates carrying one column under different dtypes cannot share
    /// it; that is a named error, found before any packing.
    #[test]
    fn a_dtype_clash_is_named() {
        let mut labelled = diatomic();
        labelled
            .get_mut("atoms")
            .unwrap()
            .insert(
                "label",
                Array1::from_vec(vec!["a".to_string(), "b".to_string()]).into_dyn(),
            )
            .unwrap();
        let mut numbered = argon();
        numbered
            .get_mut("atoms")
            .unwrap()
            .insert("label", Array1::from_vec(vec![7i32]).into_dyn())
            .unwrap();
        let targets = [Target::new(labelled, 1), Target::new(numbered, 1)];
        assert!(matches!(
            check_templates(&targets),
            Err(PackError::TemplateColumns { .. })
        ));
    }
}
