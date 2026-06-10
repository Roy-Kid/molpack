"""Assemble a ``molrs.Frame`` from packed coordinates.

molpack's numeric core packs *coordinates only* — topology (bonds/angles/
dihedrals/impropers) and per-atom metadata are dropped at the Rust boundary
because the optimizer never needs them. The wheel rebuilds a topology-complete
frame here: each input template is replayed ``count`` times onto its packed
positions, atom-index columns are offset per copy, and ``id`` / ``mol_id`` are
regenerated.

The whole exchange goes through ``Frame.to_dict()`` / ``Frame.from_dict()`` —
the one serialization boundary every frame flavour (``molrs.Frame``,
``molpy.Frame``) already shares — so this module never touches block internals.
Frames carry no concatenation logic of their own; that is molpack's job and it
lives here. Force fields are out of scope: callers merge those separately.
"""

from __future__ import annotations

from typing import Any

import numpy as np

# Topology blocks, and which of their columns are 0/1-based atom indices that
# must be offset by the per-copy atom base. Every other column is carried
# verbatim.
_INDEX_COLUMNS: dict[str, tuple[str, ...]] = {
    "bonds": ("atomi", "atomj"),
    "angles": ("atomi", "atomj", "atomk"),
    "dihedrals": ("atomi", "atomj", "atomk", "atoml"),
    "impropers": ("atomi", "atomj", "atomk", "atoml"),
}
# Overwritten from packed positions, never carried from the template.
_COORD_COLUMNS = frozenset({"x", "y", "z"})
# Regenerated from scratch per system, never carried from the template.
_REGENERATED = frozenset({"id", "mol_id", "mol"})


def _import_molrs() -> Any:
    """Import the user's installed ``molrs`` (molpack does not link its wheel)."""
    import molrs

    return molrs


def _blocks(template: Any) -> dict[str, dict[str, Any]]:
    """Normalize a template to ``{block: {column: array}}``.

    Accepts the frame flavours :class:`Target` does: a ``molrs.Frame`` /
    ``molpy.Frame`` (via ``to_dict``) or a plain ``dict`` (used as-is, with or
    without a ``"blocks"`` envelope).
    """
    if hasattr(template, "to_dict"):
        return template.to_dict()["blocks"]
    if "blocks" in template and "atoms" not in template:
        return template["blocks"]
    return template


def _stamp_box(frame: Any, box: tuple[Any, Any] | None) -> Any:
    """Set an orthorhombic periodic box on ``frame`` from ``(min, max)`` corners."""
    if box is None:
        return frame
    molrs = _import_molrs()
    lo = np.asarray(box[0], dtype=float)
    hi = np.asarray(box[1], dtype=float)
    frame.box = molrs.Box.ortho(hi - lo, origin=lo, pbc=np.array([True, True, True]))
    return frame


def _frame_from_blocks(blocks: dict[str, dict[str, np.ndarray]]) -> Any:
    """Build a ``molrs.Frame`` from a ``{block: {column: array}}`` mapping."""
    molrs = _import_molrs()
    return molrs.Frame.from_dict({"blocks": blocks})


def _empty_like(sample: np.ndarray, length: int) -> np.ndarray:
    """A ``length``-long fill array matching ``sample``'s dtype.

    Used for a metadata column a mixture template lacks: numeric → ``0``,
    text → ``""``.
    """
    if sample.dtype.kind in "US":
        return np.full(length, "", dtype=sample.dtype)
    if sample.dtype.kind == "O":
        return np.full(length, "", dtype=object)
    return np.zeros(length, dtype=sample.dtype)


def coords_only_frame(
    positions: np.ndarray,
    elements: list[str],
    box: tuple[Any, Any] | None,
) -> Any:
    """Build a coordinates-only frame (single ``atoms`` block, no topology).

    The fallback when no source templates are retained (``.inp`` script packing).
    """
    atoms = {
        "id": np.arange(1, len(positions) + 1),
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "element": np.asarray(elements, dtype=str),
    }
    return _stamp_box(_frame_from_blocks({"atoms": atoms}), box)


def topology_frame(
    positions: np.ndarray,
    templates: list[tuple[Any, int]],
    box: tuple[Any, Any] | None,
) -> Any:
    """Replay each template ``count`` times onto ``positions``, with topology.

    ``templates`` is ``[(frame, count), ...]`` in pack order; ``positions`` is
    the ``(N, 3)`` packed array, atoms laid out template-by-template, copy-by-copy.
    """
    # Normalize each template to its block mapping once: [(blocks, count), ...].
    plans = [(_blocks(frame), count) for frame, count in templates]

    expected = sum(len(blocks["atoms"]["x"]) * count for blocks, count in plans)
    if expected != len(positions):
        raise ValueError(
            f"packed positions hold {len(positions)} atoms but templates "
            f"account for {expected}; template/count order must match pack()"
        )

    # Union of carried atom columns across all templates (first-seen order),
    # each keyed to a sample array for its dtype. A template missing one of
    # these (a mixture of differing schemas) contributes a default-filled run.
    carried: dict[str, np.ndarray] = {}
    for blocks, _ in plans:
        for name, column in blocks["atoms"].items():
            if name not in _COORD_COLUMNS and name not in _REGENERATED:
                carried.setdefault(name, np.asarray(column))

    # block -> column -> list of array parts (one per template), concatenated last.
    parts: dict[str, dict[str, list[np.ndarray]]] = {}

    def add(block: str, column: str, values: np.ndarray) -> None:
        parts.setdefault(block, {}).setdefault(column, []).append(values)

    atom_base = 0  # atoms emitted so far (also the per-copy index offset)
    mol_base = 0  # molecules emitted so far
    cursor = 0  # next row in `positions`

    for blocks, count in plans:
        atoms = blocks["atoms"]
        n = len(atoms["x"])
        span = count * n
        coords = positions[cursor : cursor + span]
        cursor += span

        # Atoms: tile carried columns (default-filling any this template lacks),
        # fill coords/id/mol_id for all copies.
        for name, sample in carried.items():
            column = atoms[name] if name in atoms else _empty_like(sample, n)
            add("atoms", name, np.tile(column, count))
        add("atoms", "x", coords[:, 0])
        add("atoms", "y", coords[:, 1])
        add("atoms", "z", coords[:, 2])
        add("atoms", "id", np.arange(atom_base + 1, atom_base + span + 1))
        add(
            "atoms",
            "mol_id",
            np.repeat(np.arange(mol_base + 1, mol_base + count + 1), n),
        )

        # Topology: tile columns; shift index columns by each copy's atom base.
        for block, index_columns in _INDEX_COLUMNS.items():
            table = blocks.get(block)
            if not table:
                continue
            rows = len(next(iter(table.values())))
            copy_offset = np.repeat(np.arange(count) * n, rows) + atom_base
            for name, column in table.items():
                if name == "id":
                    continue
                tiled = np.tile(column, count)
                add(
                    block, name, tiled + copy_offset if name in index_columns else tiled
                )

        atom_base += span
        mol_base += count

    blocks_out: dict[str, dict[str, np.ndarray]] = {}
    for block, columns in parts.items():
        merged = {name: np.concatenate(chunks) for name, chunks in columns.items()}
        if block in _INDEX_COLUMNS:  # regenerate a contiguous topology id
            merged["id"] = np.arange(1, len(next(iter(merged.values()))) + 1)
        blocks_out[block] = merged

    return _stamp_box(_frame_from_blocks(blocks_out), box)
