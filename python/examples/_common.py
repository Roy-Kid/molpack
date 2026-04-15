"""Shared helpers for the molpack Python examples.

Keeps the examples self-contained: a minimal PDB loader and a small
van der Waals radius table let every example run without depending on
`molcrafts-molrs` (or any other I/O library). If you already have
`molrs` installed, replace ``read_pdb_as_arrays`` with
``molrs.read_pdb`` plus ``block.view("x")`` etc.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# Bondi (1964) van der Waals radii in Å, rounded to 2 dp. Covers the
# elements that appear in the bundled example PDBs.
VDW_RADII: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "NA": 2.27,
    "K": 2.75,
    "X": 1.50,  # fallback
}


def vdw_radius(element: str) -> float:
    """Return the Bondi VdW radius (Å) for an element symbol."""
    return VDW_RADII.get(element.upper(), VDW_RADII["X"])


# CHARMM-style ion aliases that cannot be recovered from the first letter
# of the atom name alone.
ION_ALIASES: dict[str, str] = {
    "SOD": "NA",  # sodium
    "CLA": "CL",  # chloride
    "POT": "K",  # potassium
    "CAL": "CA",  # calcium  (true element)
}


def _element_from_atom_name(name: str) -> str:
    """Best-effort element symbol from a PDB atom-name (columns 13–16).

    The PDB standard places element symbols in columns 77–78; when those
    are blank we fall back to the atom name. Convention for protein
    records: strip digits, then take the first letter — so ``"CA"``,
    ``"CB"``, ``"CG"`` all map to ``C`` (alpha / beta / gamma carbon);
    ``"OG1"`` maps to ``O``; ``"NE2"`` to ``N``. This heuristic mislabels
    two-letter elements like sodium (``"NA"``) and calcium (``"CA"``)
    when they appear under ``ATOM``/``HETATM`` without an element field
    — such PDBs should fill in columns 77–78 explicitly.
    """
    stripped = "".join(c for c in name.strip() if c.isalpha()).upper()
    if stripped in ION_ALIASES:
        return ION_ALIASES[stripped]
    return stripped[:1] if stripped else "X"


def read_pdb_as_arrays(path: str | Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Load a PDB file into (positions, radii, elements).

    Parses ``ATOM`` and ``HETATM`` records. Elements come from columns
    77–78 when present; otherwise from the stripped atom-name
    (columns 13–16). Radii are looked up from the VdW table above.

    Parameters
    ----------
    path
        Path to a PDB file.

    Returns
    -------
    positions : ndarray, shape (N, 3), dtype float64
    radii     : ndarray, shape (N,), dtype float64
    elements  : list[str]
    """
    positions: list[list[float]] = []
    elements: list[str] = []

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                element = _element_from_atom_name(line[12:16])
            positions.append([x, y, z])
            elements.append(element)

    if not positions:
        raise ValueError(f"no ATOM/HETATM records in {path}")

    pos = np.asarray(positions, dtype=np.float64)
    rad = np.asarray([vdw_radius(e) for e in elements], dtype=np.float64)
    return pos, rad, elements
