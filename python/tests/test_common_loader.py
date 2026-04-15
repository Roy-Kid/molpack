"""Tests for the examples/_common.py PDB loader helper."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent / "examples"
DATA_ROOT = HERE.parent.parent / "examples"

sys.path.insert(0, str(EXAMPLES))

from _common import (  # noqa: E402
    ION_ALIASES,
    VDW_RADII,
    _element_from_atom_name,
    read_pdb_as_arrays,
    vdw_radius,
)


class TestVdwRadius:
    def test_known_elements(self):
        assert vdw_radius("O") == pytest.approx(1.52)
        assert vdw_radius("H") == pytest.approx(1.20)
        assert vdw_radius("C") == pytest.approx(1.70)
        assert vdw_radius("N") == pytest.approx(1.55)

    def test_case_insensitive(self):
        assert vdw_radius("o") == vdw_radius("O")
        assert vdw_radius("na") == vdw_radius("NA")

    def test_ion_alias_table_covers_expected_keys(self):
        assert "SOD" in ION_ALIASES
        assert ION_ALIASES["SOD"] == "NA"
        assert ION_ALIASES["CLA"] == "CL"

    def test_unknown_element_falls_back_to_x(self):
        assert vdw_radius("Zz") == VDW_RADII["X"]


class TestElementFromAtomName:
    def test_single_letter(self):
        assert _element_from_atom_name(" O  ") == "O"
        assert _element_from_atom_name(" N  ") == "N"

    def test_protein_alpha_carbon(self):
        # "CA" in a protein context is an alpha carbon, not calcium.
        assert _element_from_atom_name(" CA ") == "C"

    def test_strips_digits(self):
        assert _element_from_atom_name(" OG1") == "O"
        assert _element_from_atom_name("1HB1") == "H"

    def test_ion_alias_sodium(self):
        assert _element_from_atom_name("SOD ") == "NA"

    def test_ion_alias_chloride(self):
        assert _element_from_atom_name("CLA ") == "CL"

    def test_empty_returns_fallback(self):
        assert _element_from_atom_name("    ") == "X"
        assert _element_from_atom_name("1234") == "X"


class TestReadPdb:
    def test_water_three_atoms(self):
        pos, rad, els = read_pdb_as_arrays(DATA_ROOT / "pack_mixture" / "water.pdb")
        assert pos.shape == (3, 3)
        assert pos.dtype == np.float64
        assert rad.shape == (3,)
        assert rad.dtype == np.float64
        assert els == ["H", "H", "O"]

    def test_radii_match_elements(self):
        _, rad, els = read_pdb_as_arrays(DATA_ROOT / "pack_mixture" / "water.pdb")
        expected = np.array([vdw_radius(e) for e in els])
        np.testing.assert_allclose(rad, expected)

    def test_sodium_ion(self):
        pos, rad, els = read_pdb_as_arrays(
            DATA_ROOT / "pack_solvprotein" / "sodium.pdb"
        )
        assert els == ["NA"]
        assert rad[0] == pytest.approx(vdw_radius("NA"))

    def test_chloride_ion(self):
        _, rad, els = read_pdb_as_arrays(
            DATA_ROOT / "pack_solvprotein" / "chloride.pdb"
        )
        assert els == ["CL"]
        assert rad[0] == pytest.approx(vdw_radius("CL"))

    def test_protein_contains_only_standard_elements(self):
        _, _, els = read_pdb_as_arrays(DATA_ROOT / "pack_solvprotein" / "protein.pdb")
        assert set(els) <= {"C", "H", "N", "O", "S"}

    def test_reads_coordinates_correctly(self):
        pos, _, _ = read_pdb_as_arrays(DATA_ROOT / "pack_mixture" / "water.pdb")
        # water.pdb:
        # HETATM    1  H   HOH ...      9.626   6.787  12.673
        # HETATM    2  H   HOH ...      9.626   8.420  12.673
        # HETATM    3  O   HOH ...     10.203   7.604  12.673
        np.testing.assert_allclose(pos[0], [9.626, 6.787, 12.673])
        np.testing.assert_allclose(pos[2], [10.203, 7.604, 12.673])

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            read_pdb_as_arrays("/nonexistent/path.pdb")

    def test_empty_pdb_raises(self, tmp_path):
        p = tmp_path / "empty.pdb"
        p.write_text("HEADER\nEND\n")
        with pytest.raises(ValueError, match="no ATOM/HETATM"):
            read_pdb_as_arrays(p)
