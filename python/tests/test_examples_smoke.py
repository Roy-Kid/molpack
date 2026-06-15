"""Smoke tests: every example script compiles and its PDB assets parse.

Full workloads (the actual pack runs) are too large for the unit-test
suite — run ``python examples/pack_<name>.py`` to exercise them.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent / "examples"
DATA_ROOT = HERE.parent.parent / "examples"

sys.path.insert(0, str(EXAMPLES))
from _common import read_pdb_as_arrays  # noqa: E402

EXAMPLE_SCRIPTS = [
    "pack_water_cube.py",
    "pack_mixture.py",
    "pack_bilayer.py",
    "pack_interface.py",
    "pack_spherical.py",
    "pack_solvprotein.py",
]

# Maps every example script to the PDB fixtures it depends on.
EXAMPLE_FIXTURES = {
    "pack_mixture.py": ["pack_mixture/water.pdb", "pack_mixture/urea.pdb"],
    "pack_bilayer.py": ["pack_bilayer/water.pdb", "pack_bilayer/palmitoil.pdb"],
    "pack_interface.py": [
        "pack_interface/water.pdb",
        "pack_interface/chloroform.pdb",
        "pack_interface/t3.pdb",
    ],
    "pack_spherical.py": ["pack_spherical/water.pdb", "pack_spherical/palmitoil.pdb"],
    "pack_solvprotein.py": [
        "pack_solvprotein/protein.pdb",
        "pack_solvprotein/water.pdb",
        "pack_solvprotein/sodium.pdb",
        "pack_solvprotein/chloride.pdb",
    ],
}


@pytest.mark.parametrize("script", EXAMPLE_SCRIPTS)
def test_example_script_parses(script):
    source = (EXAMPLES / script).read_text(encoding="utf-8")
    ast.parse(source)


@pytest.mark.parametrize("script", EXAMPLE_SCRIPTS)
def test_example_script_has_main(script):
    source = (EXAMPLES / script).read_text(encoding="utf-8")
    tree = ast.parse(source)
    fns = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert "main" in fns


@pytest.mark.parametrize(
    "script,fixtures",
    sorted(EXAMPLE_FIXTURES.items()),
)
def test_example_fixtures_exist(script, fixtures):
    for rel in fixtures:
        p = DATA_ROOT / rel
        assert p.exists(), f"{script}: missing fixture {p}"


@pytest.mark.parametrize(
    "script,fixtures",
    sorted(EXAMPLE_FIXTURES.items()),
)
def test_example_fixtures_parse(script, fixtures):
    for rel in fixtures:
        pos, rad, els = read_pdb_as_arrays(DATA_ROOT / rel)
        assert pos.shape[0] > 0
        assert pos.shape[1] == 3
        assert rad.shape[0] == pos.shape[0]
        assert len(els) == pos.shape[0]
