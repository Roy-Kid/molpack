"""StlRegion Python binding — attach path, not a packing run."""

from __future__ import annotations

from pathlib import Path

import molrs
import numpy as np
import pytest

import molpack


def _ascii_unit_cube() -> str:
    # Watertight [0,1]³ Å, 12 triangles (same winding as Rust cube_tris).
    faces = [
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),
        ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)),
        ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)),
        ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)),
        ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)),
        ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),
    ]
    lines = ["solid cube"]
    for a, b, c, d in faces:
        lines += _facet(a, b, c)
        lines += _facet(a, c, d)
    lines.append("endsolid cube")
    return "\n".join(lines) + "\n"


def _facet(a, b, c) -> list[str]:
    return [
        "  facet normal 0 0 0",
        "    outer loop",
        f"      vertex {a[0]} {a[1]} {a[2]}",
        f"      vertex {b[0]} {b[1]} {b[2]}",
        f"      vertex {c[0]} {c[1]} {c[2]}",
        "    endloop",
        "  endfacet",
    ]


def _one_atom_target():
    frame = molrs.Frame(
        {
            "atoms": {
                "x": np.array([0.5]),
                "y": np.array([0.5]),
                "z": np.array([0.5]),
            }
        }
    )
    return molpack.Target(frame, 1)


class TestStlRegion:
    def test_from_file_unit_cube(self, tmp_path: Path):
        p = tmp_path / "cube.stl"
        p.write_text(_ascii_unit_cube())
        r = molpack.StlRegion.from_file(p)
        assert "StlRegion" in repr(r)
        r2 = molpack.StlRegion.from_file(p, scale=1.0)
        assert "StlRegion" in repr(r2)

    def test_missing_path_oserror(self, tmp_path: Path):
        with pytest.raises(OSError):
            molpack.StlRegion.from_file(tmp_path / "nope.stl")

    def test_leaky_mesh_valueerror(self, tmp_path: Path):
        p = tmp_path / "tri.stl"
        p.write_text(
            "solid t\n  facet normal 0 0 0\n    outer loop\n"
            "      vertex 0 0 0\n      vertex 1 0 0\n      vertex 0 1 0\n"
            "    endloop\n  endfacet\nendsolid t\n"
        )
        with pytest.raises(ValueError):
            molpack.StlRegion.from_file(p)

    def test_no_empty_constructor(self):
        with pytest.raises(TypeError):
            molpack.StlRegion()

    def test_bad_scale(self, tmp_path: Path):
        p = tmp_path / "cube.stl"
        p.write_text(_ascii_unit_cube())
        with pytest.raises(ValueError):
            molpack.StlRegion.from_file(p, scale=0.0)

    def test_with_restraint_and_atom_restraint(self, tmp_path: Path):
        p = tmp_path / "cube.stl"
        p.write_text(_ascii_unit_cube())
        stl = molpack.StlRegion.from_file(p)
        assert not callable(getattr(stl, "f", None))
        target = _one_atom_target()
        target.with_restraint(stl)
        target.with_atom_restraint([0], stl)
