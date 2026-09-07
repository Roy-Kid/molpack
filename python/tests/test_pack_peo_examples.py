"""Smoke the three PEO scene examples: linear, mixed topology, STL cavity."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

import molpack

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
sys.path.insert(0, str(EXAMPLES))

import pack_peo_linear as linear  # noqa: E402
import pack_peo_mix as mix  # noqa: E402
import pack_peo_stl as stl  # noqa: E402


class TestLinear:
    def test_lattice_then_push_one_dimer(self):
        polymer = linear.make_linear(2, seed=42)
        state = linear.pack_linear(2, 1, 0.2, 42)
        assert state.natoms == polymer.n_atoms


class TestMix:
    def test_rejects_empty_species(self):
        with pytest.raises(ValueError, match="at least one"):
            mix.pack_mix(2, 2, 0, 1, 0.2, 1)
        with pytest.raises(ValueError, match="at least one"):
            mix.pack_mix(2, 2, 1, 0, 0.2, 1)

    def test_two_targets_one_box(self):
        lin = mix.make_linear(2, seed=42)
        star = mix.make_star(2, seed=42)
        state = mix.pack_mix(2, 2, 1, 1, 0.2, 42)
        assert state.natoms == lin.n_atoms + star.n_atoms


class TestStl:
    def test_grows_in_the_shipped_dendrite(self):
        polymer = stl.make_linear(2, seed=42)
        state = stl.pack_stl(2, 1, stl.MESH_EDGE, 42)
        assert state.natoms == polymer.n_atoms


class TestStlRegionAttach:
    def test_shipped_mesh_is_loadable(self):
        region = molpack.StlRegion.from_file(stl.MESH)
        assert "StlRegion" in repr(region)

    def test_lattice_grows_with_stl_region(self):
        cavity = molpack.StlRegion.from_file(stl.MESH)
        polymer = stl.make_linear(2, seed=1)
        target = molpack.Target(polymer.to_frame(), 1).with_restraint(cavity)
        prior = molpack.TorsionPrior.three_state_from_c_inf(stl.PEO_C_INF, stl.TET)
        grown = (
            molpack.LatticeGrow(prior)
            .with_seed(1)
            .with_tolerance(2.0)
            .with_periodic_box([0.0] * 3, [stl.MESH_EDGE] * 3)
            .with_progress(False)
            .run([target], max_loops=40)
        )
        assert grown.natoms == polymer.n_atoms
