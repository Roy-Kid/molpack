"""molrs regions as molpack restraints — the attach path and the lift.

molpack has no region class of its own: a region is a molrs object
(``Sphere``, ``Cuboid``, ``Parallelepiped``, ``HalfSpace``, ``Cylinder``,
``Ellipsoid``, ``Polyhedron``, ``SphereUnion`` or a ``&`` / ``|`` / ``~``
composition) that crosses the wheel boundary as a ``molrs.RegionRef`` capsule
and is lifted to "stay inside" by ``RegionRestraint``.
"""

from __future__ import annotations

import molrs
import numpy as np
import pytest

import molpack


def _one_atom_frame(x=0.5, y=0.5, z=0.5):
    return molrs.store.Frame(
        {
            "atoms": {
                "x": np.array([x]),
                "y": np.array([y]),
                "z": np.array([z]),
                "element": ["X"],
            }
        }
    )


class TestRegionAttach:
    def test_every_region_class_attaches(self):
        z = [0.0, 0.0, 0.0]
        regions = [
            molrs.spatial.Sphere(z, 5.0),
            molrs.spatial.Cuboid(z, [10.0, 10.0, 10.0]),
            molrs.spatial.Parallelepiped.cube(10.0, z),
            molrs.spatial.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 8.0]),
            molrs.spatial.Cylinder(z, [0.0, 0.0, 1.0], 4.0, 10.0),
            molrs.spatial.Ellipsoid(z, [5.0, 6.0, 7.0]),
            molrs.spatial.SphereUnion(
                np.array([[1.0, 1.0, 1.0], [3.0, 1.0, 1.0]]), 2.0
            ),
            ~molrs.spatial.Sphere(z, 1.0) & molrs.spatial.Cuboid(z, [10.0, 10.0, 10.0]),
        ]
        target = molpack.Target(_one_atom_frame(), 1)
        for region in regions:
            assert target.with_restraint(region) is not None, type(region)
            assert target.with_atom_restraint([0], region) is not None, type(region)
            assert molpack.GenCanPack().with_global_restraint(region) is not None

    def test_region_is_lifted_not_duck_typed(self):
        sphere = molrs.spatial.Sphere([0.0, 0.0, 0.0], 5.0)
        assert not callable(getattr(sphere, "f", None))
        assert callable(sphere._ffi_regionref_capsule)

    def test_no_molpack_geometry_classes(self):
        for name in (
            "StlRegion",
            "InsideBoxRestraint",
            "InsideSphereRestraint",
            "OutsideSphereRestraint",
            "AbovePlaneRestraint",
            "BelowPlaneRestraint",
        ):
            assert not hasattr(molpack, name), name

    def test_non_region_without_f_fg_is_a_typeerror(self):
        with pytest.raises(TypeError, match="expected a restraint"):
            molpack.Target(_one_atom_frame(), 1).with_atom_restraint([0], object())


class TestRegionPacking:
    def test_confines_inside_sphere(self):
        centre = [10.0, 10.0, 10.0]
        radius = 4.0
        ball = molrs.spatial.Sphere(centre, radius)
        target = molpack.Target(_one_atom_frame(), 8).with_restraint(ball)
        state = (
            molpack.GenCanPack()
            .with_seed(3)
            .with_tolerance(2.0)
            .with_precision(1e-4)
            .with_progress(False)
            .run([target], max_loops=200)
        )
        # Every atom centre stays inside the sphere (a soft wall: allow the
        # sub-tolerance excursion the precision permits).
        assert (
            molrs.spatial.Sphere(centre, radius + 0.5).contains(state.positions).all()
        )

    def test_void_of_a_sphere_union_is_respected(self):
        beads = np.array([[10.0, 10.0, 10.0]])
        polymer = molrs.spatial.SphereUnion(beads, 4.0)
        void = ~polymer & molrs.spatial.Cuboid([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        target = molpack.Target(_one_atom_frame(), 6).with_restraint(void)
        state = (
            molpack.GenCanPack()
            .with_seed(5)
            .with_tolerance(2.0)
            .with_precision(1e-4)
            .with_progress(False)
            .run([target], max_loops=200)
        )
        d = np.linalg.norm(state.positions - beads[0], axis=1)
        assert (d > 4.0 - 0.5).all(), d
