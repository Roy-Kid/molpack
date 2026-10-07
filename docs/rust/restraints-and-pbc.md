# Restraints and PBC

Restraints are soft penalties that guide atoms into allowed regions. They can
be attached to a target, to a subset of atoms on each target copy, or globally
on the engine entry. Geometry is a molrs [`Region`](https://docs.rs/molcrafts-molrs)
— `Sphere`, `Cuboid`, `Parallelepiped`, `HalfSpace`, `Cylinder`, `Ellipsoid`,
`Polyhedron`, `SphereUnion`, and their `AndRegion` / `OrRegion` / `NotRegion`
compositions — and molpack's one geometric restraint, `RegionRestraint`, says
"stay inside it".

## Whole-target restraints

```rust
use std::sync::Arc;
use molpack::{RegionRestraint, Target};
use molrs::core::Cuboid;
use ndarray::array;

let cube = Cuboid::new(array![0.0, 0.0, 0.0], array![40.0, 40.0, 40.0]); // origin, lengths
let target = Target::from_coords(positions, radii, 100)
    .with_restraint(RegionRestraint(Arc::new(cube)));
```

Outside a shape is its complement: `NotRegion::new(Arc::new(sphere))`.
Collective restraints match an entire species to a distribution profile; see
[Concepts](../concepts.md).

## Atom-subset restraints

Atom-subset restraints apply to selected atoms of every copy. Indices are
0-based:

```rust
use std::sync::Arc;
use molpack::{RegionRestraint, Target};
use molrs::core::HalfSpace;

// z <= 2: the half-space behind the plane through (0, 0, 2) with normal +z.
let below = HalfSpace::new([0.0, 0.0, 1.0], [0.0, 0.0, 2.0])?;
let target = Target::from_coords(positions, radii, 100)
    .with_atom_restraint(&[0, 1], RegionRestraint(Arc::new(below)));
```

If you are translating from a Packmol `.inp` `atoms ... end atoms` block,
subtract 1 from each atom index.

## Global restraints

Attach a restraint to every target through the engine entry:

```rust
use std::sync::Arc;
use molpack::{GenCanPack, PackEngine, RegionRestraint, Target};
use molrs::core::Sphere;
use ndarray::array;

let ball = Sphere::new(array![20.0, 20.0, 20.0], 30.0);
let result = GenCanPack::new()
    .with_global_restraint(RegionRestraint(Arc::new(ball)))
    .run(&[a, b], 200)?;
```

This is equivalent to cloning the same restraint onto every target before
packing.

## Periodic boxes

Periodic boundary conditions are declared on the engine entry; a region only
confines:

```rust
use molpack::{GenCanPack, PackEngine};

let engine = GenCanPack::new().with_periodic_box([0.0; 3], [30.0; 3], [true; 3]);
```

Per-axis periodicity is the third argument (`[true, true, false]` for a slab
open along z). A restraint used with periodic boundaries has to mean the same
thing in every image, and molpack refuses by name — `RestraintAcrossPeriodicAxis`
— one that does not. Two shapes do: a restraint that **confines** atoms to a
single image (a box, a sphere, a cell), and one that **repeats** along the
lattice vector (a `HalfSpace` whose plane runs parallel to it — the slab case).
A half-space across a periodic axis is refused whichever way it was built, the
`.inp` grammar's plane or a lifted molrs region. Putting a confining region in
the intersection rescues an open one: `Cuboid & HalfSpace` is bounded, so it
is usable along every axis. A triclinic cell is declared once, as a wall and as
the lattice, by `CellRestraint`:

```rust
use molpack::{CellRestraint, Target};

let cell = CellRestraint::from_lengths_angles([26.0; 3], [90.0, 90.0, 120.0], [true; 3])?;
let target = Target::from_coords(positions, radii, 100).with_restraint(cell);
```

Invalid or conflicting boxes return typed `PackError` variants.
