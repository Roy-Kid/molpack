# Restraints and PBC

Restraints are soft penalties that guide atoms into allowed regions. They can
be attached to a target, to a subset of atoms on each target copy, or globally
on the engine entry.

## Whole-target restraints

```rust
use molpack::{InsideBoxRestraint, Target};

let target = Target::from_coords(positions, radii, 100)
    .with_restraint(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [40.0, 40.0, 40.0],
        [false, false, false],
    ));
```

Built-in geometric restraints include boxes, spheres, outside-sphere
exclusions, and plane half-spaces. Collective restraints match an entire
species to a distribution profile; see [Concepts](../concepts.md).

## Atom-subset restraints

Atom-subset restraints apply to selected atoms of every copy. Indices are
0-based:

```rust
use molpack::{BelowPlaneRestraint, Target};

let target = Target::from_coords(positions, radii, 100)
    .with_atom_restraint(
        &[0, 1],
        BelowPlaneRestraint::new([0.0, 0.0, 1.0], 2.0),
    );
```

If you are translating from a Packmol `.inp` `atoms ... end atoms` block,
subtract 1 from each atom index.

## Global restraints

Attach a restraint to every target through the engine entry:

```rust
use molpack::{GenCanPack, InsideSphereRestraint, PackEngine, Target};

let result = GenCanPack::new()
    .with_global_restraint(InsideSphereRestraint::new([20.0, 20.0, 20.0], 30.0))
    .run(&[a, b], 200)?;
```

This is equivalent to cloning the same restraint onto every target before
packing.

## Periodic boxes

There are two ways to declare periodic boundary conditions.

For fully periodic boxes, set PBC on the engine entry:

```rust
use molpack::{GenCanPack, PackEngine};

let engine = GenCanPack::new().with_periodic_box([0.0; 3], [30.0; 3], [true; 3]);
```

For per-axis control, set periodic flags on `InsideBoxRestraint`:

```rust
use molpack::{InsideBoxRestraint, Target};

let slab = Target::from_coords(positions, radii, 100)
    .with_restraint(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [30.0, 30.0, 100.0],
        [true, true, false],
    ));
```

Multiple PBC declarations must agree exactly. Invalid or conflicting boxes
return typed `PackError` variants.
