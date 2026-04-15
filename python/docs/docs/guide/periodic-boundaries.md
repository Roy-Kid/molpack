# Periodic boundaries

By default the packer works under **free boundary conditions** — atoms
are not wrapped and the only geometric limits come from the restraints
you attach. Use periodic boundaries (PBC) when packing for MD input.

## Enabling PBC

Two equivalent forms:

```python
# 1. explicit min / max corners
packer = Packer().with_pbc(
    min=[0.0, 0.0, 0.0],
    max=[30.0, 30.0, 30.0],
)

# 2. origin at zero, cell lengths only
packer = Packer().with_pbc_box([30.0, 30.0, 30.0])
```

Only orthorhombic cells are supported. Triclinic / non-orthogonal
boxes are a road-map item; for now rotate / skew coordinates in
post-processing.

## Semantics

Under PBC, the pairwise distance evaluator applies minimum-image
wrapping, so atoms near opposite faces of the cell "see" each other
through the periodic images. The `tolerance` setting still applies
and is checked against the wrapped distance.

Restraints (`InsideBox`, `InsideSphere`, …) are evaluated in the
**unwrapped** frame — they describe the geometric region as defined,
regardless of the periodic cell.

## Errors

A zero-length axis, or a `max < min` pair, raises:

```python
ValueError: InvalidPBCBox(...)
```

Validation is performed in `pack()`, not in the builder, so the same
`Packer` instance can be reused with different PBC settings.

## Choosing a box

A common pattern: pack into a tight `InsideBox` restraint matching
the desired cell, then enable PBC on the packer itself so atoms near
the faces respect each other:

```python
cell = [30.0, 30.0, 30.0]

target = target.with_constraint(InsideBox([0, 0, 0], cell))
packer = Packer(tolerance=2.0).with_pbc_box(cell)
result = packer.pack([target], max_loops=200, seed=42)
```

The restraint keeps atoms inside the cell; the PBC setting ensures
distances are evaluated across the periodic faces.
