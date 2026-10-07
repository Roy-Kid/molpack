# Periodic boundaries

By default the packer works under **free boundary conditions** — atoms
are not wrapped and the only geometric limits come from the restraints
you attach. Use periodic boundaries (PBC) when packing for MD input.

## Enabling PBC

The periodic box is declared on the engine entry (Packmol's `pbc`
keyword). `with_periodic_box` is a shared builder, so it reads the same on
`GenCanPack`, `CbmcGrow` and `LatticeGrow`:

```python
from molpack import GenCanPack

packer = GenCanPack().with_periodic_box([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
```

Only orthorhombic cells are supported through this builder; a triclinic
cell is `with_cell(lengths, angles, pbc)`.

A region does not declare periodicity — it only confines. To keep every
atom centre inside the cell as well, attach a `molrs.core.Cuboid` with the same
bounds (or broadcast it with `with_global_restraint`):

```python
import molrs

cell = molrs.core.Cuboid([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])   # origin, lengths
target = target.with_restraint(cell)
```

## Semantics

Under PBC, the pairwise distance evaluator applies minimum-image
wrapping on the periodic axes, so atoms near opposite faces of the
cell "see" each other through the periodic images. The `tolerance`
setting still applies and is checked against the wrapped distance.

Regions are evaluated in the **unwrapped** frame — they describe the
solid as defined, regardless of the periodic cell. The one exception is
`molrs.core.SphereUnion` built with a `box`: its spheres are minimum-image on
the box's periodic axes, so a void computed from beads in a periodic frame
is periodic too.

## Errors

A zero-length axis on a periodic box, or `max < min` on any axis,
raises `InvalidPBCBoxError` at `run()` time:

```python
from molpack import InvalidPBCBoxError

try:
    packer.run(targets, max_loops=200)
except InvalidPBCBoxError as e:
    ...
```

The periodic box is declared in exactly one place — on the entry — so there
is no second declaration to conflict with it. Typed errors inherit from
`molpack.PackError` (which itself is a `RuntimeError` subclass), so a blanket
`except PackError` catches any packing failure.

## Choosing a box

A common pattern: one periodic cell on the entry, and the same cuboid on
every target that must stay inside it:

```python
cell_min = [0.0, 0.0, 0.0]
cell_len = [30.0, 30.0, 30.0]

box = molrs.core.Cuboid(cell_min, cell_len)
target = target.with_restraint(box)

result = (
    GenCanPack()
    .with_seed(42)
    .with_periodic_box(cell_min, [a + b for a, b in zip(cell_min, cell_len)])
    .run([target], max_loops=200)
)
```

Or broadcast the cuboid globally via `GenCanPack.with_global_restraint(box)`
when several species share the same cell.
