# Periodic boundaries

By default the packer works under **free boundary conditions** — atoms
are not wrapped and the only geometric limits come from the restraints
you attach. Use periodic boundaries (PBC) when packing for MD input.

## Enabling PBC

PBC is declared on the `InsideBoxRestraint` via the `periodic` keyword:

```python
from molpack import InsideBoxRestraint

# Fully periodic orthorhombic cell
box = InsideBoxRestraint(
    [0.0, 0.0, 0.0],
    [30.0, 30.0, 30.0],
    periodic=(True, True, True),
)
```

`periodic` is a 3-tuple of booleans — one per axis. Only orthorhombic
cells are supported.

Per-axis PBC is possible — e.g. slab geometry with in-plane PBC and
open Z:

```python
slab = InsideBoxRestraint(
    [0.0, 0.0, 0.0],
    [30.0, 30.0, 100.0],
    periodic=(True, True, False),
)
```

## Semantics

Under PBC, the pairwise distance evaluator applies minimum-image
wrapping on the periodic axes, so atoms near opposite faces of the
cell "see" each other through the periodic images. The `tolerance`
setting still applies and is checked against the wrapped distance.

Restraints themselves (`InsideBoxRestraint`, `InsideSphereRestraint`,
…) are evaluated in the **unwrapped** frame — they describe the
geometric region as defined, regardless of the periodic cell.

## System-wide PBC derivation

At `pack()` time the packer scans every restraint on every target for
a declared periodic box. The rules are:

- **Zero declarations** — non-periodic run.
- **One declaration** — its bounds define the system PBC.
- **Multiple declarations** — they must all agree (same bounds, same
  per-axis flags). Any mismatch raises `ConflictingPeriodicBoxesError`.

## Errors

A zero-length axis on a periodic box, or `max < min` on any axis,
raises `InvalidPBCBoxError` at `pack()` time:

```python
from molpack import InvalidPBCBoxError

try:
    packer.pack(...)
except InvalidPBCBoxError as e:
    ...
```

Both typed errors inherit from `molpack.PackError` (which itself is
a `RuntimeError` subclass) — a blanket `except PackError` catches any
packing failure.

## Choosing a box

A common pattern: pack into a single periodic `InsideBoxRestraint`
matching the desired cell. The restraint confines atoms softly and
simultaneously declares the PBC:

```python
cell_min = [0.0, 0.0, 0.0]
cell_max = [30.0, 30.0, 30.0]

box = InsideBoxRestraint(cell_min, cell_max, periodic=(True, True, True))
target = target.with_restraint(box)

result = Molpack().with_seed(42).pack([target], max_loops=200)
```

Or broadcast it globally via `Molpack.with_global_restraint(box)` when
several species share the same cell.
