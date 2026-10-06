# Quickstart

Pack **100 water molecules** into a **40 Å** cube. This walkthrough uses the
Python package — the shortest path from a loaded frame to a packed result. The
same model is available as a [CLI script](cli/) or the [Rust builder](rust/).

## 1. Install

```bash
pip install molcrafts-molpack
```

This installs the packing engine and pulls in `molcrafts-molrs` for the frame
type plus PDB/XYZ I/O.

## 2. Load or build a template

A **template** is one copy of the molecule you want many of — its atom
positions and element symbols, in any orientation. molpack takes the template
as a `molrs.store.Frame`, the shared MolCrafts container for atomic data, so any
loader molrs supports will do. Read one from a Protein Data Bank (PDB) file:

```python
import molrs

frame = molrs.io.read_pdb("water.pdb")
```

Or build one in memory, with no file involved — here a water molecule with the
oxygen at the origin:

```python
import molrs
import numpy as np

frame = molrs.store.Frame({
    "atoms": {
        "x": np.array([0.00, 0.96, -0.24]),
        "y": np.array([0.00, 0.00, 0.93]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
})
```

Coordinates are in ångström (Å) throughout molpack.

## 3. Define the target

A `Target` is one molecule species plus the number of copies to place. Every
mobile target needs a spatial restraint:

```python
import molrs
from molpack import Target

box = molrs.spatial.Cuboid([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])  # a molrs region
water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(box)
)
```

!!! warning "Missing restraints"
    Without a spatial restraint (or a global PBC box), initial placement has to
    invent a huge free-space region and the run can become impractical.

## 4. Pack

```python
from molpack import GenCanPack

result = (
    GenCanPack()
    .with_tolerance(2.0)
    .with_seed(42)
    .run([water], max_loops=200)
)

print(result.converged, result.natoms, result.fdist, result.frest)
packed = result.frame
```

| Field | Meaning |
|---|---|
| `converged` | Both objectives fell below the packer precision threshold |
| `fdist` | Pair-distance (overlap) violations |
| `frest` | Restraint violations |
| `frame` | Topology-complete packed `molrs.store.Frame` |

`GenCanPack` is the rigid-body entry; `CbmcGrow` grows chains instead. Each
engine runs once — `run()` consumes it, so build a new one per pack. If you
only want the frame, take `result.frame`.

## 5. Save

```python
import molrs

molrs.io.write_pdb("water_box.pdb", packed)
# or: molrs.io.write_xyz("water_box.xyz", packed)
```

## Where next

<div class="molcrafts-manual-grid molcrafts-manual-grid--cols-3">
  <a href="concepts/">
    <strong>Concepts</strong>
    <em>Targets, restraints, and packer phases.</em>
  </a>
  <a href="cli/">
    <strong>CLI</strong>
    <em>Same job as a Packmol-style `.inp`.</em>
  </a>
  <a href="python/">
    <strong>Python</strong>
    <em>Fixed solutes, PBC, collective restraints.</em>
  </a>
  <a href="rust/">
    <strong>Rust</strong>
    <em>Embed the engine in a native crate.</em>
  </a>
  <a href="packmol_parity/">
    <strong>Packmol parity</strong>
    <em>What matches Packmol, and what does not.</em>
  </a>
</div>
