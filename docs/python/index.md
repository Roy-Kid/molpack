# Python

Packmol-grade molecular packing in Rust, with a first-class Python binding.

`molpack` arranges $N$ molecule types under geometric restraints without
overlaps. The engine is a faithful port of Packmol's GENCAN-driven algorithm
(Martínez *et al.* 2009); five canonical workloads stay pinned against Packmol
reference output.

## At a glance

```python
from pathlib import Path

import molpy
from molpack import GencanPack, Target

frame = molpy.io.read_pdb(Path("water.pdb"))

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(molpy.core.Cuboid([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)

packer = GencanPack().with_tolerance(2.0).with_seed(42)
result = packer.run([water], max_loops=200)
print(f"packed {result.frame['atoms'].n_rows} atoms")
```

## Next steps

<div class="molcrafts-manual-list">
  <a href="installation/">
    <strong>Installation</strong>
    <em>pip install and verification.</em>
  </a>
  <a href="getting-started/">
    <strong>Quickstart</strong>
    <em>First pack end-to-end.</em>
  </a>
  <a href="guide/targets/">
    <strong>Targets</strong>
    <em>Templates, counts, fixed placement.</em>
  </a>
  <a href="guide/restraints/">
    <strong>Restraints</strong>
    <em>Geometric and collective restraints.</em>
  </a>
  <a href="guide/packer/">
    <strong>Packer</strong>
    <em>Builder options, diagnostics, callbacks.</em>
  </a>
  <a href="guide/growth/">
    <strong>Chain growth</strong>
    <em>Dense polymer melts via the CbmcGrow engine.</em>
  </a>
  <a href="guide/periodic-boundaries/">
    <strong>Periodic boundaries</strong>
    <em>PBC rules and validation.</em>
  </a>
  <a href="examples/">
    <strong>Examples</strong>
    <em>Packmol-equivalent workloads.</em>
  </a>
  <a href="api-reference/">
    <strong>API reference</strong>
    <em>Class-by-class summary.</em>
  </a>
</div>

## See also

- Rust crate: [`molcrafts-molpack`](https://crates.io/crates/molcrafts-molpack) —
  the underlying engine
- [`molcrafts-molpy`](https://pypi.org/project/molcrafts-molpy/) — frame I/O and
  the `Frame` data model (installed as a dependency)
