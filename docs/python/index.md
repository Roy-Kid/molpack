# Python

Packmol-grade molecular packing in Rust, with a first-class Python binding.

`molpack` arranges $N$ molecule types under geometric restraints without
overlaps. The engine is a faithful port of Packmol's GENCAN-driven algorithm
(Martínez *et al.* 2009); five canonical workloads stay pinned against Packmol
reference output.

## At a glance

```python
import molrs
from molpack import GencanPack, Target

frame = molrs.io.read_pdb("water.pdb")

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(molrs.core.Cuboid([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)

packer = GencanPack().with_tolerance(2.0).with_seed(42)
result = packer.run([water], max_loops=200)
print(f"packed {result.frame['atoms'].nrows} atoms")
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
    <em>Dense polymer melts via the CbmcGrow entry.</em>
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
- [`molcrafts-molrs`](https://pypi.org/project/molcrafts-molrs/) — frame I/O and
  the `Frame` data model (installed as a dependency)
