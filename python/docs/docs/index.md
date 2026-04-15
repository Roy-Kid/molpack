# molpack — Python

Packmol-grade molecular packing in Rust, with Python bindings.

`molpack` produces a non-overlapping arrangement of $N$ molecule types
with user-specified copy counts under geometric restraints. The engine
is a faithful port of Packmol's GENCAN-driven three-phase algorithm
(Martínez *et al.* 2009); correctness is pinned against Packmol's
reference output for five canonical workloads.

## At a glance

```python
import molrs
from molpack import InsideBox, Molpack, Target

frame = molrs.read_pdb("water.pdb")

water = (
    Target("water", frame, count=100)
    .with_restraint(InsideBox([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)

result = Molpack(tolerance=2.0).pack([water], max_loops=200, seed=42)
print(f"converged={result.converged}  natoms={result.natoms}  fdist={result.fdist:.4f}")
```

## Where next

- [Installation](installation.md) — pip install and verification.
- [Getting Started](getting-started.md) — first pack end-to-end.
- User Guide:
  [Targets](guide/targets.md) ·
  [Restraints](guide/restraints.md) ·
  [Packer](guide/packer.md) ·
  [Periodic boundaries](guide/periodic-boundaries.md)
- [Examples](examples.md) — five Packmol-equivalent workloads.
- [API Reference](api-reference.md) — class-by-class summary.

## Related

- Rust crate: [`molcrafts-molpack`](https://crates.io/crates/molcrafts-molpack)
  — the underlying engine. All algorithmic details are documented there.
- [`molcrafts-molrs`](https://pypi.org/project/molcrafts-molrs/) —
  companion package for file I/O (PDB, XYZ, …) and the `Frame` data
  model. Pass a `molrs.Frame` directly to `Target`; no manual array
  extraction needed.

[targets]: guide/targets.md
