# Python Guide

If you are using `molpack` from Python, you are still using the same Rust
engine. The difference is the feel: the Python binding is friendlier to scripts,
notebooks, and frame-oriented data pipelines.

## What This Section Is For

Use the Python pages when you want:

- a quick path from `pip install` to a first packed system;
- a clearer mental model for `Target`, `Molpack`, and the built-in restraints;
- a binding-focused reference without leaving the main docs site.

## The Shape of a Python Packing Job

Most Python jobs look like this:

```python
import molrs
from molpack import InsideBoxRestraint, Molpack, Target

frame = molrs.read_pdb("water.pdb")

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)

result = Molpack().with_tolerance(2.0).with_seed(42).pack([water], max_loops=200)
print(result.converged, result.fdist, result.frest)
```

The pattern stays the same no matter how fancy the system gets:

1. load or build a frame;
2. wrap it in one or more `Target` objects;
3. attach restraints;
4. run the packer;
5. inspect diagnostics before trusting the coordinates.

## Where to Go Next

- [Install](../install.md) for package setup and source builds.
- [Getting Started](../getting_started.md) for the first end-to-end run.
- [Targets](targets.md) for counts, centering, fixed placements, and rotation bounds.
- [Restraints](restraints.md) for whole-target, atom-subset, and custom geometry.
- [Packer](packer.md) for tuning knobs, handlers, and result diagnostics.
- [Periodic Boundaries](periodic_boundaries.md) for orthorhombic PBC rules.
- [API Reference](api_reference.md) for a binding-focused symbol map.

## Takeaway

The Python binding is not a separate design. It is the same packing model with a
lighter surface. Once you learn the `Target` plus `Molpack` split, most of the
rest falls into place quickly.
