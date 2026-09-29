# Restraints

A restraint is a soft penalty every atom of a target — or a chosen subset —
pays when it leaves where it should be. molpack knows two families: the
**region lift**, one geometric restraint that says "stay inside this molrs
region", and seven **collective** restraints — six distribution-matching
([further down](#collective-distribution-matching-restraints)) plus
`SelfSeparation` ([after those](#collective-separation-restraint)). All
attach the same way, via `.with_restraint()`.

## Regions are molrs objects

Geometry is not molpack's. A region is a molrs solid with a signed distance
to its boundary — `distance(points)` is negative inside, positive outside —
and every shape describes its *inside*. Outside, shells and voids are
compositions: `~`, `&`, `|`. There is no "outside sphere" class; it is
`~molrs.Sphere(...)`.

| molrs class     | Constructor                                   | Meaning |
|-----------------|-----------------------------------------------|---------|
| `Sphere`        | `(center, radius)`                            | solid sphere |
| `Cuboid`        | `(origin, lengths)`                           | axis-aligned box from its minimum corner |
| `Parallelepiped`| `(h, origin)`, `.cube(a, origin)`, `.ortho(lengths, origin)` | triclinic cell volume |
| `HalfSpace`     | `(normal, point)`                             | the side of the plane the normal points *away* from |
| `Cylinder`      | `(base, axis, radius, length)`                | finite capped cylinder |
| `Ellipsoid`     | `(center, semi_axes)`                         | axis-aligned ellipsoid |
| `Polyhedron`    | `(mesh: TriMesh)`                             | solid bounded by a watertight triangle mesh |
| `SphereUnion`   | `(centers, radii, box=None)`                  | union of spheres, minimum image on the box's periodic axes |

Length-3 arguments accept plain lists. Lengths are in the same unit as the
coordinates (Å for a molpack run).

```python
import molrs
from molpack import Target

box    = molrs.Cuboid([0, 0, 0], [40, 40, 40])          # origin, lengths
ball   = molrs.Sphere([0, 0, 0], 20.0)
shell  = ball & ~molrs.Sphere([0, 0, 0], 10.0)
above  = ~molrs.HalfSpace([0, 0, 1], [0, 0, 5.0])        # z >= 5
below  = molrs.HalfSpace([0, 0, 1], [0, 0, 20.0])        # z <= 20
water  = Target(frame, count=500).with_restraint(box)
```

molpack lifts the region to a quadratic exterior penalty,
`scale · max(0, distance)²`, which is zero inside and on the boundary. The
test is the **atom centre**, not the van der Waals ball and not the molecule
COM. A region reaches this wheel as a `molrs.RegionRef` capsule — no data is
marshalled, and both wheels must share one molrs minor line.

How each packing entry uses a region is **not** inferred from the shape:

- `GenCanPack` — soft quadratic wall on atom centres (`frest`).
- `CbmcGrow` — hard reject on `propose`; `force_place` may leave atoms
  outside and counts `degraded`.
- `LatticeGrow` — diamond sites outside the region are blocked
  (Region ∩ lattice). An empty intersection is a named error. Decorated
  hydrogens may still sit slightly outside; chain `GenCanPack.with_restart`.

### A cavity from a mesh

`molrs.io.read_stl` reads an ASCII or binary STL into a `TriMesh`;
`TriMesh.scaled` converts the file's unit; `molrs.Polyhedron` is the solid
the mesh bounds. The mesh must be watertight (every edge shared by exactly
two faces) — an open or self-touching mesh is a `ValueError`, because
parity cannot decide inside from outside on it.

```python
cavity = molrs.Polyhedron(molrs.io.read_stl("cavity.stl").scaled(4.18))
target = Target(frame, n).with_restraint(cavity)
```

`python/examples/pack_peo_mesh.py` grows 200 × EO25 with `LatticeGrow`
inside the shipped dendrite `examples/pack_peo/dendrite.stl` and stops
after the grow: at melt density inside a cavity a rigid-body push-off can
only resolve overlap through the wall.

### A void from atoms

The solvent-accessible void of a bead cloud needs no mesh. One sphere per
bead of radius `bead radius + probe radius` is the solvent-accessible
volume; its complement is where a probe centre may go:

```python
polymer = molrs.SphereUnion(centers, 0.5 * sigma + 1.0, box=frame.box)
void = ~polymer
target = Target(peo, n).with_restraint(void)
```

`python/examples/pack_peo_void.py` does this for a LAMMPS data file: the
bonded atoms are the polymer, the rest is solvent and is dropped, and PEO
threads the channels the solvent left.

## Periodic boxes

The periodic box is the entry's declaration, `with_periodic_box(min, max)`;
a region only confines. See
[Periodic boundaries](periodic-boundaries.md) for the full semantics and
validation rules.

## Collective (distribution-matching) restraints

Where a region penalises **each atom** against a solid, a collective
restraint sees **every copy of the target at once** and drives the
species' spatial *distribution* toward a target profile (via a squared 1-D
Wasserstein penalty). Six built-ins cover Gaussian, exponential, and
arbitrary-tabulated priors along either a plane or a radius — e.g. a
Gaussian slab centred at `z = 20`:

```python
from molpack import GaussianPlane

slab = GaussianPlane(normal=[0, 0, 1], offset=0.0, strength=1.0, mu=20.0, sigma=3.0)
target = Target(frame, count=200).with_restraint(slab)
```

The full list and constructor signatures are in the
[API reference](../api-reference.md#collective-distribution-matching-restraints).

## Collective (separation) restraint

The pair term keeps atoms from overlapping, but nothing in it stops a
species from piling all its copies into one corner — it cannot tell two
copies of one species from a copy of each of two. `SelfSeparation` is the
missing bound: no two copies closer than `d_min`, centre to centre.

```python
import molrs
from molpack import SelfSeparation, Target

ions = (
    Target(frame, count=27)
    .with_restraint(molrs.Cuboid([0, 0, 0], [40, 40, 40]))
    .with_restraint(SelfSeparation(10.0))  # ions stay 10 Å apart
)
```

Distances use the minimum image, so copies on opposite sides of a periodic
boundary are pushed apart like any other neighbours. Nothing checks that
the request fits: if it cannot, the run does not converge and says so
through `frest`.

## Stacking multiple restraints

Apply several restraints to the same target by chaining `.with_restraint()`,
or compose the regions first — the two are equivalent for regions:

```python
target = (
    Target(frame, count=500)
    .with_name("water")
    .with_restraint(molrs.Cuboid([0, 0, 0], [40, 40, 40]) & ~molrs.Sphere([20, 20, 20], 5.0))
)
```

Each call attaches an independent restraint. All active restraints are
evaluated at every optimizer step.

## Scopes

A restraint can be applied at two scopes:

- **Whole target** — `target.with_restraint(r)` — penalises every atom.
- **Atom subset** — `target.with_atom_restraint([0, 1, 2], r)` —
  penalises only the listed atoms (0-based indices).

Example — a bilayer: pin heads above z=12, tails below z=2:

```python
lipid = (
    Target(frame, count=20)
    .with_name("lipid")
    .with_restraint(molrs.Cuboid([0, 0, 0], [40, 40, 14]))
    .with_atom_restraint([0, 1],   ~molrs.HalfSpace([0, 0, 1], [0, 0, 12.0]))
    .with_atom_restraint([30, 31],  molrs.HalfSpace([0, 0, 1], [0, 0, 2.0]))
)
```

## Global restraints

To apply one restraint to every target in a pack, attach it on the
engine entry:

```python
packer = (
    GenCanPack()
    .with_global_restraint(molrs.Cuboid([0, 0, 0], [40, 40, 40]))
)
```

Semantically equivalent to calling `.with_restraint(r)` on every
target, but avoids the duplication.

## Custom restraints

Pass any object implementing `f(x, scale, scale2) -> float` and
`fg(x, scale, scale2) -> (float, (gx, gy, gz))`. See the
`Restraint` Protocol in `molpack` for the full contract. Reach for this
when the penalty is not "stay inside a region"; a new *shape* is a molrs
region, not a custom restraint.

```python
class SphereRestraint:
    def __init__(self, center, radius):
        self.c = np.asarray(center)
        self.r = radius
    def f(self, x, scale, scale2):
        d = np.linalg.norm(np.asarray(x) - self.c) - self.r
        return scale2 * d * d if d > 0 else 0.0
    def fg(self, x, scale, scale2):
        rel = np.asarray(x) - self.c
        d = float(np.linalg.norm(rel))
        over = d - self.r
        if over <= 0:
            return 0.0, (0.0, 0.0, 0.0)
        factor = 2 * scale2 * over / d
        return scale2 * over * over, tuple(factor * rel)

target = Target(frame, count=10).with_restraint(SphereRestraint([0,0,0], 6.0))
```

## Semantics

Every restraint contributes a continuously differentiable penalty
$f_{\text{rest}}(\mathbf{x})$ that is zero inside the allowed region
and rises quadratically outside. The aggregate objective minimised by
the packer is:

$$
U(\mathbf{x}) = f_{\text{dist}}(\mathbf{x}) + f_{\text{rest}}(\mathbf{x})
$$

where $f_{\text{dist}}$ is the pairwise distance-violation sum for the
user-specified `tolerance`. Convergence is declared when both fall
below `precision`.

!!! note "Restraints vs hard constraints"
    All built-in restraints are *soft penalties* — the optimizer may
    momentarily produce a violating configuration while searching. Hard
    geometric constraints (frozen placement, rotation bounds) are set
    on the `Target` directly via `fixed_at` and `with_rotation_bound`.
