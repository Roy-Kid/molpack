# Chain Growth

The default packing algorithm treats every molecule as a rigid body: the
template's shape is fixed, and the optimizer searches over translations
and rotations. That model is right for water, urea, lipids, proteins —
species whose conformation is settled before packing begins. For a dense
polymer melt it is wrong, and wrong in geometry rather than in numerics:
at melt density every point of the box lies inside the pervaded volume
of many chains at once, so a valid structure consists of *interpenetrating*
coils. No downhill search over rigid placements reaches such a state —
in the PEO benchmark backing this feature, the rigid path stops
converging at roughly **a third of the real melt density** and no amount
of extra iterations recovers it. A melt must be grown into place, one
torsion at a time, inside the final box.

`CbmcGrow` is that second algorithm, and it is a peer of the default,
not a preprocessor for it: growth consumes the same radii, tolerance,
and restraints, is judged by the same `fdist` / `frest` objective, and
returns the same `State`. You choose the algorithm by choosing the
entry — `GenCanPack` places rigid bodies, `CbmcGrow` grows chains — and
every target in that call is handled by it. molpack never infers the
algorithm from the molecule and never silently falls back from one to
the other.

A pack that needs both (a polymer electrolyte, say) is two runs, written
out explicitly: grow the chains, then pack the small molecules around
the grown matrix held fixed. `Target.fixed_from(result)` is the joint —
see [Staging a mixed pack](#staging-a-mixed-pack) below.

## A complete growth pack

```python
import molrs
from molpack import CbmcGrow, Target, TorsionPrior

# The template must carry its bond graph — growth consumes the
# template's chemistry, not its shape.
chain = molrs.io.read_pdb("peo_chain.pdb")

prior = TorsionPrior.three_state_from_c_inf(5.5, 1.9106)   # PEO C∞, tetrahedral θ
peo = Target(chain, count=25)

result = (
    CbmcGrow(prior)
    .with_density(0.5)      # g/cm³ — the box is cubic, periodic, final-sized
    .with_seed(42)
    .run([peo], max_loops=60)
)
print(result.converged, result.softened)
```

The torsion prior is `CbmcGrow`'s one mandatory constructor argument;
the growth knobs (`with_trials`, `with_retract`, `with_relax`, …) are
builders on the entry, alongside the shared ones (`with_density`,
`with_seed`, `with_tolerance`, …).

A target that cannot be grown — no bond graph, fewer than 3 atoms, a
`fixed_at` placement, or no box — makes `run()` raise `ValueError`
naming the problem and suggesting `GenCanPack` where that is the right
fix.

## Staging a mixed pack

Growth and rigid placement do not mix inside one call. Run them in
order, and turn the first result into a single fixed obstacle for the
second:

```python
from molpack import CbmcGrow, GenCanPack, InsideBoxRestraint, Target

grown = CbmcGrow(prior).with_density(0.5).with_seed(42).run([peo], max_loops=60)

# Same cell for stage two — read it off the grown frame, or declare the
# lengths you sized the melt to.
cell = InsideBoxRestraint([0.0, 0.0, 0.0], [l, l, l], periodic=(True, True, True))

result = (
    GenCanPack()
    .with_seed(42)
    .run([Target.fixed_from(grown), salt.with_restraint(cell)], max_loops=200)
)
```

`Target.fixed_from(result)` wraps a whole `State` as one fixed
target whose coordinates are kept verbatim, so the second run places
only the new species and never disturbs the grown chains.

## Branched trees and rings

Build the chemistry with molrs (SMILES + conformer) and molpy
`PolymerBuilder` — do not invent coordinates. A 4-arm star is a
tetrafunctional core plus EO arms (`build_star`); a macrocycle is
`build_ring`. Both `CbmcGrow` and `LatticeGrow` consume the **bond
graph**: a tree is legal for either grower; a cycle raises
`RingTemplate` on both, so `pack_ring` then picks rigid `GenCanPack`.

`python/examples/pack_peo_topo.py` `pack_star` is one explicit pick:
`LatticeGrow` @ 2.0 Å (occupancy guard on) then caller-side
`GenCanPack.with_restart` @ 2.0 Å — the same lattice-then-push-off
shape as the melt sample below. `CbmcGrow` remains a peer tree grower;
its reduced-EV (0.6 Å) then 2.0 Å push-off recipe stays on the
`CbmcGrow` path and is not copied onto `LatticeGrow`.

```bash
python python/examples/pack_peo_linear.py 8 8 0.5 42
python python/examples/pack_peo_mix.py 4 2 4 4 0.5 42
python python/examples/pack_peo_topo.py star 4 8 0.5 42
python python/examples/pack_peo_topo.py ring 6 8 0.4 42
python python/examples/pack_peo_stl.py 4 4 30 42
```

`pack_peo_mix.py` puts two topologies in **one** `LatticeGrow.run` (linear
`Target` + 4-arm star `Target`, density-sized box). `pack_peo_stl.py` is
the mesh-cavity scene: attach `StlRegion.from_file` and grow with
`LatticeGrow` — diamond sites outside the mesh are blocked
(Region ∩ lattice), then `GenCanPack.with_restart`.

## The torsion prior is load-bearing

`CbmcGrow` has exactly one mandatory argument, and it is the one that
decides the physics. Torsion angles are the free variables of growth;
whatever distribution they are drawn from becomes the chain statistics
of the product. Sampling them uniformly gives the freely-rotating chain,
whose characteristic ratio is C∞ = 2.0 — but real PEO melts have
C∞ ≈ 5.5, and since chain size scales as √C∞, uniform sampling
undershoots the melt radius of gyration by about 40%. Crowding at melt
density does not repair this: excluded-volume screening corrects the
scaling exponent, not the prefactor. That is why the prior has no
default — `TorsionPrior.uniform()` exists as a negative control, and
using it must be a visible decision.

The practical calibration is a single scalar. Given a target
characteristic ratio (from the literature or from ⟨R²⟩₀/M) and the
backbone valence angle, `TorsionPrior.three_state_from_c_inf(c_inf,
theta_rad)` solves the trans fraction of a trans/gauche± three-state
prior in closed form — for PEO with tetrahedral angles,
`three_state_from_c_inf(5.5, 1.9106)` yields a trans fraction of about
0.645. When you have a measured torsion distribution instead, pass it
directly as `TorsionPrior.states([(angle_rad, weight), ...])`; and
`TorsionPrior.template(kappa)` spreads draws around the template's own
torsion values. Priors are geometric data, never force-field objects —
if you want force-field-derived weights, compute them outside and pass
the numbers in.

## Density-sized boxes

Growth needs the box to be at its final volume from the first atom:
there is no compression stage, so the target density must hold from the
start. Rather than computing box lengths by hand,
`CbmcGrow.with_density(rho)` sizes a cubic, fully periodic box from the
total mass of **all** targets at `rho` g/cm³. Masses resolve from
element symbols; targets without usable elements (coarse-grained beads)
declare theirs with `Target.with_mass(amu)`. A density combined with an
explicit `with_periodic_box`, or a mass that cannot be resolved, raises
`ValueError` — molpack does not guess. An explicit periodic box works
too, if you prefer to control the geometry yourself.

## All-atom versus coarse-grained templates

Two per-target settings encode the difference between atomistic and CG
chains; both default to the all-atom convention. Mixed AA/CG packs are
one run — each target carries its own table and its own angle prior.

Bond angles: in an all-atom chain they are stiff coordinates, so the
default `AnglePrior.template()` copies them verbatim from the template.
In a CG model the angle potential is soft — it is usually *fitted* to
reproduce a persistence length — so copying template angles would freeze
the chain stiffness at whatever the template happened to be. CG chains
instead use `AnglePrior.wlc_from_c_inf(c_inf)`, a worm-like-chain angle
prior calibrated from the target characteristic ratio (for Kremer–Grest
melts, `wlc_from_c_inf(1.76)`). The calibration assumes uniform
torsions; combining a WLC angle prior with a non-uniform torsion prior
double-counts stiffness.

Intramolecular skip table: atom pairs close along the chain must be
exempt from the hard core, because their distances are governed by
bonds, angles, and the torsion prior — a hard core applied to 1-4 pairs
would reject every gauche state. The table lives on the target, not the
engine: `Target.with_special_bonds([0, 0, 0, 1])` is the default
(Cassandra depth 3: 1-2/1-3/1-4 exempt, 1-5+ scored). CG templates
conventionally use a shallower table — `[0, 1]` (depth 1) or
`[0, 0, 1]` (depth 2) — on that target only.

All-atom chains with explicit hydrogen keep the default depth-3 table
and shrink hydrogen with `Target.with_atom_radius`. Hydrogen's packing
radius is a first-class per-atom setting, not a reason to deepen the
skip table. On the 2026-09-04 dp5 PEO melt, depth 3 plus H = 0.85 Å
finished in 142 rounds / 0.3 s.

```python
h = [i for i, e in enumerate(chain["atoms"]["element"]) if e == "H"]
aa = Target(chain, count=25).with_atom_radius(h, 0.85)   # default [0, 0, 0, 1]
cg = Target(beads, count=40).with_special_bonds([0, 0, 1])  # depth 2
result = CbmcGrow(prior).with_density(0.5).run([aa, cg], max_loops=60)
```

## Reading `softened`

Growth's guarantee is constructive: a candidate placement that violates
the hard core or a restraint is rejected, never penalized, so a
successfully grown structure has `fdist == 0` by construction rather
than by convergence. When a region of the box becomes so crowded that a
chain dead-ends repeatedly even after retracting and regrowing, the
solver's last resort is to shrink the hard core — and every one of those
shrinks increments `State.softened`. Each unit therefore records
one relaxation of the constructive guarantee; a grown structure only
counts as converged when the count is zero at full tolerance, and on the
rigid-body path it is always zero.

Softening is not a dead end, but the remedy is something you ask for —
an explicit second stage, never a hidden fallback. Feed the **same free
targets** to a seeded `GenCanPack`:

```python
grown = CbmcGrow(prior).with_density(0.9).run([chain], max_loops=60)
pushed = GenCanPack().with_restart(grown).with_seed(7).run([chain], max_loops=60)
```

The seeded run continues on the very same state — zero coordinate
conversion — through the GENCAN phases, driving the remaining contact
violations out by rigid-body descent (the classic slow push-off). Each
link reports honestly: the grow result keeps its `softened` count so you
can see the guarantee was relaxed, and the seeded run's `converged`
tells you whether the push-off restored the full tolerance.

## Reading `intra`

`State.fdist` is intermolecular only. Same-copy contacts are
classified into `result.intra.scored` and `result.intra.exempted` (Å,
minimum image; an empty class is `+∞`). Exempted pairs are the ones the
target's special-bonds table skipped (1-2/1-3/1-4 at the default);
scored pairs are intramolecular contacts the hard core was supposed to
keep.

If scored intramolecular contacts dominate after a depth-3 all-atom run
with shrunk hydrogens, the residual is the signal to deepen that
target's table as an escape hatch — `with_special_bonds([0, 0, 0, 0, 0, 1])`
exempts out to 1-6 — not a reason to change the engine.

## What growth delivers — and what it does not

The product is a **geometric starting structure**: no contacts below the
declared tolerance, chain statistics governed by your priors, and
homogeneous density. It is not an equilibrium Boltzmann ensemble —
greedy Rosenbluth selection has a known, mild bias (grown chains come
out slightly compact in crowded systems), and equilibration belongs to
the downstream MD following the classic generate → push-off →
equilibrate pipeline. Validate the product by measuring it: Rg against
the unperturbed value, the internal-distance curve, density uniformity —
not by expecting equilibrium statistics from a constructor.

## Melt density: the lattice grower

Past ρ ≈ 0.5 the continuum grower grinds: its candidates are proposed in
continuous space and die by retraction. `LatticeGrow` maps the box onto a
diamond lattice where the same trans/gauche± RIS states are *exact* lattice
moves and excluded volume is an O(1) site check — a melt-density system
generates in milliseconds, then every atom is rebuilt from the template's
exact internal coordinates:

```python
grown = (
    LatticeGrow(prior)          # same mandatory torsion prior
    .with_seed(42)
    .with_density(1.1)
    .run([Target(frame, 200)], max_loops=60)
)
pushed = GenCanPack().with_restart(grown).with_seed(42).run([Target(frame, 200)], max_loops=200)
```

The lattice decides only the torsion sequence; bond lengths and angles are
the template's, bit-exact. Residual contacts (hydrogens, decoration drift)
are reported honestly in `fdist` and belong to the seeded push-off.

This repository's polymer-melt benchmark claim cap is ρ = 1.2 g/cm³; the
lattice sample above uses 1.1. `with_density` itself has no algorithm
upper bound.
