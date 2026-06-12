# Extensibility case: a measured head-group profile, prototyped in Python, run in Rust

This note designs a real-research packing problem around molpack's new
**profile-distribution restraints** and tells it as the paper's extensibility
arc: *someone wants something the input grammar can't express, prototypes it in
Python in minutes, and — once the physics is right — runs the same idea as
native Rust with no change to the driver.* It is the restraint-side companion to
Section 4's relaxer (chain-folding) case, and is backed by two runnable scripts
in this directory:

- [`pack_profile_monolayer.py`](./pack_profile_monolayer.py) — the prototype.
- [`profile_speed_py_vs_rust.py`](./profile_speed_py_vs_rust.py) — the Python-vs-Rust timing.

## 1. The problem the grammar can't express

Packmol-style restraints describe **hard regions** (inside/outside a box,
sphere, plane …) and therefore fill them **uniformly**. A common interfacial
task does not fit that mould: you have measured, from neutron or X-ray
reflectometry, *where* a surfactant's polar head group sits relative to the
water surface — a Gaussian band centred at depth `μ` with width `σ` — and you
want an MD **starting** configuration whose heads already follow that profile.
A bare `inside box` spreads the heads uniformly through the slab, so the run
then burns a long equilibration pulling them into place. The same shape need
recurs across asymmetric bilayer leaflets, electrode–electrolyte double layers
(exponential decay), and polymer-blend interfaces (erf/tanh steps).

molpack's **profile-distribution restraints** express exactly this. A restraint
biases a selected site toward a target distribution `ρ*(ξ)` by Boltzmann
inversion,

```
U(ξ) = −kT · ln( ρ*(ξ) / ρ₀ ),
```

composing a **geometric coordinate** `ξ` (planar / radial / cylindrical /
region-distance) with a **distribution** `ρ*` (Gaussian / erf / tanh /
exponential / tabulated). The Gaussian is the simplest member: inverting
`ρ*(z) ∝ exp(−(z−μ)²/2σ²)` gives the harmonic well
`U(z) = (kT/2σ²)(z−μ)²` — a spring of stiffness `kT/σ²` along the plane normal.

## 2. Prototype in Python (minutes, no recompile)

A practitioner who wants "heads on a Gaussian" writes the restraint as an
ordinary Python object exposing `f` (energy) and `fg` (energy + gradient), and
drops it into the **same slot** a built-in restraint occupies — there is no
plugin wrapper. The whole restraint is ten lines:

```python
class PlanarGaussian:
    def __init__(self, normal, point, mu, sigma, kt=1.0):
        n = np.asarray(normal, float); self.n = n / np.linalg.norm(n)
        self.x0 = np.asarray(point, float)
        self.mu = mu; self.k = kt / (sigma * sigma)       # spring kT/σ²
    def _xi(self, x):
        return float(self.n @ (np.asarray(x, float) - self.x0))
    def f(self, x, scale, scale2):
        d = self._xi(x) - self.mu
        return scale * 0.5 * self.k * d * d
    def fg(self, x, scale, scale2):
        d = self._xi(x) - self.mu
        g = scale * self.k * d * self.n
        return scale * 0.5 * self.k * d * d, (float(g[0]), float(g[1]), float(g[2]))
```

Attaching it to the head site (a *subset* of atoms — biasing the whole molecule
would distort it) and packing is the normal API:

```python
target = (molpack.Target(frame, count=64).with_name("surfactant")
          .with_restraint(molpack.InsideBoxRestraint([0,0,0], [40,40,30]))
          .with_atom_restraint([30, 31], PlanarGaussian([0,0,1], [0,0,0], mu=6, sigma=2)))
packed = molpack.Molpack().with_seed(20240612).pack_with_report([target], max_loops=400)
```

Packing 64 palmitoyl chains (`palmitoil.pdb`, heads = atoms 31–32) with this
prototype reproduces the target band:

```
target  : mu=6.00 Å  sigma=2.00 Å
realised: mean=6.01 Å  std=0.40 Å  (n=128 head sites)   100% of heads within mu±2σ
```

The head **centre** lands on `μ` to 0.01 Å. The realised width is *tighter* than
`σ` because packing is energy **minimisation**, not thermal sampling — a finite-
temperature MD run broadens it back toward `σ`; the restraint sets where the
band sits, the ensemble sets its width. This is the spec's "soft, not exact"
property (§6.6): a soft bias *approaches* `ρ*` while competing with the overlap
term.

## 3. Mature, then sink to Rust (same driver, native speed)

Once the prototype has earned its place, the identical physics is already a
first-class **native** restraint. No new Python class is needed: the same
program reaches it through the `profile` script keyword, with the head bias a
single line inside the existing `atoms … end atoms` sub-block:

```
structure palmitoil.pdb
  number 64
  inside box 0 0 0 40 40 30
  atoms 31 32
    profile gaussian plane 0. 0. 1.  0. 0. 0.  mu 6. sigma 2. density
  end atoms
end structure
```

Driven from Python with `molpack.load_script(...)`, the restraint now runs
inline in the compiled hot loop — no Python callback per evaluation. molpack's
native lowering uses `kT = 1.0`, so the Rust `U(z) = (kT/2σ²)(z−μ)²` is the
*same* harmonic well the prototype computes; hardening into Rust additionally
gives a density floor / energy cap (§6.3, `U_max = −kT·ln(ρ_min/ρ₀) ≈ 13.8`) and
the shell-volume Jacobian for radial/cylindrical geometries (§6.2) that the
minimal prototype omits — exactly the robustness a one-off script skips and a
production restraint must have.

## 4. The speed gap

The same restraint, both ways, driven from the same Python program — same
structure / count / box / `μ`,`σ` / seed / tolerance / `max_loops`, the only
difference being *where the restraint executes*. The Python prototype pays a
per-site Python call (GIL + marshalling) on every objective and gradient
evaluation; the native keyword runs inline.

How visible that overhead is depends on *what fraction of the work is the
restraint*. On the monolayer above the bias touches only 2 of 42 atoms per
chain, so it is a sliver of each evaluation — the shared overlap kernel
dominates both paths and the prototype is already within ~1.4× of native. The
honest way to isolate the restraint cost is a system where the biased site **is**
the molecule: an electrode-style diffuse layer of **single-atom ions**, each
biased toward a Gaussian band at `μ` from the plane, so the restraint is
evaluated for every atom. `profile_speed_py_vs_rust.py` runs exactly this and
measures (Apple M-series, median of 3, `max_loops = 50`, 25×25×40 Å box):

| ions | Python prototype (s) | native Rust (s) | speed-up | realised band (mean±std, Å) |
|-----:|---------------------:|----------------:|---------:|:----------------------------|
|  250 | 15.9                 | 1.8             | **8.8×** | py 20.00±0.73 / rs 19.99±0.66 |
|  500 | 51.0                 | 9.3             | **5.5×** | py 20.00±1.39 / rs 19.98±1.37 |
| 1000 | 108.0                | 8.9             | **12.1×**| py 20.01±2.64 / rs 19.98±2.60 |

The two paths realise the **same** ion band (functional parity, not bitwise —
same restraint math, separate random streams: every mean lands on `μ = 20` and
the widths agree to ~0.03 Å). The Python prototype is consistently **5–12×
slower**; the ratio is non-monotonic because the optimizer's stochastic path
(`movebad` relocation + convergence) varies run to run, but the per-atom Python
callback overhead is always several-fold. The lesson for the user is the useful
one: prototype a restraint in Python for free when it touches a few sites; when
it touches every atom, the same code sunk to Rust is several-fold faster — and
the driver script never changes.

## 5. Why this belongs in the paper

It is the restraint-interface twin of the relaxer case: a constraint the input
grammar cannot express, written in a few lines of Python that occupy the same
slot as a built-in, validated against the realised structure, and then run as
native code with the driver untouched. The profile-distribution family turns
"match this measured distribution" — reflectometry profiles, asymmetric
leaflets, double layers — into a single programmable restraint, packed in one
pass, and the prototype-to-production path costs the user nothing but deleting
their Python class.
