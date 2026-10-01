# Extending the Crate

Tutorials for writing your own `AtomRestraint` / `Region` / `Handler` types,
plus the in-loop optimizer seam and the `Stage` seam every packing algorithm
implements. Every extension trait in this crate follows the same shape
(direction-3 rule — see [`concepts`](crate::concepts)):

> `pub trait X` + N concrete `pub struct` types implementing it.
> User types `impl X` identically. No built-in/plugin type-level
> distinction.

## Custom `AtomRestraint`

Goal: pull atoms toward a target plane with a quadratic attractive
well.

### Step 1 — define the struct

```rust
use molrs::types::F;
# use molpack::AtomRestraint;

#[derive(Debug, Clone, Copy)]
pub struct PlaneTether {
    pub normal: [F; 3],
    pub offset: F,
    pub k: F,
}
# impl AtomRestraint for PlaneTether {
#     fn f(&self, _x: &[F; 3], _s: F, _s2: F) -> F { 0.0 }
#     fn fg(&self, _x: &[F; 3], _s: F, _s2: F, _g: &mut [F; 3]) -> F { 0.0 }
# }
```

- `pub` fields — users construct with `PlaneTether { normal, offset,
  k }`, no builder.
- `Debug` required because [`AtomRestraint`](crate::AtomRestraint) has a
  `Debug` supertrait bound (so `Target`'s derived `Debug` keeps
  working).

### Step 2 — implement `AtomRestraint`

```rust
# use molrs::types::F;
# use molpack::AtomRestraint;
# #[derive(Debug)]
# pub struct PlaneTether { pub normal: [F; 3], pub offset: F, pub k: F }
impl AtomRestraint for PlaneTether {
    fn f(&self, pos: &[F; 3], _scale: F, _scale2: F) -> F {
        let d = self.normal[0] * pos[0]
              + self.normal[1] * pos[1]
              + self.normal[2] * pos[2]
              - self.offset;
        0.5 * self.k * d * d
    }
    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let d = self.normal[0] * pos[0]
              + self.normal[1] * pos[1]
              + self.normal[2] * pos[2]
              - self.offset;
        g[0] += self.k * d * self.normal[0];
        g[1] += self.k * d * self.normal[1];
        g[2] += self.k * d * self.normal[2];
        self.f(pos, scale, scale2)
    }
}
```

Three contracts that all restraints must obey:

1. **Gradient accumulates with `+=`.** Multiple restraints may touch
   the same atom.
2. **`fg` returns the value.** The hot path uses the returned value
   for the `fdist`/`frest` accumulation — don't return `0.0` just
   because the caller might discard it.
3. **Scale/scale2 usage is your choice.** Linear-penalty restraints
   typically use `scale`; quadratic-penalty ones use `scale2`. Your
   tether uses its own `k` — ignore both knobs if you prefer.

### Step 3 — write a gradient test

```no_run
# use molrs::types::F;
# use molpack::AtomRestraint;
# #[derive(Debug)] pub struct PlaneTether { pub normal: [F; 3], pub offset: F, pub k: F }
# impl AtomRestraint for PlaneTether {
#     fn f(&self, _x: &[F; 3], _s: F, _s2: F) -> F { 0.0 }
#     fn fg(&self, _x: &[F; 3], _s: F, _s2: F, _g: &mut [F; 3]) -> F { 0.0 }
# }
#[test]
fn plane_tether_gradient_matches_fd() {
    let r = PlaneTether { normal: [0.0, 0.0, 1.0], offset: 5.0, k: 2.0 };
    let x = [1.0, 2.0, 7.0];
    let mut g = [0.0; 3];
    let _ = r.fg(&x, 1.0, 1.0, &mut g);
    let h: F = 1e-5;
    for k in 0..3 {
        let mut xp = x; xp[k] += h;
        let mut xm = x; xm[k] -= h;
        let fd = (r.f(&xp, 1.0, 1.0) - r.f(&xm, 1.0, 1.0)) / (2.0 * h);
        assert!(
            (g[k] - fd).abs() < 1e-4,
            "axis {k}: analytic={}, fd={}", g[k], fd,
        );
    }
}
```

Convention: `ε = 1e-5`, tolerance `1e-3` (looser if the restraint has
kinks).

### Step 4 — use it

```no_run
# use molrs::types::F;
# use molpack::AtomRestraint;
# #[derive(Debug, Clone, Copy)] pub struct PlaneTether { pub normal: [F; 3], pub offset: F, pub k: F }
# impl AtomRestraint for PlaneTether {
#     fn f(&self, _x: &[F; 3], _s: F, _s2: F) -> F { 0.0 }
#     fn fg(&self, _x: &[F; 3], _s: F, _s2: F, _g: &mut [F; 3]) -> F { 0.0 }
# }
use std::sync::Arc;
use molpack::{RegionRestraint, Target};
use molrs::spatial::region::Cuboid;
use ndarray::array;
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);

let cube = Cuboid::new(array![0.0, 0.0, 0.0], array![40.0, 40.0, 40.0]);
let target = Target::from_coords(pos, rad, 100)
    .with_restraint(RegionRestraint(Arc::new(cube)))
    .with_restraint(PlaneTether { normal: [0.0, 0.0, 1.0], offset: 20.0, k: 1.0 });
```

The region lift `RegionRestraint` and the user `PlaneTether` take the same
code path — direction-3 in action.

### In Python — the same restraint, duck-typed

The interfaces are also exposed to Python through a *duck-typed* protocol: a
restraint is any object exposing `f` and `fg`, and the packer consumes it on the
same code path as a built-in. That makes Python the natural place to *prototype*
a restraint the input grammar can't express — a few lines, no recompile — and
then, once it earns its place, swap in the native equivalent with the driver
script unchanged.

A duck-typed restraint is any object exposing `f` and `fg`. The Python protocol
differs from the Rust trait in exactly one way: `fg` **returns**
`(energy, (gx, gy, gz))` instead of accumulating into a `&mut [F; 3]` — the
binding does the `+=` for you. The contracts are otherwise identical: `fg`
returns the energy, the gradient is the chain rule `(dU/dξ)·∇ξ`, a linear-energy
penalty rides `scale` (a quadratic one `scale2`). A missing `f` or `fg` is
rejected at attach time with a `TypeError`, and an exception raised inside `fg`
propagates back out of `run`. The `PlaneTether` above, duck-typed in Python:

```python
import numpy as np

class PlaneTether:
    """Pull a site toward the plane z = offset (a per-site geometric bias)."""
    def __init__(self, normal, offset, k):
        n = np.asarray(normal, float)
        self.n, self.offset, self.k = n / np.linalg.norm(n), offset, k

    def f(self, x, scale, scale2):
        d = float(self.n @ np.asarray(x, float)) - self.offset
        return scale * 0.5 * self.k * d * d

    def fg(self, x, scale, scale2):
        d = float(self.n @ np.asarray(x, float)) - self.offset
        g = scale * self.k * d * self.n
        return scale * 0.5 * self.k * d * d, (float(g[0]), float(g[1]), float(g[2]))
```

Attach it with `with_atom_restraint` (a site subset) or `with_restraint` (every
copy), and pack.

> **A per-site field cannot reproduce a distribution.** It is tempting to build a
> per-site penalty from a target density `ρ*` (e.g. Boltzmann inversion
> `U = −kT·ln ρ*`) to make sites *follow* `ρ*`. This fails under packing's energy
> **minimisation**: `∑ᵢ U(xᵢ)` is minimised by driving *every* site to the single
> minimum of `U` — the mode of `ρ*` — so the sites collapse onto the peak instead
> of spreading over the distribution. To drive a whole species onto a target
> profile, use a **collective** restraint (`Target.with_collective_restraint`,
> e.g. the built-in `TabulatedPlane`), whose penalty is a function of the entire
> group and whose gradient couples the copies, so the fixed point is *empirical
> distribution = target*.

## Custom `Region`

A region is molrs vocabulary: implement
[`molrs::spatial::region::Region`] and lift it with
[`RegionRestraint`](crate::RegionRestraint). Goal: a conical region with
apex at origin, axis along +z, half-angle 30°.

```rust
use molrs::spatial::region::Region;
use molrs::types::{F, FNx3};
use ndarray::Array2;

#[derive(Debug, Clone, Copy)]
pub struct Cone {
    pub apex: [F; 3],
    pub axis: [F; 3],
    pub half_angle_cos: F,
}

impl Region for Cone {
    fn bounds(&self) -> FNx3 {
        // Unbounded along the axis; ±∞ is the honest answer.
        let mut b = Array2::zeros((3, 2));
        for d in 0..3 {
            b[[d, 0]] = F::NEG_INFINITY;
            b[[d, 1]] = F::INFINITY;
        }
        b
    }
    fn distance(&self, x: &[F; 3]) -> F {
        let dx = x[0] - self.apex[0];
        let dy = x[1] - self.apex[1];
        let dz = x[2] - self.apex[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-12 { return 0.0; }
        let axis_dot =
            (dx * self.axis[0] + dy * self.axis[1] + dz * self.axis[2]) / r;
        // Not Euclidean, but the sign and the gradient direction are right,
        // which is all the lift needs.
        self.half_angle_cos - axis_dot
    }
    // The default `distance_grad` is a finite difference — fine for a
    // prototype; override analytically for hot-path use (below).
}
```

### Compose with built-ins

```no_run
# use molrs::spatial::region::Region;
# use molrs::types::{F, FNx3};
# use ndarray::Array2;
# #[derive(Debug, Clone, Copy)]
# pub struct Cone { pub apex: [F; 3], pub axis: [F; 3], pub half_angle_cos: F }
# impl Region for Cone {
#     fn bounds(&self) -> FNx3 { Array2::zeros((3, 2)) }
#     fn distance(&self, _x: &[F; 3]) -> F { 0.0 }
# }
use std::sync::Arc;
use molpack::{RegionRestraint, Target};
use molrs::spatial::region::{AndRegion, Sphere};
use ndarray::array;
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);

let cone = Cone {
    apex: [0.0; 3],
    axis: [0.0, 0.0, 1.0],
    half_angle_cos: (std::f64::consts::PI / 6.0).cos(),
};
let sphere = Sphere::new(array![0.0, 0.0, 0.0], 10.0);
let region = AndRegion::new(Arc::new(cone), Arc::new(sphere));

let target = Target::from_coords(pos, rad, 100)
    .with_restraint(RegionRestraint(Arc::new(region)));
```

`AndRegion` / `OrRegion` / `NotRegion` take `Arc<dyn Region + Send + Sync>`:
one algebra for Rust, for Python (`&`, `|`, `~`), and for a region handed
across the wheel boundary.

### Analytic gradient override

For hot-path use, override `distance_grad` analytically. The cone above:

```rust
# use molrs::spatial::region::Region;
# use molrs::types::{F, FNx3};
# use ndarray::Array2;
# #[derive(Debug, Clone, Copy)]
# pub struct Cone { pub apex: [F; 3], pub axis: [F; 3], pub half_angle_cos: F }
# impl Region for Cone {
#     fn bounds(&self) -> FNx3 { Array2::zeros((3, 2)) }
#     fn distance(&self, _x: &[F; 3]) -> F { 0.0 }
fn distance_grad(&self, x: &[F; 3]) -> [F; 3] {
    let dx = x[0] - self.apex[0];
    let dy = x[1] - self.apex[1];
    let dz = x[2] - self.apex[2];
    let r2 = dx * dx + dy * dy + dz * dz;
    let r = r2.sqrt();
    if r < 1e-12 { return [0.0; 3]; }
    let axis_dot = (dx * self.axis[0] + dy * self.axis[1] + dz * self.axis[2]) / r;
    let inv_r = 1.0 / r;
    // distance = cos(α) - axis_dot ⇒ grad = -∂axis_dot/∂x
    [
        -(self.axis[0] * inv_r - axis_dot * dx * inv_r * inv_r),
        -(self.axis[1] * inv_r - axis_dot * dy * inv_r * inv_r),
        -(self.axis[2] * inv_r - axis_dot * dz * inv_r * inv_r),
    ]
}
# }
```

Then finite-difference check it — same pattern as the `AtomRestraint`
test.

## Custom `Handler`

Goal: a handler that writes a CSV row per step so you can plot the
objective evolution.

```no_run
use std::fs::File;
use std::io::{BufWriter, Write};
use molpack::{F, Handler, PackContext, StepInfo};

pub struct CsvHandler { writer: BufWriter<File> }

impl CsvHandler {
    pub fn new(path: &str) -> std::io::Result<Self> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "phase,loop_idx,fdist,frest,improvement_pct")?;
        Ok(Self { writer: w })
    }
}

impl Handler for CsvHandler {
    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        let _ = writeln!(
            self.writer,
            "{},{},{},{},{}",
            info.phase.phase,
            info.loop_idx,
            info.fdist,
            info.frest,
            info.improvement_pct,
        );
    }
}
```

Handler notes:

- **`on_step` is the only required method.** Everything else has a
  default no-op.
- **`sys` is `&PackContext`, never `&mut`.** Handlers cannot mutate
  packer state — bind an in-loop optimizer (next section) if you need to.
- **Multiple handlers run in registration order.** Register your CSV
  handler before `ProgressHandler` to get a row on every step,
  vice-versa otherwise.
- **`should_stop` is polled every iteration.** Return `true` to break
  the outer loop early. Useful for time budgets or custom convergence
  criteria.
- **Every step names its stage.** `info.stage` is a `StageInfo` — `index`,
  `total`, `name` — identifying the packing algorithm that emitted the step
  (stages are the subject of the `Stage` section below). A run driven by one
  engine entry has one stage, so it reports `index = 0` and `total = 1`. The
  `on_stage_start` / `on_stage_end` hooks bracket a whole stage the way
  `on_phase_start` / `on_phase_end` bracket one GENCAN phase; both are
  no-op defaults, and a single-stage run calls neither.

## Custom in-loop optimizer

Everything above places molecules as **rigid bodies**: the packer moves and
turns a copy, but the copy's internal shape — its *conformer* — stays frozen at
whatever the template said. A flexible molecule often cannot satisfy its
restraints in any single rigid pose; it has to change shape *while* it is being
placed.

That is what an **in-loop optimizer** is for. Once per outer packing iteration,
the packer hands each selected copy to an object that may rewrite that copy's
coordinates, and keeps the result only if the packing objective did not get
worse. The seam is molrs's `Optimizer` trait — one method, which relaxes a
`Frame` (molrs's atomic-data container) in place:

```text
pub trait Optimizer: Send + Sync {
    fn run(&mut self, frame: &mut Frame) -> Result<OptReport, String>;
}
```

molpack ships one implementation, `TorsionMcOptimizer`: Metropolis
Monte-Carlo sampling of rotations about a molecule's rotatable bonds.
("Metropolis Monte-Carlo" means propose a random move, then accept it with
probability `min(1, exp(-ΔE / T))` — downhill moves always pass, uphill ones
sometimes do, which is what lets the search escape a bad local shape.) molrs
ships `LBFGS`, a force-field minimizer. Anything else implementing the trait
drops into the same slot.

The trait, its molpack implementation, and the `GenCanPack::with_optimizer`
binder are always compiled. Binding a molrs force-field optimizer such as
`LBFGS` needs molrs's `ff` module, which molpack's `ff` feature forwards.

![Flexible chains packed inside a spherical cavity](assets/images/paper-confinement-sphere.png)

Confinement like this is a molpack extension workflow, not a Packmol parity
claim: an in-loop optimizer reshapes flexible chains while the packer places
them, so one engine can fit them inside a cavity no rigid pose would clear.
Two runnable programs drive `TorsionMcOptimizer` exactly this way —
`cargo run --release --example pack_adsorption` and
`cargo run --release --example pack_translocation`.

### Step 1 — implement `Optimizer`

Goal: a *jiggle* optimizer that proposes random rigid translations of the group
it is handed and keeps the ones that lower a score. A real optimizer scores
with a force field; this one just pulls the group toward the origin, so the
mechanics stay visible.

Two conventions the packer relies on:

- **Coordinates arrive and leave through the `Frame`.** Read them with
  `Frame::coords` — an `N × 3` array in ångström — and write them back with
  `Frame::set_coords`.
- **Frozen atoms are flagged.** When the binding asks for the local
  environment, the `Frame` also carries neighbouring atoms that must not move.
  They are marked by a boolean `atoms.free` column; a missing column means
  every atom is free.

```rust
use molrs::Frame;
use molrs::optimize::{OptReport, Optimizer};
use molrs::types::{F, FNx3};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

pub struct JiggleOptimizer {
    /// Trial moves per call.
    pub steps: usize,
    /// Largest displacement per axis, in ångström.
    pub max_delta: F,
    rng: SmallRng,
}

impl JiggleOptimizer {
    pub fn new(steps: usize, max_delta: F, seed: u64) -> Self {
        Self { steps, max_delta, rng: SmallRng::seed_from_u64(seed) }
    }
}

/// Toy score (Å²): squared distance of the free atoms from the origin.
fn score(coords: &FNx3, free: &[bool]) -> F {
    coords
        .rows()
        .into_iter()
        .zip(free)
        .filter(|(_, free)| **free)
        .map(|(r, _)| r.dot(&r))
        .sum()
}

impl Optimizer for JiggleOptimizer {
    fn run(&mut self, frame: &mut Frame) -> Result<OptReport, String> {
        let mut best = frame.coords().map_err(|e| e.to_string())?;
        let n = best.nrows();
        let free: Vec<bool> = match frame.get("atoms").and_then(|a| a.get_bool("free")) {
            Some(col) if col.len() == n => col.iter().copied().collect(),
            _ => vec![true; n],
        };

        let max_delta = self.max_delta;
        let mut best_score = score(&best, &free);
        let mut accepted = 0usize;

        for _ in 0..self.steps {
            let d = [
                self.rng.random_range(-max_delta..max_delta),
                self.rng.random_range(-max_delta..max_delta),
                self.rng.random_range(-max_delta..max_delta),
            ];
            let mut trial = best.clone();
            for i in 0..n {
                if !free[i] {
                    continue; // never move a frozen environment atom
                }
                for k in 0..3 {
                    trial[[i, k]] += d[k];
                }
            }
            let trial_score = score(&trial, &free);
            if trial_score < best_score {
                best = trial;
                best_score = trial_score;
                accepted += 1;
            }
        }

        frame.set_coords(best.view()).map_err(|e| e.to_string())?;
        Ok(OptReport {
            converged: accepted > 0,
            n_steps: self.steps,
            final_energy: best_score,
            final_fmax: 0.0,
        })
    }
}
```

### Step 2 — bind it to a selection

`GenCanPack::with_optimizer` takes two arguments: an `OptimizeSelect` saying
which targets the optimizer sees and how, then the optimizer itself.
`OptimizeSelect::per_copy(names)` hands over one copy at a time;
`OptimizeSelect::joint(names)` hands over every copy of the named targets as a
single movable group. `.with_environment(rcut)` additionally includes every
atom within `rcut` ångström as frozen context, so a chain folds against its
real neighbours rather than empty space.

```no_run
# use molrs::{Frame, optimize::{OptReport, Optimizer}, types::F};
# struct JiggleOptimizer;
# impl JiggleOptimizer { fn new(_: usize, _: F, _: u64) -> Self { Self } }
# impl Optimizer for JiggleOptimizer {
#     fn run(&mut self, _: &mut Frame) -> Result<OptReport, String> { unimplemented!() }
# }
# let targets: Vec<molpack::Target> = Vec::new();
use molpack::{GenCanPack, OptimizeSelect, PackEngine};

let result = GenCanPack::new()
    .with_tolerance(2.0)
    .with_optimizer(
        OptimizeSelect::per_copy(["chain"]).with_environment(8.0),
        JiggleOptimizer::new(25, 0.5, 42),
    )
    .run(&targets, 200)?;
# Ok::<(), molpack::PackError>(())
```

Optimizer notes:

- **Bindings run in the all-type phase only.** That is the phase whose
  coordinate vector holds every molecule, so per-copy indexing is unambiguous.
- **Non-harm is the packer's job, not yours.** After write-back molpack
  re-evaluates the packing objective and reverts your conformer if it got
  worse. A `run` that makes things worse wastes time; it cannot break a pack.
- **A binding is matched by target name.** `Target::with_name("chain")` is what
  `OptimizeSelect::per_copy(["chain"])` looks up; an unmatched name is skipped
  with a log warning, which is the usual cause of "my optimizer never ran".
- **Copies diverge.** Under `per_copy`, each copy's reference conformer is
  relaxed independently, so copies of one target stop sharing a shape.

## Custom `Stage` — a whole packing algorithm

Everything above plugs into an algorithm that already exists: a restraint
changes *what* the packer tries to satisfy, a handler watches it work, an
in-loop optimizer reshapes molecules while it runs. Replacing the algorithm
itself — how molecules get from nothing to a non-overlapping arrangement — is
what the `Stage` seam is for.

A **stage** is one packing algorithm behind four methods. molpack ships three:
`GenCanStage` (rigid-body descent on the shared objective), `GrowStage`
(configurational-bias chain growth) and `LatticeStage` (a self-avoiding walk on
a diamond lattice). They are peers — no stage reaches into another stage's
driver — and each is judged afterwards by the same objective, so none of them
can grade its own work.

```text
pub trait Stage: Send {
    fn name(&self) -> &'static str;
    fn requires(&self) -> Requires;
    fn guarantees(&self) -> Guarantees;
    fn run(&mut self, state: &mut PackState, targets: &[Target],
           budget: &Budget, handlers: &mut [Box<dyn Handler>])
           -> Result<StageOutcome, PackError>;
}
```

Of the four arguments to `run`, `state` is the run's geometry — the subject of
the next section — while `targets` are the molecule types the caller passed to
the engine's `run`, `budget` is the caller's iteration allowance, and
`handlers` are the observers to notify while you work.

`run` itself returns a `Result`, not a bare outcome: a stage that cannot do
its job fails with a named [`PackError`](crate::PackError), never by
returning `Ok` with `converged: false`. The two report different things —
`converged = false` is the honest "I ran to completion and did not reach my
own criterion", while `Err` is "I could not complete the run at all" — and
collapsing the second into the first would let an outright failure through
disguised as a merely unconverged run. The lifecycle propagates that `Err` out
of [`PackEngine::run`](crate::PackEngine::run) on the same path as its own
validation errors: the failing stage gets no `on_stage_end`, the run gets no
`on_finish`, and no half-built result is assembled. molpack's three built-in
stages always return `Ok`.

### The state a stage is handed

`run` takes a `&mut PackState`, which carries three things:

- a [`PackContext`](crate::PackContext) — the run's mutable geometry: per-atom
  radii, the restraint pool, the cell list, and `xcart`, the **lab-frame**
  coordinates, meaning where each atom actually sits in the packing cell;
- a [`RigidView`](crate::RigidView) — the **rigid placement vector**: for every
  free molecule a centre of mass and three Euler angles (the three-parameter
  description of that molecule's orientation). It is `6 · nmol` numbers, the
  centre-of-mass block first, then the Euler block;
- a [`Placed`](crate::Placed) marker — `Placed::None` (nothing placed yet) or
  `Placed::All` (every free molecule has a placement).

Take the first two apart with `state.rigid_split_mut()`, which hands back
`(&mut PackContext, &mut RigidView)` as two disjoint borrows.

Two coordinate frames meet inside the context. `ctx.coor` holds each copy's
**reference conformer** — its atoms measured from that copy's own centre of
mass with no orientation applied, so it describes the molecule's shape and
nothing else. `ctx.xcart` holds the lab-frame positions. The two are related by
`xcart = com + R(euler) · coor`, where `R(euler)` is the rotation matrix built
from the three Euler angles.

That relation decides where a stage must leave its answer. A stage writes
placements into the `RigidView`; once `run` returns, the lifecycle rebuilds
`xcart` from them with `RigidView::write_xcart`. A stage that instead works
directly in lab-frame coordinates — both growth stages do, because a chain is
grown atom by atom — must capture them back into placements with
`RigidView::capture_from_xcart` before returning. Anything left only in `xcart`
is overwritten.

### What `requires` and `guarantees` declare

Both are one-field records over the `Placed` marker above:
`Requires::new(Placed::All)` says "hand me a state whose molecules are already
placed"; `Guarantees::new(Placed::All)` says "when I return, every free
molecule has a placement". They are declarations, not checks — the lifecycle
advances the state's marker to whatever the stage *guaranteed*, never by
inspecting what it actually did. Declaring them truthfully is what lets a
caller decide which stages may legally follow which.

Both types are `#[non_exhaustive]`, so build them with `::new` rather than a
struct literal; that is what lets a second precondition be added later without
breaking your code.

### Re-entrancy: never consume your own configuration

A stage may be run **more than once** on an evolving state. The rule that falls
out of that is short: the second `run` must have every capability the first
had. The only thing a call may consume is scratch it created itself.

The contract is not decorative. The rigid-body stage used to move its bound
in-loop optimizers out of `self` on the first call; from the second call on it
kept running, just *degraded* — with no error and no name for what it had lost.
Borrow your configuration, do not take it.

### `StageOutcome` carries no verdict

`run` returns `Ok(StageOutcome::new(converged, degraded))` on a successful
run, and those two numbers are all a stage reports: whether it met its own
convergence criterion, and how
many times it had to relax a constructive guarantee (growth's hard-core
softening rungs; always `0` on the rigid-body path, which relaxes nothing). The
run's verdict — the largest inter-molecular overlap `fdist` and the largest
restraint violation `frest` — is read off the context *after* `run` returns, by
the shared objective. That is the one-ruler rule: no algorithm reports its own
score.

### A minimal stage

```rust
use molpack::{
    Budget, F, Guarantees, Handler, PackError, PackState, Placed, Requires, Stage, StageOutcome,
    Target,
};

/// Nudges every molecule by a fixed offset. Not useful — just the smallest
/// thing that is still a stage.
struct ShakeStage {
    /// Displacement added to each centre of mass, in ångström.
    step: [F; 3],
}

impl Stage for ShakeStage {
    fn name(&self) -> &'static str {
        "shake"
    }

    /// It moves placements that already exist; it never creates them.
    fn requires(&self) -> Requires {
        Requires::new(Placed::All)
    }

    /// It leaves every molecule placed, because it only moved them.
    fn guarantees(&self) -> Guarantees {
        Guarantees::new(Placed::All)
    }

    fn run(
        &mut self,
        state: &mut PackState,
        _targets: &[Target],
        _budget: &Budget,
        _handlers: &mut [Box<dyn Handler>],
    ) -> Result<StageOutcome, PackError> {
        // `self.step` is read, never taken — re-entrancy contract.
        let (_ctx, x) = state.rigid_split_mut();
        for i in 0..x.nmol() {
            let com = x.com(i);
            x.set_com(
                i,
                [
                    com[0] + self.step[0],
                    com[1] + self.step[1],
                    com[2] + self.step[2],
                ],
            );
        }
        // Placements only: the lifecycle expands them into `xcart`.
        Ok(StageOutcome::new(true, 0))
    }
}
```

### Wiring it into a run

A stage does not run itself. The lifecycle around it — validation, restraint
broadcast, context construction, handler bracketing, the lab-frame rebuild,
frame assembly — lives in exactly one place: [`Pipeline`](crate::pipeline::Pipeline)'s
[`run`](crate::PackEngine::run), reached through the
[`PackEngine`](crate::PackEngine) trait. An **entry** is a type that plugs into
that lifecycle by implementing [`StageFactory`](crate::StageFactory); its
`stages` method is the one thing an entry must supply: the stage(s) it drives,
built from the run's resolved [`EngineSetup`](crate::pipeline::EngineSetup).
`PackEngine` then adds the shared `with_*` builders plus a one-line `run` that
hands the entry to the lifecycle as its sole stage source
(`Pipeline::single(self).run(targets, max_loops)`) — the same code path a
hand-composed pipeline runs, not a second implementation to keep in step.

```rust
# struct ShakeStage { step: [molpack::F; 3] }
# impl molpack::Stage for ShakeStage {
#     fn name(&self) -> &'static str { "shake" }
#     fn requires(&self) -> molpack::Requires { molpack::Requires::new(molpack::Placed::All) }
#     fn guarantees(&self) -> molpack::Guarantees { molpack::Guarantees::new(molpack::Placed::All) }
#     fn run(&mut self, _state: &mut molpack::PackState, _targets: &[molpack::Target],
#            _budget: &molpack::Budget, _handlers: &mut [Box<dyn molpack::Handler>])
#            -> Result<molpack::StageOutcome, molpack::PackError> {
#         Ok(molpack::StageOutcome::new(true, 0))
#     }
# }
use molpack::pipeline::EngineSetup;
use molpack::{F, Handler, PackEngine, PackError, PackSettings, Pipeline, Stage, StageFactory};

pub struct ShakePack {
    settings: PackSettings,
    handlers: Vec<Box<dyn Handler>>,
    step: [F; 3],
}

impl ShakePack {
    pub fn new(step: [F; 3]) -> Self {
        Self { settings: PackSettings::default(), handlers: Vec::new(), step }
    }
}

impl StageFactory for ShakePack {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }

    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        std::mem::take(&mut self.handlers)
    }

    /// The factory's single contribution: which stage(s) this run drives.
    fn stages(&mut self, _setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        Ok(vec![Box::new(ShakeStage { step: self.step })])
    }
}

impl PackEngine for ShakePack {
    fn settings_mut(&mut self) -> &mut PackSettings {
        &mut self.settings
    }
    fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>> {
        &mut self.handlers
    }

    /// One line: hand the entry to the crate's one lifecycle body as its
    /// sole stage source — every preset writes exactly this.
    fn run(
        self,
        targets: &[molpack::Target],
        max_loops: usize,
    ) -> Result<molpack::State, PackError> {
        Pipeline::single(self).run(targets, max_loops)
    }
}
```

`ShakePack` now has every shared builder (`with_seed`, `with_tolerance`,
`with_handler`, `with_global_restraint`, …) and the terminal
`run(&targets, max_loops)`, exactly like `GenCanPack`, because all of them are
provided methods on `PackEngine`. `EngineSetup` (the resolved targets, cell and
context shape) is the one thing `stages` reads to build its stage(s) from —
there is no separate hook for pre-loading a placement vector before a stage
runs. Instead, each stage's own `run` does that as its own prelude: the
shipped GENCAN stage installs its grid, and — only if it inherited a seed via
`with_restart` — copies the seed's placements in, before deciding whether to
skip `initial()` and continue from what is already there.

### Composing stages

A single-entry `run` is `Pipeline::single(self)` under the hood, so the same
lifecycle also runs a hand-built sequence of stages directly, with no entry
type of your own. [`Pipeline::new`](crate::pipeline::Pipeline::new) builds an
empty pipeline; [`with_stage`](crate::pipeline::Pipeline::with_stage) appends
one factory at a time — a chain grower feeding a rigid-body packer, for
example:

```no_run
use molpack::grow::TorsionPrior;
use molpack::{CbmcGrow, GenCanPack, PackEngine, Pipeline, Target};
# let targets: Vec<Target> = Vec::new();
let result = Pipeline::new()
    .with_stage(CbmcGrow::new(TorsionPrior::Uniform))
    .with_stage(GenCanPack::new())
    .with_seed(42)
    .with_periodic_box([0.0; 3], [30.0; 3], [true; 3])
    .run(&targets, 100)?;
# Ok::<(), molpack::PackError>(())
```

`with_stage` treats what a factory carries two ways, so nothing is silently
dropped. Its **handlers are adopted** — appended to the pipeline's own set in
stage order, so they go on to observe every stage in the run, not just the one
that carried them in. Its **shared settings are refused** the moment any knob
is off its default (`PackError::PresetSettingsInsidePipeline`, naming the
offending knob) — tolerance, precision, seed and the cell are one ruler shared
by the whole run, and two stages each bringing their own would leave the
objective with no single ruler to read. That is why the example above sets
`with_seed` and `with_periodic_box` on the `Pipeline` itself rather than on
`CbmcGrow` or `GenCanPack`: those knobs belong to the run, never to one stage
inside it.

### Custom invariants and combinators

Composing stages linearly, as above, covers a pipeline whose stages each run
once. Two situations need more than that: repeating a body of stages until
some condition is met, and refusing to accept a stage's exit until a property
of the resulting state has been verified. molpack answers both with an
[`Invariant`](crate::Invariant) trait plus two `Stage` combinators,
[`Pipeline::with_repeat`](crate::pipeline::Pipeline::with_repeat) and
[`Pipeline::with_guarded`](crate::pipeline::Pipeline::with_guarded).

An [`Invariant`](crate::Invariant) is a property of a
[`PackState`](crate::PackState) a caller can demand — "restraints are
satisfied within tolerance", say — checked *after* a stage has already run.
It is never a second algorithm competing with the stage itself:

```text
pub trait Invariant: Send {
    fn name(&self) -> &'static str;
    fn layer(&self) -> Layers;
    fn check(&self, state: &PackState) -> Vec<Violation>;
}
```

`check` reads the **shared objective's own verdict** off the state — the same
`fdist` / `frest` numbers [`StageOutcome`](crate::StageOutcome) deliberately
does not carry — and reports every way the property is broken as a
[`Violation`](crate::Violation), or an empty `Vec` when it holds. An invariant
that computed its own second measure of "is this satisfied" would be exactly
the mistake the one-ruler rule exists to catch, so reuse the objective's
numbers rather than re-deriving them; molpack's own
[`RestraintsSatisfied`](crate::RestraintsSatisfied) does this by comparing
`state.ctx().frest` against a tolerance. `layer` names which rung of the
repair-cost ladder ([`Layers`](crate::Layers), documented where it is
consumed) a violation sits on — from a wrong bond graph that nothing
downstream can repair, to bond lengths and angles a downstream force-field
minimization fixes for free — which is what tells a caller whether *rerunning
the same stage* can plausibly help at all.

The smallest possible invariant, useful as a placeholder while wiring up a
guarded stage:

```rust
use molpack::{Invariant, Layers, PackState, Violation};

#[derive(Debug, Clone, Copy)]
pub struct AlwaysSatisfied;

impl Invariant for AlwaysSatisfied {
    fn name(&self) -> &'static str {
        "always-satisfied"
    }

    fn layer(&self) -> Layers {
        Layers::EMPTY
    }

    fn check(&self, _state: &PackState) -> Vec<Violation> {
        vec![]
    }
}
```

Two combinators build on `Invariant`, and are themselves
[`Stage`](crate::Stage)s — the pipeline needs no branch for either, because
chain-checking, handler bracketing and the run's final verdict read them
exactly as they read `GenCanPack`.

- [`with_repeat(body, until)`](crate::pipeline::Pipeline::with_repeat) runs a
  `Vec<Box<dyn StageFactory>>` body repeatedly:
  [`Until::Passes(n)`](crate::Until::Passes) stops after exactly `n` passes
  (`Passes(0)` contributes **no stage at all**, never a silently clamped
  single pass), and [`Until::Converged`](crate::Until::Converged)
  stops the first time a pass ends with the body's last stage reporting its
  own convergence — unbounded by construction, so a body that never converges
  repeats until a handler asks the run to stop. Each pass *continues* from
  where the last one left off, the same `Placed::All` continuation a seeded
  run uses, so `Repeat` around a chain-growth stage feeding a rigid-body one
  is a real "connect, then refine, then connect again" recipe rather than `n`
  independent packs.
- [`with_guarded(stage, invariants, on_violation)`](crate::pipeline::Pipeline::with_guarded)
  runs `stage`, then checks every invariant against the state it left. On the
  first broken one, [`OnViolation`](crate::OnViolation) has exactly two
  answers, never a third: `Fail` returns
  [`PackError::InvariantViolated`](crate::PackError::InvariantViolated)
  immediately; `Rerun { max }` reruns **the same stage** up to `max` more
  times and, if it is still broken once that budget is spent, reports
  `converged: false` (logged, not silently dropped) instead of looping
  forever. `Rerun { max: 0 }` is `Fail`.

That "same stage or named failure" rule is not an oversight: it is the
project's law that the *user* picks the packing method — molpack never
guesses on their behalf. A guard that quietly swapped in a different
algorithm when one failed would be exactly that guess, so `OnViolation` has
no such arm and never will.

```no_run
use molpack::grow::TorsionPrior;
use molpack::{
    CbmcGrow, GenCanPack, Invariant, OnViolation, PackEngine, Pipeline, RestraintsSatisfied,
    StageFactory, Target, Until,
};
# let targets: Vec<Target> = Vec::new();

let body: Vec<Box<dyn StageFactory>> = vec![
    Box::new(CbmcGrow::new(TorsionPrior::Uniform)),
    Box::new(GenCanPack::new()),
];
let invariants: Vec<Box<dyn Invariant>> = vec![Box::new(RestraintsSatisfied::new(1e-3))];

let result = Pipeline::new()
    .with_repeat(body, Until::Converged)
    .with_guarded(GenCanPack::new(), invariants, OnViolation::Rerun { max: 2 })
    .with_seed(42)
    .with_periodic_box([0.0; 3], [30.0; 3], [true; 3])
    .run(&targets, 100)?;
# Ok::<(), molpack::PackError>(())
```

Stage notes:

- **Pick your own name.** `name()` is what every `StepInfo` a handler sees
  carries in `info.stage.name`; keep it short and lowercase — the three
  built-ins report `"gencan"`, `"growth"` and `"lattice"`.
- **`targets` is the one source of chemistry.** They are the same objects the
  caller passed to `run`, so a stage never needs a second copy of the molecule
  data.
- **`budget` is an allowance, not a schedule.** Its two fields, `max_loops` and
  `precision`, say how much work the caller is willing to pay for and how small
  the violation maxima must get; how you spend it is your algorithm's business,
  so your stage should document how it reads `max_loops`.
- **Stages are geometry only.** The seam takes no force field and performs no
  chemistry perception; priors, radii and restraints all arrive as user data on
  the targets.

## Testing discipline

molpack keeps **one** test tier: unit tests in a `#[cfg(test)]` module next to
the code that owns the behaviour. There is no `tests/` directory, no benchmark
suite and no regression harness — an end-to-end packing scenario belongs in
`examples/`, not in the test suite.

| Kind | Location | Convention |
|---|---|---|
| Unit test | `#[cfg(test)] mod tests` in the same file | One `#[test]` fn per behavior |
| Large test body | child module (`src/grow/tests/`, `src/pipeline/tests.rs`) | Still `--lib`, still owned by that module |
| Gradient finite-difference | alongside the unit test | ε=1e-5, tol=1e-3 |

Run all:

```bash
cargo test --lib --all-features
cargo test --doc --all-features
```

Rules:

- Every new `AtomRestraint` gets an FD gradient test.
- Every new `Region` gets a boolean-algebra + signed-distance sign
  test and (for hot-path use) an analytic-gradient FD test.

## Common pitfalls

- **Gradient sign.** Every `AtomRestraint` accumulates `∂penalty/∂x`.
  Optimizer negates for descent. If your molecules fly out of the
  region, the gradient has the wrong sign — penalty should point
  toward the violation boundary.
- **Rotation convention.** Single-atom tests pass with both LEFT and
  RIGHT Euler multiplication; multi-atom tests don't. Always test
  Euler changes with ≥ 2 atoms.
- **`Cell<f64>` is not `Sync`.** Use `AtomicU64` +
  `f64::to_bits` / `from_bits` for interior mutability in
  `Send + Sync` contexts.
- **0-based atom indexing.** `Target::with_atom_restraint` uses
  Rust-native 0-based indices. `&[0, 1]` selects the first two atoms.
  Packmol `.inp` files use 1-based — subtract 1 at the parse boundary.
- **An in-loop optimizer only runs in the all-type phase.** If you expected it
  during per-type pre-compaction, you will see no effect there.
- **`radscale` is phase-dependent.** Don't hard-code atomic radii —
  always go through `sys.radius[i]`. The crate-internal `evaluate_unscaled`
  (`src/context/pack_state.rs`) temporarily swaps `radius` with `radius_ini`
  so the numbers it reports are unscaled.
- **PBC boxes must be valid.** Zero-length axis returns
  `PackError::InvalidPBCBox`.

## Contributing flow

1. Write a failing test.
2. Implement until it passes.
3. Run the full gate:
   ```bash
   cargo test --all-features
   cargo clippy -- -D warnings
   cargo fmt --all --check
   ```
