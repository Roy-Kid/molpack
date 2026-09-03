# Extending the Crate

Tutorials for writing your own `AtomRestraint` / `Region` / `Handler` types,
plus the in-loop optimizer seam. Every extension trait in this crate follows
the same shape (direction-3 rule — see [`concepts`](crate::concepts)):

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
use molpack::{InsideBoxRestraint, Target};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);

let target = Target::from_coords(pos, rad, 100)
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3]))
    .with_restraint(PlaneTether { normal: [0.0, 0.0, 1.0], offset: 20.0, k: 1.0 });
```

Built-in `InsideBoxRestraint` and user `PlaneTether` take the same
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
> e.g. the built-in `ProfileMatch`), whose penalty is a function of the entire
> group and whose gradient couples the copies, so the fixed point is *empirical
> distribution = target*.

## Custom `Region`

Goal: a conical region with apex at origin, axis along +z,
half-angle 30°.

```rust
use molrs::types::F;
# use molpack::Region;

#[derive(Debug, Clone, Copy)]
pub struct ConeRegion {
    pub apex: [F; 3],
    pub axis: [F; 3],
    pub half_angle_cos: F,
}

impl Region for ConeRegion {
    fn contains(&self, x: &[F; 3]) -> bool {
        self.signed_distance(x) <= 0.0
    }
    fn signed_distance(&self, x: &[F; 3]) -> F {
        let dx = x[0] - self.apex[0];
        let dy = x[1] - self.apex[1];
        let dz = x[2] - self.apex[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-12 { return 0.0; }
        let axis_dot =
            (dx * self.axis[0] + dy * self.axis[1] + dz * self.axis[2]) / r;
        self.half_angle_cos - axis_dot
    }
    // Default FD gradient is OK for prototypes. Override analytically
    // for hot-path use — see below.
}
```

### Compose with built-ins

```no_run
# use molrs::types::F;
# use molpack::Region;
# #[derive(Debug, Clone, Copy)]
# pub struct ConeRegion { pub apex: [F; 3], pub axis: [F; 3], pub half_angle_cos: F }
# impl Region for ConeRegion {
#     fn contains(&self, _x: &[F; 3]) -> bool { true }
#     fn signed_distance(&self, _x: &[F; 3]) -> F { 0.0 }
# }
use molpack::{InsideSphereRegion, RegionExt, RegionRestraint, Target};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);

let cone = ConeRegion {
    apex: [0.0; 3],
    axis: [0.0, 0.0, 1.0],
    half_angle_cos: (std::f64::consts::PI / 6.0).cos(),
};
let sphere = InsideSphereRegion::new([0.0; 3], 10.0);
let region = cone.and(sphere);

let target = Target::from_coords(pos, rad, 100)
    .with_restraint(RegionRestraint(region));
```

[`RegionExt::and`](crate::RegionExt::and) / `or` / `not` come from a
blanket impl on every `Region`. The resulting type
`And<ConeRegion, InsideSphereRegion>` is static-dispatch — no heap.

### Analytic gradient override

For hot-path use, override `signed_distance_grad` analytically. The
cone above:

```rust
# use molrs::types::F;
# use molpack::Region;
# #[derive(Debug, Clone, Copy)]
# pub struct ConeRegion { pub apex: [F; 3], pub axis: [F; 3], pub half_angle_cos: F }
# impl Region for ConeRegion {
#     fn contains(&self, _x: &[F; 3]) -> bool { true }
#     fn signed_distance(&self, _x: &[F; 3]) -> F { 0.0 }
fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
    let dx = x[0] - self.apex[0];
    let dy = x[1] - self.apex[1];
    let dz = x[2] - self.apex[2];
    let r2 = dx * dx + dy * dy + dz * dz;
    let r = r2.sqrt();
    if r < 1e-12 { return [0.0; 3]; }
    let axis_dot = (dx * self.axis[0] + dy * self.axis[1] + dz * self.axis[2]) / r;
    let inv_r = 1.0 / r;
    // signed_distance = cos(α) - axis_dot ⇒ grad = -∂axis_dot/∂x
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

## Custom in-loop optimizer (feature `ff`)

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
binder all live behind molpack's `ff` Cargo feature, which pulls in molrs's
force-field module.

![Flexible chains packed inside a spherical cavity](assets/images/paper-confinement-sphere.png)

Confinement like this is a molpack extension workflow, not a Packmol parity
claim: an in-loop optimizer reshapes flexible chains while the packer places
them, so one engine can fit them inside a cavity no rigid pose would clear.
Two runnable programs drive `TorsionMcOptimizer` exactly this way —
`cargo run --release --example pack_adsorption --features ff` and
`cargo run --release --example pack_translocation --features ff`.

### Step 1 — implement `Optimizer`

Goal: a *jiggle* optimizer that proposes random rigid translations of the group
it is handed and keeps the ones that lower a score. A real optimizer scores
with a force field; this one just pulls the group toward the origin, so the
mechanics stay visible.

Two conventions the packer relies on:

- **Coordinates arrive and leave through the `Frame`.** Read them with
  `molrs::ff::potential::extract_coords` — a flat `[x₀, y₀, z₀, x₁, …]` buffer
  in ångström — and write them back with `write_coords`.
- **Frozen atoms are flagged.** When the binding asks for the local
  environment, the `Frame` also carries neighbouring atoms that must not move.
  They are marked by a boolean `atoms.free` column; a missing column means
  every atom is free.

```rust
use molrs::Frame;
use molrs::ff::potential::{extract_coords, write_coords};
use molrs::optimize::{OptReport, Optimizer};
use molrs::types::F;
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
fn score(coords: &[F], free: &[bool]) -> F {
    (0..free.len())
        .filter(|&i| free[i])
        .map(|i| {
            coords[3 * i] * coords[3 * i]
                + coords[3 * i + 1] * coords[3 * i + 1]
                + coords[3 * i + 2] * coords[3 * i + 2]
        })
        .sum()
}

impl Optimizer for JiggleOptimizer {
    fn run(&mut self, frame: &mut Frame) -> Result<OptReport, String> {
        let mut best = extract_coords(frame)?;
        let n = best.len() / 3;
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
                trial[3 * i] += d[0];
                trial[3 * i + 1] += d[1];
                trial[3 * i + 2] += d[2];
            }
            let trial_score = score(&trial, &free);
            if trial_score < best_score {
                best = trial;
                best_score = trial_score;
                accepted += 1;
            }
        }

        write_coords(frame, &best)?;
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

The next snippet is marked `ignore`, so rustdoc shows it without compiling it:
`with_optimizer` and `OptimizeSelect` exist only in an `ff` build, and doctests
run on the default feature set.

```ignore
use molpack::{GenCanPack, OptimizeSelect, PackEngine};

let result = GenCanPack::new()
    .with_tolerance(2.0)
    .with_optimizer(
        OptimizeSelect::per_copy(["chain"]).with_environment(8.0),
        JiggleOptimizer::new(25, 0.5, 42),
    )
    .run(&targets, 200)?;
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

## Testing discipline

| Kind | Location | Convention |
|---|---|---|
| Unit test | `#[cfg(test)] mod tests` in the same file | One `#[test]` fn per behavior |
| Integration test | `tests/<name>.rs` | `use molpack::{…};` only public API |
| Gradient finite-difference | alongside unit test | ε=1e-5, tol=1e-3 |
| Regression vs Packmol | `tests/examples_batch.rs` (`#[ignore]`) | Run with `--ignored --release` |

Run all:

```bash
cargo test --all-features
cargo test --release --test examples_batch -- --ignored
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
  always go through `sys.radius[i]`. `evaluate_unscaled` temporarily
  swaps `radius` with `radius_ini` for user-facing numbers.
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
