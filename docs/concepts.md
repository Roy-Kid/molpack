# Core Concepts

This chapter defines each abstraction in the crate in one place.
Cross-link to the types for full API details.

## AtomRestraint

An [`AtomRestraint`](crate::AtomRestraint) is a **soft penalty** applied per
atom: `f(x, scale, scale2) -> F` and `fg(x, scale, scale2, g) -> F`.
It contributes to the packing objective and — in all current
implementations — derives from Packmol's `comprest.f90` / `gwalls.f90`.

```text
pub trait AtomRestraint: Send + Sync + std::fmt::Debug {
    fn f (&self, x: &[F; 3], scale: F, scale2: F) -> F;
    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F;
    fn is_parallel_safe(&self) -> bool { true }
    fn name(&self) -> &'static str { std::any::type_name::<Self>() }
}
```

The crate ships 14 concrete `*Restraint` structs (one per Packmol
`kind` 2..=15), each holding its own semantically-named geometric
fields. User types `impl Restraint` sit in the same type slot —
there is no `Builtin*` wrapper in the public API. See
[`extending`](crate::extending) for a tutorial.

### Gradient convention

`fg` accumulates the TRUE gradient (∂penalty/∂x) INTO `g` with `+=`.
Do not overwrite: many restraints may touch the same atom. The
optimizer negates for descent.

### Two-scale contract

Packmol convention (mirrored in the port):

- Linear penalties — kinds 2, 3, 6, 7, 10, 11 (box / cube / plane) —
  consume `scale`.
- Quadratic penalties — kinds 4, 5, 8, 9, 12, 13, 14, 15 (sphere /
  ellipsoid / cylinder / gaussian) — consume `scale2`.

Each `impl Restraint` picks one internally. User-defined restraints
may ignore both knobs and use their own stiffness coefficient as an
instance field.

## CollectiveRestraint

A [`CollectiveRestraint`](crate::restraint::Restraint) is a
**group-level** penalty — unlike [`AtomRestraint`](#atomrestraint), which sees one atom
at a time and contributes an independent external field Σᵢ U(xᵢ), a
collective restraint sees *every* copy of a species at once and returns a
single penalty whose gradient is **coupled across the whole group**:

```text
pub trait Restraint: Send + Sync + std::fmt::Debug {
    fn f (&self, coords: &[[F; 3]], scale: F, scale2: F) -> F;
    fn fg(&self, coords: &[[F; 3]], scale: F, scale2: F, grads: &mut [[F; 3]]) -> F;
    fn is_parallel_safe(&self) -> bool { true }
    fn name(&self) -> &'static str { std::any::type_name::<Self>() }
}
```

`coords` and `grads` have equal length — one entry per atom in the group.
The gradient convention mirrors [`AtomRestraint`](#atomrestraint): `fg` accumulates INTO
`grads[i]` with `+=`.

### Why collective?

A per-atom field built from a target density ρ\* (e.g. Boltzmann inversion
U = −kT·ln ρ\*) is minimised by driving *every* site to the single minimum
of U — the mode of ρ\* — so the sites collapse onto the peak instead of
spreading over the distribution. A collective penalty over the whole
species that matches the empirical distribution to a target profile via
the squared 1-D Wasserstein (sorted-CDF) distance has the correct fixed
point: *empirical distribution = target*.

### Geometry × distribution cross-product

Every member of this family matches a target distribution of a scalar
reaction coordinate ξ defined by a **geometry**, via the Wasserstein
engine. The two axes are orthogonal:

- **Geometry** — maps Cartesian coordinates to ξ and scatters ∂L/∂ξ back:
  `plane` (ξ = signed distance to a plane → slab profile), `point`
  (ξ = distance to a centre → radial profile).
- **Distribution** — the target quantile function q(p) = F⁻¹(p):
  Gaussian, exponential, or tabulated (arbitrary user-supplied profile).

Concrete types are the cross product, named `<Distribution><Geometry>`:

| | Gaussian | Exponential | Tabulated |
|---|---|---|---|
| Plane | `GaussianPlane` | `ExponentialPlane` | `TabulatedPlane` |
| Point | `GaussianPoint` | `ExponentialPoint` | `TabulatedPoint` |

Attach with `Target::with_collective_restraint(r)`. A `TabulatedPlane`
with a histogram from a target simulation is the one-line way to drive a
species toward an experimentally-observed density profile.

## Region

A [`Region`](crate::Region) is a **geometric predicate** with a signed
distance function:

```text
pub trait Region: Send + Sync + std::fmt::Debug {
    fn contains(&self, x: &[F; 3]) -> bool;
    fn signed_distance(&self, x: &[F; 3]) -> F;
    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] { /* default FD */ }
    fn bounding_box(&self) -> Option<Aabb> { None }
}
```

Regions compose via the zero-cost combinators
[`And`](crate::And) / [`Or`](crate::Or) / [`Not`](crate::Not), with
analytic chain-rule gradients (max / min / negate). The
[`RegionExt`](crate::RegionExt) trait gives every `Region` ergonomic
`.and(...)` / `.or(...)` / `.not()` methods.

Any `Region` lifts to a `Restraint` via
[`RegionRestraint<R>`](crate::RegionRestraint):

```text
penalty(x) = scale2 * max(0, signed_distance(x))²
```

Use `Region` when you want compositional geometry (intersection /
union / complement). Use `Restraint` directly when you want a specific
penalty shape (linear vs quadratic, custom stiffness, multi-atom).

## In-loop optimizer (feature `ff`)

Rigid-body packing never changes a molecule's internal shape. An **in-loop
optimizer** does: it rewrites a copy's **reference geometry** — the conformer
all of that copy's Cartesian positions are generated from — between outer
GENCAN calls. Use cases: torsion Monte-Carlo sampling for flexible chains
(propose a rotation about a rotatable bond, accept it with probability
`min(1, exp(-ΔE / T))`), or force-field minimisation of a strained conformer.

The seam is molrs's `Optimizer` trait, one method over a `Frame`:

```text
pub trait Optimizer: Send + Sync {
    fn run(&mut self, frame: &mut Frame) -> Result<OptReport, String>;
}
```

An implementation is attached to the *engine*, not the target, by naming which
targets it applies to:

- `OptimizeSelect::per_copy(names)` — every copy of every named target,
  independently. `OptimizeSelect::joint(names)` — all of them as one movable
  group. `.with_environment(rcut)` adds nearby atoms as frozen context.
- `GenCanPack::with_optimizer(select, optimizer)` binds the pair. Bindings are
  resolved by target name once at `run()` entry and fire in the all-type phase
  only, where the coordinate vector covers every molecule.
- After each call molpack re-evaluates the packing objective and reverts the
  conformer if it got worse, so an optimizer cannot damage a pack.

Because each copy is relaxed on its own, copies of one target start identical
and then diverge — there is no `count == 1` restriction.

Built-in: `TorsionMcOptimizer` (Metropolis torsion sampling with
self-avoidance). All the names in this section — `OptimizeSelect`,
`TorsionMcOptimizer`, `with_optimizer` — require the `ff` Cargo feature, so
they are written in plain code font here rather than linked.

## Handler

A [`Handler`](crate::Handler) is an observer invoked at well-defined
lifecycle points:

```text
pub trait Handler: Send {
    fn on_start        (&mut self, ntotat, ntotmol)       {}
    fn on_initialized  (&mut self, sys: &PackContext)     {}
    fn on_step         (&mut self, info: &StepInfo, sys);   // required
    fn on_phase_start  (&mut self, info: &PhaseInfo)      {}
    fn on_phase_end    (&mut self, info, report: &PhaseReport) {}
    fn on_inner_iter   (&mut self, iter, f, sys)          {}
    fn on_stage_start  (&mut self, info: &StageInfo)      {}
    fn on_stage_end    (&mut self, info: &StageInfo, outcome: &StageOutcome, sys) {}
    fn on_finish       (&mut self, sys: &PackContext)     {}
    fn should_stop     (&self) -> bool                    { false }
}
```

Observer contract: `sys` is always `&PackContext`, never `&mut`.
Handlers cannot modify packer state — bind an in-loop optimizer if you need to.

Built-ins: [`NullHandler`](crate::NullHandler),
[`ProgressHandler`](crate::ProgressHandler),
[`EarlyStopHandler`](crate::EarlyStopHandler),
[`XYZHandler`](crate::XYZHandler).

### Which stage a callback came from

A **stage** is one packing algorithm behind the crate's packing seam — the
[`Stage`](crate::Stage) trait, whose implementors are `GencanStage`,
`GrowStage` and `LatticeStage`; [`extending`](crate::extending) walks through
writing one. Every `StepInfo` names the stage that emitted it in `info.stage`,
a [`StageInfo`](crate::handler::StageInfo) with three fields: `index` (0-based
position of the stage in the run), `total` (how many stages the run has), and
`name` (the stage's own [`Stage::name`](crate::Stage::name), e.g. `"gencan"`).
A run driven by one engine entry has one stage, so it reports `index = 0` and
`total = 1`.

`on_stage_start` and `on_stage_end` bracket a whole stage the way
`on_phase_start` / `on_phase_end` bracket one GENCAN phase. Both are provided
no-ops, and a single-stage run calls neither; they are the seam a caller that
chains stages itself brackets each stage with. `on_stage_end`'s
[`StageOutcome`](crate::StageOutcome) deliberately carries no verdict, only
`converged` and `degraded` (how many times the stage had to relax a
constructive guarantee). The violation maxima are read off `sys`, the
post-stage `PackContext` — the same place `on_finish` reads them — so the
shared objective stays the only ruler.

[`Stage::run`](crate::Stage::run) itself is fallible: it returns
`Result<StageOutcome, PackError>`, not a bare `StageOutcome`. A stage that
cannot do its job fails with a named [`PackError`](crate::PackError), never
by reporting `converged: false` — that flag means only "I ran to completion
and did not reach my own criterion". A stage that returns `Err` gets no
`on_stage_end` at all; the run propagates the error instead, on the same path
as its own validation errors.

## Pipeline

A [`Pipeline`](crate::pipeline::Pipeline) is what actually runs a sequence of
stages. [`Pipeline::new().with_stage(a).with_stage(b)`](crate::pipeline::Pipeline::with_stage)
composes as many stages as a run needs, and
[`Pipeline::single(engine)`](crate::pipeline::Pipeline::single) wraps one
[`PackEngine`](crate::PackEngine) the same way — which is why every preset's
`run` is one line, `Pipeline::single(self).run(targets, max_loops)`. There is
exactly one lifecycle in the crate: a hand-composed pipeline and a preset's own
run are the same code path, never two implementations that could drift apart.

What a pipeline tracks between stages is the [`Placed`](crate::Placed) marker
from the Stage explanation above: it starts at `Placed::None`, and after each
stage returns it advances to whatever that stage's
[`Guarantees`](crate::Guarantees) declared — never by inspecting what the
stage actually did. A stage whose [`Requires`](crate::Requires) needs
`Placed::All` where the marker is still `Placed::None` is a named error,
[`PackError::StageOrder`](crate::PackError::StageOrder), raised before any
handler is notified and before any stage runs.

The run's verdict — `fdist`, `frest`, `converged` — is read off the shared
[`PackContext`](crate::PackContext) after the *last* stage returns, never
assembled from what the individual stages self-reported: the same one-ruler
rule the Stage explanation describes for a single algorithm, applied across a
whole chain. `on_stage_start` / `on_stage_end` bracket each stage in turn, so a
handler watching a two-stage run sees both brackets fire and `info.stage.index`
move from `0` to `1` partway through, while `on_start` / `on_finish` still
bracket only the run as a whole, once.

A stage source handed to `with_stage` may carry handlers and shared settings of
its own. Its handlers are **adopted** — appended to the pipeline's own set, in
stage order, so they go on to watch every later stage too. Its
[`PackSettings`](crate::PackSettings) are **refused** the moment any knob
differs from the default
([`PackError::PresetSettingsInsidePipeline`](crate::PackError::PresetSettingsInsidePipeline),
naming the offending knob): tolerance, precision, seed and the cell are one
ruler for the whole run, and two stages each bringing their own would leave
that ruler ambiguous. Declare shared knobs on the `Pipeline` itself instead.

## Invariants and combinators

Composing stages linearly, as above, covers a pipeline whose stages each run
once. Two situations need more than that: repeating a body of stages until
some condition holds, and refusing to accept a stage's exit until a property
of the resulting state has been verified.

An [`Invariant`](crate::Invariant) is a property of a
[`PackState`](crate::PackState) that is checked *after* a stage has already
run — never a second algorithm competing with the stage itself:

```text
pub trait Invariant: Send {
    fn name(&self) -> &'static str;
    fn layer(&self) -> Layers;
    fn check(&self, state: &PackState) -> Vec<Violation>;
}
```

`check` reads the shared objective's own verdict off the state (`frest`, for
the built-in [`RestraintsSatisfied`](crate::RestraintsSatisfied)) and reports
every broken property as a [`Violation`](crate::Violation) — never a second,
independently computed metric; that would be exactly the mistake the
one-ruler rule exists to catch. `layer` names where a violation sits on the
repair-cost ladder, [`Layers`](crate::Layers): a six-rung bit set ordered from
the almost-unrepairable down to the cheapest to fix.

- **L0 connectivity** — which atoms are bonded to which; nothing downstream
  can repair a wrong bond graph, so it is fixed once, when the template is
  read, and never again.
- **L1 topological state** — knots, entanglement, catenation; undoing one
  needs a chain to pass through itself, which no local move or minimizer can
  do.
- **L2 chain statistics** — end-to-end distance, radius of gyration,
  orientation; fixing these means re-growing the chain, reptation-scale
  motion far beyond a packing run.
- **L3 density** — density and its homogeneity; repairable only by moving
  whole molecules between regions, global and slow but mechanical.
- **L4 local overlaps** — overlaps between neighbouring atoms; the classic
  push-off, removed by a short descent on the shared objective.
- **L5 local geometry** — bond lengths and angles; the cheapest rung, fixed
  for free by the user's own force field in the first steps of minimization.

[`Pipeline::with_repeat(body, until)`](crate::pipeline::Pipeline::with_repeat)
runs a body of stages repeatedly: [`Until::Passes(n)`](crate::Until::Passes)
stops after exactly `n` passes (`Passes(0)` contributes no stage at all,
never a silently clamped single pass), and
[`Until::Converged`](crate::Until::Converged) stops the first time a pass's
last stage reports its own convergence. Each pass
*continues* from where the previous one left off — the same `Placed::All`
continuation a seeded run uses — rather than packing again from nothing.

[`Pipeline::with_guarded(stage, invariants, on_violation)`](crate::pipeline::Pipeline::with_guarded)
runs `stage`, then checks every invariant against the state it left.
[`OnViolation`](crate::OnViolation) answers a broken invariant two ways,
never a third: `Fail` returns a named
[`PackError::InvariantViolated`](crate::PackError::InvariantViolated)
immediately; `Rerun { max }` reruns **the same stage** up to `max` more times
and, once that budget is spent still broken, reports `converged: false`
rather than looping forever. `Guarded` never switches to a different
algorithm to work around a violation — the user picks the packing method,
and molpack does not guess on their behalf.

Both combinators are themselves [`Stage`](crate::Stage) implementations, so
`Pipeline` needs no branch for either: chain-checking, handler bracketing and
the run's final verdict read a `Repeat` or a `Guarded` exactly as they read
`GenCanPack`.

## Objective

The [`Objective`](crate::objective::Objective) trait abstracts over
what GENCAN sees. `PackContext` implements it; synthetic test
objectives (Rosenbrock / Booth / Beale) can implement it to exercise
the optimizer in isolation.

```text
pub trait Objective {
    fn evaluate(&mut self, x: &[F], mode: EvalMode, g: Option<&mut [F]>) -> EvalOutput;
    fn fdist(&self) -> F;
    fn frest(&self) -> F;
    fn ncf(&self) -> u32;
    fn ncg(&self) -> u32;
    fn reset_eval_counters(&mut self);
    fn bounds(&self, l: &mut [F], u: &mut [F]);
}
```

GENCAN (`pgencan`, `gencan`, `tn_ls`, `spg`, `cg`) takes `&mut dyn
Objective` rather than `&mut PackContext` — the optimizer is
decoupled from the packing state.

## Target

A [`Target`](crate::Target) describes one molecule type:

- Input coordinates + centered reference coordinates.
- Van der Waals radii, element symbols, copy count, name.
- Its attached restraints (per-target + per-atom-subset).
- Its attached collective restraints (one penalty over the whole species).
- Optional fixed placement (Euler + translation).
- Optional Euler-angle bounds (`with_rotation_bound(Axis, Angle, Angle)`).
- Optional per-copy mass override (`with_mass`) for density-sized boxes.
- Intramolecular skip table (`with_special_bonds`; default depth-3
  `[0, 0, 0, 1]`). All-atom explicit hydrogen keeps that table and
  shrinks hydrogen via `with_atom_radius`.
- Optionally built from a previous run's output as one fixed obstacle
  ([`Target::fixed_from(&result)`](crate::Target::fixed_from)) — the
  chaining primitive for staged packs.

The packing algorithm is *not* a target property: you pick it by picking
the engine entry ([`GenCanPack`](crate::GenCanPack) or
[`CbmcGrow`](crate::CbmcGrow)), and every target in that call is packed
by it. An unsupported target/entry combination is a named error, never a
silent fall-back.

Targets are snapshotted at `run()` entry — mutating a `Target` after
passing it to an engine has no effect.

## PackEngine and its entries

[`PackEngine`](crate::PackEngine) is the shared lifecycle: one entry type
per algorithm, all with the same builders and the same terminal verb.
[`GenCanPack`](crate::GenCanPack) is rigid-body GENCAN descent;
[`CbmcGrow`](crate::CbmcGrow) is configurational-bias chain growth.

```text
GenCanPack::new()
    .with_log_level(...)
    .with_handler(...)
    .with_global_restraint(...)  // broadcast to every target
    .with_periodic_box(min, max) // or via periodic InsideBoxRestraint
    .run(&[targets], max_loops)  // -> State
```

Every tuning knob (`with_tolerance`, `with_precision`,
`with_inner_iterations`, `with_seed`, `with_avoid_overlap`, …) has a
Packmol-matching default, so `GenCanPack::new().run(&targets, max_loops)`
is a complete call. You only set a knob to *change* its default — e.g.
`with_avoid_overlap(false)` to let solvent seed inside a fixed solute
(on by default), or `with_seed(n)` to pick a different RNG stream (the
default seed is Packmol's `1_234_567`).

Every setter consumes and returns `self`, and so does `run` — an engine
is one-shot by construction, which is what makes it impossible to lose
its handler set on a second call. To stage two algorithms, run the first
entry and feed its output to the second as a fixed matrix:

```text
let grown = CbmcGrow::new(prior).with_density(0.9).run(&[chain], 60)?;
let full  = GenCanPack::new().run(&[Target::fixed_from(&grown), solvent], 200)?;
```

Both entries return the same [`State`](crate::State) —
`frame`, `fdist`, `frest`, `converged`, `degraded`, `intra`.

## PackContext

[`PackContext`](crate::PackContext) is the single owner of mutable
packing state — coordinates, cell lists, restraint pool, rotation
buffers, counters. All optimizer / movebad / handler code paths take
`&mut PackContext` (for writers) or `&PackContext` (for observers).

Structure (`molpack/src/context/`):

- `ModelData` — topology and inputs (immutable after init).
- `RuntimeState` — borrowed telemetry view over `PackContext`.
- `WorkBuffers` — scratch arrays (xcart, gxcar, radiuswork).

Users rarely touch `PackContext` directly — it reaches them through
handler callbacks and the in-loop optimizer bridge. Power users implementing a
custom `Objective` against synthetic test problems will interact with it.

## Scope equivalence law

```text
engine.with_global_restraint(r)
    ≡  for t in targets { t.with_restraint(r.clone()) }
```

There is no separate "global-restraint" storage path in `PackContext`.
The broadcast happens inside `PackEngine::run()`; each target receives an
`Arc::clone` of every global restraint (refcount bump, not a deep
copy).

Per-atom-subset scope is a method-argument pair
`(indices: &[usize], restraint: impl Restraint)`, not a wrapper
struct. There is no `AtomRestraint` public type.

## Restraint versus Constraint

- **Restraint** = soft penalty (violable; pays energy cost).
- **Constraint** = hard constraint (must satisfy; SHAKE / RATTLE /
  LINCS / Lagrange-multiplier mechanisms).

Packmol implements all 15 "constraints" as soft penalties
(`scale * max(0, d)` or `scale2 * max(0, d)²`). Honest naming ⇒
`Restraint`. This crate does not currently define a `Constraint`
trait; adding hard constraints is future work.

## Direction-3 extension pattern

Every extension trait in this crate follows the same shape:

1. Public trait: `pub trait X`.
2. N concrete `pub struct` types that `impl X`, each holding its own
   semantically-named fields.
3. User types `impl X` identically — zero type-level distinction from
   built-ins.

Forbidden in the public API:

- `Builtin*` / `Native*` / `Packmol*` prefixed wrapper types.
- Tagged-union enums that package N built-ins as a single exposed
  type.
- Builder pattern (`X::new().add(...).add(...)`).
- Injection of composition operators into the main trait
  (composition lives on separate traits, e.g. `Region` vs `Restraint`).
- Wrapper types for per-atom-subset scope (that's a method-argument
  pair, not a type).

If you need internal AoS performance structures (e.g. tagged unions
for hot-path match dispatch), they go behind `pub(crate)` and opt-in
via a crate-private hook — invisible to users.
