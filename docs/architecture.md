# Architecture

Developer-oriented view of the crate. Read [`concepts`](crate::concepts)
first for the abstractions this chapter assumes.

This page covers four things, in order:

1. [Module map](#module-map) — where everything lives
2. [Data flow](#data-flow) — how values travel from user input to packed frame
3. [Algorithms](#algorithms) — pseudo-code for the three nested loops
4. [Hot path](#hot-path-objective-evaluation) — what one objective evaluation does
5. [Invariants and conventions](#invariants-and-conventions) — load-bearing rules

## Module map

```text
src/
├── lib.rs              public re-exports + rustdoc chapters
├── pipeline/           the lifecycle body — the only one in the crate
│   ├── mod.rs          Pipeline + the five-part run() (validate → space/state →
│   │                   chain check → run each stage → assemble)
│   ├── engine.rs       StageFactory + PackEngine traits, EngineSetup
│   ├── combinators.rs  Repeat (Until) + Guarded (OnViolation) — combinators
│   │                   are Stage impls, so mod.rs needs no branch for either
│   └── bracket.rs      the handler bracket: adopted-handler tagging, open_bracket, close_bracket
├── entry/              settings + space + result — no lifecycle, never imports pipeline/
│   ├── mod.rs          PackSettings + LogSpec
│   ├── setup.rs        density / pbc / cell resolution + restraint broadcast
│   └── result.rs       PackResult + the verbatim Placements a run hands the next one
├── target.rs           Target — molecule type + per-molecule restraints + fixed_from
├── restraint/          AtomRestraint trait + geometric/ and collective/ impls
├── region.rs           Region trait + And/Or/Not + RegionRestraint
├── handler.rs          Handler trait + LogLevel + 4 built-in observers
├── objective.rs        compute_f / compute_g / compute_fg + Objective impl
├── context/            PackContext = single owner of mutable packing state
│   ├── pack_context.rs
│   ├── pack_state.rs   PackState — context + Placed + RigidView; evaluate_unscaled
│   ├── rigid_view.rs   RigidView — the 6·ntotmol COM + Euler placement vector
│   ├── build.rs        build_context — context + CSR restraint pool from Targets
│   ├── model.rs        immutable topology + inputs
│   ├── state.rs        mutable per-iteration state
│   └── work_buffers.rs scratch arrays (xcart, gxcar, …)
├── constraints/        EvalMode / EvalOutput facade
├── gencan/             rigid-body path: entry + bound-constrained optimizer
│   ├── entry.rs        GenCanPack — the rigid-body engine entry
│   ├── solver.rs       GencanStage — GENCAN on the Stage seam
│   ├── phases.rs       run_phase / run_iteration
│   ├── mod.rs          pgencan / gencan / tn_linesearch
│   ├── cg.rs           conjugate-gradient inner solve
│   └── spg.rs          spectral projected gradient fallback
├── stage.rs            Stage seam — the interface every packing algorithm implements
├── invariant.rs        Layers (L0–L5 repair-cost ladder) + Invariant trait +
│                       Violation + RestraintsSatisfied — consumed by combinators.rs::Guarded
├── grow/               chain-growth path, peer of the GENCAN path
│   ├── entry.rs        CbmcGrow — the chain-growth engine entry (honest verdicts)
│   ├── lattice/        LatticeStage — diamond-lattice SAW for melt density
│   │                   (entry.rs LatticeGrow entry, saw.rs walk,
│   │                   decorate.rs template rebuild, config.rs leaf)
│   ├── config.rs       GrowConfig / GrowError (leaf — no target/entry imports)
│   ├── prior.rs        TorsionPrior / AnglePrior + C∞ calibration
│   ├── internal.rs     template bond graph → internal-coordinate tree
│   ├── field.rs        OverlapField — cell-listed hard-core / soft-shell probe
│   ├── driver.rs       GrowStage round loop (seeding, retraction, softening)
│   └── moves.rs        propose / commit / retract / relax primitives
├── optimizer/          in-loop conformation optimizers (ff feature)
├── initial.rs          initial random placement + restmol pre-fit
├── movebad.rs          worst-molecule perturbation heuristic
├── euler.rs            Euler angles ↔ rotation matrices
├── frame.rs            PackContext ↔ molrs::Frame conversions
├── assemble.rs         packed coords + targets → topology-complete Frame
├── validation.rs       post-pack correctness check
├── script/             .inp parser + lowering to Targets
└── bin/molpack/        CLI front-end (cli feature)
```

### Dependency direction

```text
                       lib.rs
                          │
                     pipeline/  (the lifecycle — depends on everything below)
                          │
    ┌────────┬────────┬──┴──────┬──────────┬──────────┐
    ▼        ▼        ▼         ▼          ▼          ▼
  entry/   target   initial    gencan     movebad    handler
    │        │                 grow/
    │        └────────┐         │
    ▼                 ▼         ▼
    └───────────► context/PackContext
                            │
                            ▼
                       objective.rs   ← hot path
                            │
                            └── constraints/  (EvalMode facade)
```

`target` / `restraint` / `region` are pure data — no driver imports.
`pipeline/` is the only module that imports everything else; `entry/`
shrank to settings + space + result and imports nothing from `pipeline/` —
the arrow points one way, `pipeline/` reads `entry/`, never the reverse.
`objective` is the narrow waist through which all per-atom work flows.

The chain-growth path enters at the same level as `gencan`: a preset's
`stages()` ([`StageFactory`](crate::StageFactory)) builds one `Stage` from
the `stage` seam, [`Pipeline::run`](crate::PackEngine::run) drives it, and
`grow/` (the `GrowStage`) consumes the same `PackState` and is judged by the
same objective — it never calls the GENCAN internals.

A **stage** is one interchangeable packing algorithm behind four methods:
`name`, `requires`, `guarantees`, and `run`. The middle two declare the
shape of the state the stage needs on the way in and promises on the way
out, as a `Placed` marker — `Placed::None` (nothing placed yet) or
`Placed::All` (every free molecule has a placement). `run` receives a
`PackState`: the run's `PackContext`, plus that marker, plus the rigid
placement vector `RigidView` (three centre-of-mass and three Euler values
per free molecule). `run` returns `Result<StageOutcome, PackError>`: a stage
that cannot do its job fails with a named error, never by disguising failure
as `converged = false`. On the `Ok` side, `StageOutcome` is only what
the stage alone knows — whether it met its own convergence criterion, and
how many times it had to relax a constructive guarantee. The run's
`fdist` / `frest` verdict is read off the state afterwards, so no
algorithm grades its own paper.

A [`Pipeline`](crate::pipeline::Pipeline) composes several stages behind one
call — `Pipeline::new().with_stage(CbmcGrow::new(prior)).with_stage(GenCanPack::new()).run(..)`
runs chain growth, then rigid-body push-off, in one lifecycle, continuing
from the first stage's placements rather than re-placing from scratch. Still
out of scope: parallel or branching stage graphs (v1 is a linear sequence
only). For two genuinely independent packs, run each to completion and hand
the earlier result to the next as a **fixed** obstacle via
`Target::fixed_from(&result)`.

## Data flow

```text
USER INPUTS                 ─→  Target / PackEngine builders
  Frame, count, restraints,       (GenCanPack | CbmcGrow | Pipeline)
  handlers, tolerance, seed
                            ─→  PackEngine::run()      (one line per preset)
                            ─→  Pipeline::run()          the lifecycle body
                                a. broadcast global → per-target restraints
                                b. snapshot every Target
                                c. build PackContext, wrap into PackState
                                     ModelData (immutable topology)
                                     RuntimeState (x, coor, radius)
                                     WorkBuffers (xcart, gxcar, scratch)
                                d. flatten restraints → CSR pool
                                e. per stage: Stage::run(state, targets, …)
                                     GENCAN: install grid, then initial
                                     placement — or continue from a seed

PER-ITERATION                ─→  evaluate(x, mode, &mut g)
  (inside a stage's own loop —      → expand_molecules: x → xcart
   GencanStage for the rigid path)  → restraint penalties per atom
  reads f / g via                   → cell list + pair penalties
  &mut dyn Objective                → project gradient back: gxcar → g
                                     returns f_total, fdist, frest

OUTPUT                       ─→  PackResult
                                  frame, plus converged, fdist,
                                  frest, softened
```

Three rules govern this flow:

- **`PackContext` owns mutable state.** GENCAN, movebad, handlers, and the
  phase driver all take `&mut PackContext` (writers) or `&PackContext`
  (observers). No other module owns mutable state across iterations.
- **`Arc<dyn Restraint>` for polymorphic storage.** Cheap clone (refcount
  bump) into the per-atom CSR pool. The hot path does one virtual call
  per restraint per atom.
- **GENCAN is decoupled.** `gencan/pgencan` takes `&mut dyn Objective`,
  not `&mut PackContext`. Synthetic objectives (Rosenbrock, Booth, Beale)
  exercise the optimizer in isolation.

### Coordinate layout

The optimizer variable vector `x` packs centers of mass and Euler angles:

```text
x = [com₀(3), com₁(3), …, comₙ(3),  eul₀(3), eul₁(3), …, eulₙ(3)]
length = 6 · ntotmol
```

Cartesian atom positions `xcart: Vec<[F; 3]>` of length `ntotat` are
expanded each evaluation:

```text
xcart[icart_for(i, m, a)] = com_m + R(eul_m) · ref_coords[i, a]
```

where `i` is molecule type, `m` is copy index, `a` is atom index.

## Algorithms

One lifecycle loop over stages, then three nested loops inside the
rigid-body stage.

### The lifecycle: `Pipeline::run()` (one call)

```text
fn run(targets, max_loops):
    validate inputs (non-empty, valid PBC, atoms > 0)
    broadcast settings.global_restraints → each target's molecule_restraints
    resolve packing space; build PackContext, wrap into PackState
    check the stage chain (empty list / bad order / a preset's non-default settings)
    handlers.on_start
    for stage in stages:                  // one stage for a preset's own run
        state.invalidate_geometry_cache()
        handlers.on_stage_start
        outcome := stage.run(state, targets, budget, handlers)?  // named PackError
                                                                  // skips on_stage_end + on_finish
        state.set_placed(stage.guarantees().placed)
        handlers.on_stage_end
        if handlers.should_stop(): break
    rebuild xcart from the final rigid view; handlers.on_finish
    assemble Frame into PackResult (+ converged / fdist / frest / softened)
```

Every preset's `PackEngine::run` is one line —
`Pipeline::single(self).run(targets, max_loops)` — so `GenCanPack::run()`,
`CbmcGrow::run()` and `LatticeGrow::run()` all resolve to the loop above. It
lives once, in `src/pipeline/mod.rs`, never duplicated per entry.

### Outer: `GencanStage::run()` (one stage)

```text
fn run(state, targets, budget, handlers):
    if state already placed, or this stage carries a seed:
        install the box + cell grid
    if this stage carries a seed:
        inject the seed's placements verbatim; state.placed := All
    push_off := state.placed == All
    if push_off:
        write_xcart(state)                        // continue from what is there
    else:
        run init_passes of restmol():              // geometric pre-fit, no pair kernel
            for each free target type:
                place molecules randomly inside their restraints
                relax restraint penalties only
    handlers.on_initialized
    for phase in 0 ..= ntype:
        if phase < ntype:
            comptype[i] := (i == phase)  // PER-TYPE pre-compaction
        else:
            comptype[i] := true          // ALL-TYPES main phase
        report := run_phase(phase, max_loops, …)
        if report.error_phase: break
    return Ok(StageOutcome::new(converged, 0))
```

The preamble — box/grid install, seed injection, the `initial()`-vs-push-off
choice — and the phase loop both live in `GencanStage::run`
(`src/gencan/solver.rs`); nothing above the stage boundary decides when
`initial()` (and the `movebad` heuristic it configures) runs. `CbmcGrow`'s
`GrowStage` (`src/grow/driver.rs`) is a peer stage under the same lifecycle,
with the phase loop above replaced by its own growth round loop.

Why per-type pre-compaction first: if every type optimizes simultaneously
from a random start, cross-type interference traps the solver in shallow
minima. Compacting one type at a time inside its own restraint volume
gives the all-types phase a much better seed.

### Middle: `run_phase` (one phase)

```text
fn run_phase(phase_id, max_loops):
    handlers.on_phase_start(phase_info)
    radscale := discale            // start with inflated radii (default 1.1)
    // Quick-exit: if the unscaled objective is already below precision,
    // skip the whole phase.
    if evaluate_unscaled(sys, x).below(precision): return Converged
    for loop_idx in 0 .. max_loops:
        result := run_iteration(loop_idx, radscale, optimizer_bindings)
        radscale := decay(radscale)            // → 1.0 over the phase
        handlers.on_step(step_info, sys)
        if result.converged: return Converged
        if handlers.should_stop(): return EarlyStop
    return MaxLoops
```

`radscale` starts at `discale` (1.1) and decays toward 1.0 over the
phase. This soft-starts the pair penalty: the optimizer first sees
slightly oversized atoms (easier to push apart) and tightens to true
tolerance as the phase progresses.

### Inner: `run_iteration` (one outer step)

```text
fn run_iteration(loop_idx, radscale, optimizer_bindings):
    // 1. Movebad — relocate the K worst molecules.
    if movebad enabled:
        identify atoms with largest restraint + pair penalty
        perturb their COM/Euler within init_box_half_size
    // 2. In-loop optimizers (feature `ff`) — all-type phase only, so that
    //    COM/Euler indexing covers every molecule.
    for binding in optimizer_bindings:
        assemble a Frame per selection (moving copies + frozen neighbours)
        binding.optimizer.run(&mut frame)
        map the displacement back into each copy's reference conformer
        revert the copy if the packing objective got worse
    // 3. GENCAN — bound-constrained quasi-Newton solve.
    pgencan(x, &mut sys, params, precision)
        // Internally: tn_linesearch → CG inner solve → SPG fallback,
        // each step calls sys.evaluate(x, mode, g).
    // 4. Convergence check on the unscaled objective.
    f_unscaled := evaluate_unscaled(sys, x)
    fimp := percentage improvement vs previous loop
    converged := fdist < precision AND frest < precision
    return { converged, fimp, fdist, frest }
```

GENCAN itself runs three nested solvers:

```text
pgencan: project x onto bounds, then call gencan
gencan:  truncated-Newton outer; calls tn_linesearch
tn_ls:   conjugate-gradient line search; SPG fallback if CG stalls
```

Each leaf step calls `sys.evaluate(x, mode, &mut g)` — the hot path.

## Hot path: objective evaluation

`PackContext::evaluate` is invoked O(10³–10⁴) times per `run()`.
Performance lives here.

```text
evaluate(x, mode, g) dispatches by mode:
    FOnly        → compute_f
    GradientOnly → compute_g
    FAndGradient → compute_fg
    RestMol      → compute_fg (init phase, pair kernel skipped)
```

`compute_fg` is the canonical path — it does five steps:

```text
1. expand_molecules(x):
       for each molecule type t, copy m, atom a:
         xcart[icart] := com_t,m + R(eul_t,m) · ref_coords[t, a]

2. accumulate_constraint_value_and_gradient (per atom icart):
       range := iratom_offsets[icart] .. iratom_offsets[icart + 1]
       for &irest in iratom_data[range]:
           f += sys.restraints[irest].fg(xcart[icart], scale, scale2,
                                          &mut grad_xcart[icart])
       // Linear penalties consume `scale`; quadratic consume `scale2`.

3. insert_atom_in_cell (per atom):
       linked-list bucket atoms into cells
       cell side ≈ 2 × max_radius × radscale

4. accumulate_pair_fg (or _parallel under rayon):
       for each non-empty cell c:
         for each neighbor cell c′ in 13-cell stencil:
           for each (i ∈ c, j ∈ c′):
             d  := pbc_distance(xi, xj)
             σ  := (rᵢ + rⱼ) · radscale
             if d < σ:
                 penalty := (σ − d)²
                 grad_xcart[i] += d penalty / d xi
                 grad_xcart[j] += d penalty / d xj

5. project_cartesian_gradient:
       for each molecule m, atom a:
         g_com[m]   += grad_xcart[icart]
         g_euler[m] += Jᵀ(eul_m, ref_a) · grad_xcart[icart]
       // J = ∂xcart/∂eul, derived once per molecule from R(eul).
```

Cost breakdown: steps 1–3 are O(N_atoms); step 4 is
O(N_atoms × neighbor_avg) ≈ O(N_atoms × 32) and dominates wall time on
realistic workloads. Step 4 is the rayon parallelization point
(`accumulate_pair_fg_parallel`), reducing into per-atom gradient slots
via `AtomicU64` (since `Cell<f64>` is not `Sync`).

The `Arc<dyn Restraint>` virtual call in step 2 measured at +0.22% e2e
versus the prior monomorphic dispatch — a negligible cost for the
flexibility of user-defined restraints.

## Invariants and conventions

**Gradient accumulation.** `Restraint::fg` accumulates the true
gradient (∂penalty/∂x) into `g` with `+=`. Optimizer negates for descent.
Multiple restraints may touch one atom, so never overwrite.

**Two-scale contract.** Linear penalties (Packmol kinds 2/3/6/7/10/11)
consume `scale`; quadratic penalties (kinds 4/5/8/9/12/13/14/15) consume
`scale2`. Each `impl Restraint` picks one internally.

**Rotation convention.** `R_new = δR · R_old` (LEFT multiplication).
Single-atom tests cannot detect LEFT/RIGHT bugs — always test with
≥ 2 atoms.

**Coordinate layout.** GENCAN's `x` is `[com₀..n, eul₀..n]` of length
`6·ntotmol`. Cartesian atom positions `xcart` are `Vec<[F; 3]>` of length
`ntotat`.

**Thread safety.** All trait objects are `Send + Sync`. Interior
mutability inside parallel reductions uses `AtomicU64` with
`f64::to_bits` / `f64::from_bits` — `Cell<f64>` is not `Sync`.

**Scope equivalence.**

```text
engine.with_global_restraint(r)
  ≡  for t in targets: t.with_restraint(r.clone())
```

There is no separate global-restraint storage path. The broadcast at
`Pipeline::run()` entry is the implementation.

**Restraint vs Constraint.** Packmol implements all 15 "constraints" as
soft penalties. Naming reflects mechanism, not user intent → `Restraint`.

**Direction-3 extension pattern.** Every extension trait follows the
same shape: public trait, N concrete `pub struct` impls, user types
`impl Trait` identically. No `Builtin*` / `Native*` wrappers in the
public API.

**`init1` short-circuit.** Set during the initial geometric pre-fit.
Skips the pair kernel — the restraint-only objective is enough to get
atoms into their regions before pair conflicts matter.

## Cheatsheet

| Question | Where to look |
|---|---|
| How is one restraint's penalty computed for one atom? | `restraint/*::f` / `*::fg` |
| Where does `with_global_restraint` broadcast? | `entry/setup.rs::broadcast_global_restraints` |
| Where is the per-atom CSR pool built? | `context/build.rs::build_context` (CSR build loop) |
| How are `x` ↔ Cartesian coords expanded? | `objective.rs::expand_molecules`, `euler.rs::eulerrmat` |
| Where is the pair-overlap kernel? | `objective.rs::accumulate_pair_fg_parallel` |
| What does the initial pre-fit do? | `initial.rs::initial`, `initial.rs::restmol` |
| How is precision-based termination tested? | `gencan/mod.rs::packmolprecision` |
| What does `movebad` do? | `movebad.rs::movebad` |
| How is torsion MC wired in? | `optimizer/torsion_mc.rs::TorsionMcOptimizer::run`, called from `optimizer/mod.rs::run_optimizer_bindings` (feature `ff`) |
| Where does periodic boundary wrap apply? | `context/pack_context.rs::pbc_distance` |
