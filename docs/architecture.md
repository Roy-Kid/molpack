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
├── packer.rs           Molpack builder + pack() driver + phase / iteration loops
├── target.rs           Target — molecule type + per-molecule restraints
├── restraint.rs        Restraint trait + 14 concrete *Restraint structs
├── region.rs           Region trait + And/Or/Not + RegionRestraint
├── relaxer.rs          Relaxer / RelaxerRunner + TorsionMcRelaxer
├── handler.rs          Handler trait + 4 built-in observers
├── objective.rs        compute_f / compute_g / compute_fg + Objective impl
├── context/            PackContext = single owner of mutable packing state
│   ├── pack_context.rs
│   ├── model.rs        immutable topology + inputs
│   ├── state.rs        mutable per-iteration state
│   └── work_buffers.rs scratch arrays (xcart, gxcar, …)
├── constraints/        EvalMode / EvalOutput facade
├── gencan/             bound-constrained quasi-Newton optimizer
│   ├── mod.rs          pgencan / gencan / tn_linesearch
│   ├── cg.rs           conjugate-gradient inner solve
│   └── spg.rs          spectral projected gradient fallback
├── initial.rs          initial random placement + restmol pre-fit
├── movebad.rs          worst-molecule perturbation heuristic
├── euler.rs            Euler angles ↔ rotation matrices
├── cell.rs             cell-list neighbor lookup
├── frame.rs            PackContext ↔ molrs::Frame conversions
├── validation.rs       post-pack correctness check
├── script/             .inp parser + lowering to Targets
├── api/                builder facade re-exports
└── bin/molpack/        CLI front-end (cli feature)
```

### Dependency direction

```text
                       lib.rs
                          │
                       packer.rs  (driver — depends on everything below)
                          │
       ┌──────────┬──────┴──────┬──────────┬──────────┐
       ▼          ▼             ▼          ▼          ▼
    target     initial        gencan     movebad    handler
       │          │             │          │
       │          └────────┐    │          │
       ▼                   ▼    ▼          ▼
 restraint + region    context/PackContext
                            │
                            ▼
                       objective.rs   ← hot path
                            │
                            └── constraints/  (EvalMode facade)
```

`target` / `restraint` / `region` are pure data — no driver imports.
`packer` is the only module that imports everything else. `objective`
is the narrow waist through which all per-atom work flows.

## Data flow

```text
USER INPUTS                 ─→  Target / Molpack builders
  Frame, count, restraints,
  handlers, tolerance, seed
                            ─→  pack() entry
                                a. broadcast global → per-target restraints
                                b. snapshot every Target
                                c. build PackContext
                                     ModelData (immutable topology)
                                     RuntimeState (x, coor, radius)
                                     WorkBuffers (xcart, gxcar, scratch)
                                d. flatten restraints → CSR pool
                                e. initial placement → x[0..6·ntotmol]

PER-ITERATION                ─→  evaluate(x, mode, &mut g)
  GENCAN reads x,                  → expand_molecules: x → xcart
  reads f / g via                  → restraint penalties per atom
  &mut dyn Objective              → cell list + pair penalties
                                  → project gradient back: gxcar → g
                                  returns f_total, fdist, frest

OUTPUT                       ─→  PackResult
                                  positions (Frame),
                                  converged, fdist, frest,
                                  per-phase report
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

Three nested loops drive the packer.

### Outer: `pack()` (one call)

```text
fn pack(targets, max_loops):
    validate inputs (non-empty, valid PBC, atoms > 0)
    broadcast Molpack.global_restraints → each target's molecule_restraints
    split targets into free / fixed
    build PackContext
    run init_passes of restmol():        // geometric pre-fit, no pair kernel
        for each free target type:
            place molecules randomly inside their restraints
            relax restraint penalties only
    handlers.on_start, handlers.on_initialized
    for phase in 0 ..= ntype:
        if phase < ntype:
            comptype[i] := (i == phase)  // PER-TYPE pre-compaction
        else:
            comptype[i] := true          // ALL-TYPES main phase
        report := run_phase(phase, max_loops, …)
        if report.error_phase: break
    handlers.on_finish
    build PackResult { frame, converged, fdist, frest, per-phase reports }
```

Why per-type pre-compaction first: if every type optimizes simultaneously
from a random start, cross-type interference traps the solver in shallow
minima. Compacting one type at a time inside its own restraint volume
gives the all-types phase a much better seed.

### Middle: `run_phase` (one phase)

```text
fn run_phase(phase_id, max_loops):
    handlers.on_phase_start(phase_info)
    radscale := discale            // start with inflated radii (default 1.1)
    relax_runners := build relaxer runners for this phase
    // Quick-exit: if the unscaled objective is already below precision,
    // skip the whole phase.
    if evaluate_unscaled(sys, x).below(precision): return Converged
    for loop_idx in 0 .. max_loops:
        result := run_iteration(loop_idx, radscale, relax_runners)
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
fn run_iteration(loop_idx, radscale, runners):
    // 1. Movebad — relocate the K worst molecules.
    if movebad enabled:
        identify atoms with largest restraint + pair penalty
        perturb their COM/Euler within init_box_half_size
    // 2. Relaxers — update reference geometry per type (count == 1 only).
    for (type, runner) in runners:
        runner.on_iter(ref_coords, f_current, &mut evaluate, rng)
        if accepted: write back new ref_coords
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

`PackContext::evaluate` is invoked O(10³–10⁴) times per `pack()` run.
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
versus the prior monomorphic dispatch — well inside the perf budget.

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
molpack.with_global_restraint(r)
  ≡  for t in targets: t.with_restraint(r.clone())
```

There is no separate global-restraint storage path. The broadcast at
`pack()` entry is the implementation.

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
| How is one restraint's penalty computed for one atom? | `restraint.rs::*::f` / `*::fg` |
| Where does `with_global_restraint` broadcast? | `packer.rs::pack` (top of fn) |
| Where is the per-atom CSR pool built? | `packer.rs::pack` (CSR build loop) |
| How are `x` ↔ Cartesian coords expanded? | `objective.rs::expand_molecules`, `euler.rs::eulerrmat` |
| Where is the pair-overlap kernel? | `objective.rs::accumulate_pair_fg_parallel` |
| What does the initial pre-fit do? | `initial.rs::initial`, `initial.rs::restmol` |
| How is precision-based termination tested? | `gencan/mod.rs::packmolprecision` |
| What does `movebad` do? | `movebad.rs::movebad` |
| How is torsion MC wired in? | `relaxer.rs::TorsionMcRelaxer::on_iter` |
| Where does periodic boundary wrap apply? | `context/pack_context.rs::pbc_distance` |
