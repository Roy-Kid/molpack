# gencan-anneal-budget — adaptive GENCAN inner-iteration budget during radius annealing

**Status:** SUPERSEDED — root cause found elsewhere; workaround below was never implemented and is **not** needed.
**Origin:** investigation of `examples/pack_solvprotein` running ~5.5× slower than packmol (2026-06-12).
**Resolution:** 2026-06-12 — see *Resolution* immediately below. Fixed in `src/initial.rs` by porting
Packmol's `avoid_overlap` fixed-atom rejection (the one piece of `initial.f90` the port had dropped).

---

## Resolution (the actual root cause — read this first)

The slowdown was **not** in GENCAN, the objective, the gradient, CG, the radscale schedule, or
floating-point summation order. It was a **missing piece of the initialization**: Packmol's
`avoid_overlap` (default **on**, `initial.f90:386-423`) rejects any random molecule COM that lands
within a ±1-cell stencil of a **fixed** atom, retrying up to 20 times. molpack's `initial.rs`
placement loop only checked the region constraint (`restmol(false)`) and **omitted the fixed-atom
rejection entirely** — the `avoid_overlap` keyword was even parsed (`parser.rs`) but never wired to
the packer.

Consequence on solvprotein: the 4,246-atom protein occupies ~18% of the 50 Å packing sphere, so
~15-20% of the 16,500 waters were seeded **inside the protein**. That:
1. inflated the water-phase **initial** objective ~2× (≈6.7e5 vs Packmol's ≈3.3e5) — a
   *seed-independent* factor (it is a fixed fraction of the volume, not RNG luck), and
2. created deep, stiff overlaps that drove the truncated-Newton CG into immediate
   negative-curvature truncation (1-4 iters instead of ~10), so each outer loop made far less
   progress and the radscale anneal had to grind from 1.1 down to 1.0 to escape.

**The earlier diagnosis below was wrong on two counts:** "radscale=1.1 is near-infeasible at the
protein surface" (Packmol solves it in *one* GENCAN call to fdist=0 — it is perfectly feasible once
the waters start *outside* the protein), and "all of init verified faithful" (the fixed-atom
rejection was the gap). Summation-order (Track 2) was a red herring.

### Fix

`src/initial.rs`: port the `avoid_overlap` branch into the random-placement loop — after drawing a
random COM, scan the ±1-cell stencil of `sys.latomfix`; if any neighbor cell holds a fixed atom,
mark `overlap` and retry (bounded by `MAX_GUESS_TRY = 20`), matching `initial.f90:396-423` exactly.
Gated by a new `Molpack::with_avoid_overlap(bool)` builder (default `true`), wired from the parsed
script in `src/script/build.rs`. When there are no fixed atoms the new code path is byte-identical
to the old one, so the 4 fixed-atom-free examples are provably unchanged.

### Evidence (four-way A/B, seed 1234567, `diag_sp` / `/usr/bin/time`)

| init behavior | molpack | packmol |
|---|---|---|
| `avoid_overlap` **on** (default) | **2.9s** | 5.9s |
| `avoid_overlap` **off** | 33.8s | 78.0s |

Disabling avoidance reproduces the stall on *both* codes (molpack 11×, packmol 13×); the symptom is
identical (water overlap plateaus at ≈2.87 at radscale=1.1, escaping only when radscale anneals to
1.0). molpack-on is now ~2× **faster** than packmol — consistent with molpack already beating
packmol on the other 4 examples; the missing init rejection was the sole cause of the outlier.

### Regression coverage / no-regression check

- `tests/packer.rs::avoid_overlap_reduces_work_around_fixed_solute` — fixed solute + solvent
  converges in strictly fewer outer loops with avoidance on (deterministic, seed 1234567: 5 vs 9).
- Re-timed all 5 shipped examples via `diag_sp`: mixture 0.69s, interface 0.13s, bilayer 10.8s,
  spherical 19.1s, solvprotein **2.9s** — all converge (overlap ≤ precision); the 4 non-fixed cases
  are unchanged within noise (their code path does not change).
- **Separate bug found & fixed while validating:** `tests/examples_batch` was failing solvprotein
  not from packing (the packer reaches `fdist=0`) but from validation: `validation.rs::expand_targets`
  enumerated molecules **free-first** while `PackResult::positions()` now returns **declared order**
  (the in-flight frame-topology work). A `fixed` solute declared first thus had its coordinates
  sliced into the wrong molecules, flagging its own ~1.5 Å bonds as inter-molecule overlaps
  (solvprotein: 2501 false pairs). Fixed `expand_targets` to iterate declared order; `examples_batch`
  now passes all 5. Guard: `validation::tests::fixed_target_declared_first_skips_its_internal_contacts`.
- The adaptive-`maxit` workaround proposed below is therefore **unnecessary** and was not built.

---

## Goal

Eliminate the `pack_solvprotein` performance outlier. molpack is already faster than packmol
on 4 of 5 shipped examples, but `pack_solvprotein` — dense water packed around a large
*fixed* protein — runs **5.5× slower** (31.9s vs 5.8s). Make the GENCAN inner-iteration
budget (`maxit`) **adaptive to the radius-annealing phase**: a small budget while the
effective radius scale `radscale > 1.0` (where molpack's GENCAN stalls on the near-infeasible
inflated problem), and the full packmol budget (`20`) for the final solve at `radscale == 1.0`.
Measured effect: solvprotein **31.9s → 8.9s (1.5×)** with **better** final overlap (0 vs
6.9e-3), neutral-to-faster on every other example, output still converged with overlap ≤
precision and zero restraint violation.

## Non-goals

- **Not** closing the underlying per-iteration parity gap (molpack's truncated-Newton makes
  ~30–40% less progress per iteration than packmol on the stiff inflated problem). That is a
  separate, higher-effort follow-up (see *Risks / open questions → Track 2*); this spec only
  removes the wasteful grinding the gap causes.
- **Not** touching the objective, gradient, CG, line search, SPG, radscale schedule, or
  movebad — all verified faithful (see *Diagnosis*).
- **Not** changing default packing *quality* — the final solve keeps `maxit = 20`.

## Public surface

- **Rust:** add a builder `Molpack::with_anneal_maxit(n: usize) -> Self` (default `5`) and a
  field `anneal_iterations: usize` on `Molpack` (src/packer.rs). In `run_iteration`, select
  the per-call budget:
  ```rust
  let maxit_eff = if *radscale > 1.0 { self.anneal_iterations } else { self.inner_iterations };
  ```
  Build `GencanParams { maxit: maxit_eff, maxfc: maxit_eff * 10, .. }` per outer step instead
  of once before the loop.
- **CLI:** optional `--anneal-maxit <n>` flag mapping to the builder (src/bin/molpack/). Omit
  if we prefer zero new CLI surface; the default already fixes the shipped examples.
- **Python:** mirror as `Molpack.with_anneal_maxit` / a `anneal_maxit=` kwarg in the PyO3
  wrapper (python/src/) if the Rust builder is exposed there; otherwise no change.
- **`.inp`:** none. (packmol's `maxit` keyword still maps to `inner_iterations`/final-solve
  budget, preserving script compatibility.)

## Module placement

Entirely within the packer driver: `src/packer.rs` (`Molpack` struct, `pack`, `run_iteration`).
No change to `src/gencan/` (the optimizer stays a faithful port; only the *caller's* per-call
budget changes). CLI/Python surface is a thin pass-through. This respects the layering rule
that GENCAN is packmol-faithful and tuning lives in the packer.

## Numerical contract

- **Final-solve parity preserved:** at `radscale == 1.0` the call still uses `maxit = 20`,
  `maxfc = 200`, matching packmol's `easygencan`. Only annealing-phase calls (radscale>1) are
  shortened — these never produce the final configuration, only intermediate ones the anneal
  discards anyway.
- **Correctness invariant (must hold on every example):** `converged == true`, final
  `fdist ≤ precision` (0.01), final `frest ≤ precision`. Validated on all 5 examples at
  `anneal_maxit = 5`.
- **Determinism unchanged:** same seed (1234567) → same RNG stream; the change only caps inner
  iterations, no new randomness.
- **Mechanism:** the radscale anneal triggers on `fimp < 2` (improvement between outer loops).
  A smaller annealing-phase `maxit` makes each outer loop's `fimp` smaller sooner, so radscale
  reaches 1.0 in fewer outer loops — cutting the wasted grinding at radscale=1.1.

## Test plan

- **Unit/integration (`cargo test -p molcrafts-molpack --lib --tests`)** — must stay green; add
  a `tests/packer.rs` case asserting `with_anneal_maxit` is honored and a solvprotein-like
  fixed-solute + dense-solvent case converges in fewer outer loops than the flat-20 baseline.
- **Regression (`cargo test --release --test examples_batch -- --ignored`)** — must pass
  unchanged (same final overlap tolerance).
- **Benchmark gate (manual, documented in this spec):** the 5-example table below; require
  solvprotein < 2× packmol and no >10% regression elsewhere.
- **Python (`cd python && maturin develop --release && pytest`)** — if the builder is exposed,
  add a smoke test.

## Doc plan

- Update the perf/architecture note describing the radscale anneal to mention the adaptive
  budget and *why* (dense solvation around a fixed solute).
- If CLI/Python surface is added: documenter sync table entries for `--anneal-maxit` /
  `with_anneal_maxit`.

## Risks / open questions

- **Parity deviation:** packmol uses `maxit = 20` for all phases. This is a deliberate
  molpack-specific tuning. Mitigation: only the *annealing* phase deviates; final solve is
  byte-for-byte packmol behavior. Option to gate behind a `strict_parity` mode if required.
- **Tuning constant `5`:** chosen from a sweep (maxit 5/8/10/15/20 → 9.2/15.4/15.0/26.7/34.2s
  on solvprotein); 5 was best and converged with overlap 0. Could be auto-derived but a fixed
  default is simpler and validated.
- **Track 2 (follow-up, separate spec):** close the per-iteration direction-quality gap so
  molpack solves the radscale=1.1 problem in one call like packmol (→ ~1×). Decisive next
  experiment: shared-config A/B — drive packmol's GENCAN from molpack's post-init coords, then
  diff the analytic gradient element-wise (suspect: cell-list pair-reduction summation order,
  Rust iterator-sum vs Fortran sequential loop, compounding via the trust-region feedback).

---

## Diagnosis (evidence captured from the investigation)

### Benchmarks (idle machine, seed 1234567, `/usr/bin/time -p`)

| example | molecules | packmol | molpack maxit=20 | molpack maxit=5 | quality |
|---|---|---|---|---|---|
| mixture | 1,400 | 1.31s | 0.68s (0.5×) | 0.58s | overlap 0 |
| interface | 1,219 | 0.44s | 0.09s (0.2×) | 0.15s | overlap 0 |
| bilayer | 1,100 | 21.1s | 10.5s (0.5×) | 8.2s | overlap 9e-3 |
| spherical | 18,234 | 32.5s | 17.9s (0.55×) | 17.7s | overlap 0 |
| **solvprotein** | 16,551 | 5.83s | **31.9s (5.5×)** | **8.94s (1.5×)** | overlap 0 |

molpack is faster than packmol on 4/5 (incl. spherical, *larger* than solvprotein). The
outlier's distinguishing feature: dense water packed around a **large fixed solute** (4,246-atom
protein), which makes the radscale=1.1 inflated problem near-infeasible at the protein surface.

### Why solvprotein stalls (per-call GENCAN trace, water phase, radscale=1.1)

molpack runs `maxit=20` and hits the cap (`inform=7`) every outer loop, plateauing:

| outer loop | inform | f_entry→f_exit | radscale |
|---|---|---|---|
| 0 | 7 | 6.7e5→1.3e4 | 1.1 |
| 1–6 | 7 | 1.3e4→…→7.85e3 (grind) | 1.1 |
| 7 | 7 | 5.1e4→3e-3 | **1.0** |
| 8 | 0 | 3e-3→8e-5 (converged) | 1.0 |

~6 outer loops (≈140 TN iterations) are wasted grinding at radscale=1.1 before the anneal
escapes to 1.0. Per-iteration vs packmol on the **same** call:

| iter | packmol f | molpack f | packmol CG | molpack CG |
|---|---|---|---|---|
| 7 | 1.03e4 | 3.36e4 | 3 | 3 |
| 8 | **1.34e3** | 2.52e4 | 2 | 2 |
| 18 | 1.46 (solved) | 1.36e4 | 10 | 2 |

Same CG depth per Newton step, but packmol's steps reduce f far more (iter 7→8: 7.7× vs 1.3×).
packmol solves radscale=1.1 in **one 18-iteration call**; molpack needs ~7× more iterations and
relies on the radscale anneal to escape.

### Verified faithful — NOT the bug (do not re-investigate)

Line-by-line vs `/Users/roykid/work/packmol`: objective `fparc` (incl. `short_radius`/`short_tol`,
off by default both sides), gradient `gparc` (`fscale·4·(datom−tol)`), FD Hessian-vector step
`max(sterel·‖x‖∞, steabs)/‖d‖∞` (sterel=1e-7, steabs=1e-10 both), CG recurrence
(`β=rnorm2/rnorm2prev`, dnorm2/dtr updates, force-descent), CG convergence test, line search
(γ=1e-4, β=0.5, σ1=0.1, σ2=0.9, mininterp=4, maxextrap=100), SPG, outer loop, radscale schedule
(`(fimp<2)||(fdist<prec&&fimp<10) → radscale×0.9`), movebad gating (`radscale==1 && fimp≤10`),
and all constants (maxit=20, maxfc=200, epsgpsn=1e-6, delmin=1e-2, ncomp=50, discale=1.1,
precision=0.01). No preconditioner on either side.

### Ruled out (measured)

- **Higher maxit** → monotonically slower (solvprotein: 50→69s, 100→75s, 200→114s).
- **Parallelism (rayon `--threads`)** → negative scaling (1→4→8 = 21→25→33s); pair-loop too
  fine-grained; parallel reduction also perturbs the trajectory (overlap drifts).
- **delmin / seed / nloop** → already at parity (delmin 1e-2 confirmed = packmol effective via
  easygencan override; seed default 1234567; nloop default 200·ntype — both aligned 2026-06-12).
