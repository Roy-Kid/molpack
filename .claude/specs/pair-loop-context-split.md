---
title: pair-loop-context-split — separate what the pair kernel reads from what it writes
status: draft
created: 2026-07-25
---

# pair-loop-context-split

## Summary

Split `PackContext` so the pair kernel's read set and its write set are distinct
types. Today they are the same 60-field struct, which forces two things this
spec undoes: five near-identical traversal loops, and — on the one attempt to
share them — a measured **+81%** on `compute_f`.

**The split is the deliverable, on its own merit.** A 60-field bag threaded
through the hot path is the architecture defect; the five duplicated traversals
are its most visible symptom, and the +81% from sharing them without splitting
first is what happens when the symptom is treated instead.

Performance is measured at every step and reported, but a regression is a
follow-up optimisation task against the corrected structure — not grounds to
restore the bag. What is never acceptable is landing one silently: the number
goes in the commit message either way.

## Domain basis

`pair_term` reads six things: `xcart`, `atom_props`, `short_radius`,
`short_radius_scale`, `latomnext`, and the `any_short_radius` /
`any_fixed_atoms` summary flags. Its callers write two: `fdist_atom` and
`work.gxcar`. Both live in `PackContext`, so the kernel takes `&PackContext` —
a shared borrow of everything — and a caller wanting `&mut sys.fdist_atom`
alongside it does not typecheck. Each of `fparc` / `gparc` / `fgparc` /
`fparc_stats` / `fgparc_into` therefore inlines the same ten-line chain walk;
they differ only in what they do with a contribution.

### What the first attempt got wrong

Commit 291cea0 introduced a `PairView` holding the six read slices, so the
write targets stayed borrowable, and moved the walk into one generic
`walk_chain` taking a sink closure. It was correct — 256 tests, including
serial/parallel equivalence — and unusable:

| bench | before | after |
|---|---|---|
| `pair_kernel/compute_f` | 3.864 µs | **6.979 µs** (+81%) |
| `pair_kernel/compute_fg` | 5.620 µs | 6.565 µs (+17%) |
| `pack_end_to_end` | 896.3 µs | 1033.6 µs (+15%) |

Reverted in 60eeca7. Two candidate causes, neither confirmed, and the next
attempt has to distinguish them before committing to a design:

1. **Fat value in the inner loop.** `PairView` is ~88 bytes of slices, passed
   into a function called once per pair, where the old code passed one pointer.
   Switching to `&PairView` did not visibly recover it, but that was measured on
   a shared login node with 8–16 µs spreads — i.e. not measured.
2. **Aliasing.** The view holds shared slices of the context while the sink
   closure captures `&mut fdist_atom` / `&mut work.gxcar` of the same context.
   If LLVM cannot prove they are disjoint it must reload the slice pointers
   every iteration; the old code did its reads and writes through one `&mut sys`
   reborrow, which it could reason about.

Cause 2 is the one the design has to answer, and splitting the struct answers it
directly: if the read set and the write set are *different objects*, there is
nothing to disambiguate.

## Design

### Two structs, one owner

```rust
/// Everything the pair kernel reads. Never written during a traversal.
pub struct PairInputs {
    xcart: Vec<[F; 3]>,
    atom_props: Vec<AtomProps>,
    short_radius: Vec<F>,
    short_radius_scale: Vec<F>,
    latomnext: Vec<u32>,
    any_short_radius: bool,
    any_fixed_atoms: bool,
}

/// Everything a traversal accumulates into.
pub struct PairOutputs {
    fdist_atom: Vec<F>,
    gxcar: Vec<[F; 3]>,
}

pub struct PackContext {
    inputs: PairInputs,
    outputs: PairOutputs,
    // …the rest, unchanged
}
```

The kernel takes `&PairInputs`; a traversal takes `(&PairInputs, &mut
PairOutputs)`. Those are provably disjoint at the type level, so no closure
needs to capture a borrow of the same object the reads come from, and the
`walk_chain`-with-sink shape becomes available without the aliasing question.

### Then the traversal

Unchanged from the reverted attempt: one `walk_chain<const GRAD, const
VIOLATION>` plus a `scatter_pair_gradient` helper, with each caller supplying a
sink. Land it as a **separate commit** so a regression can be bisected to the
split or the sharing.

### What stays put

`PackContext` keeps its other ~50 fields. This spec is not a general
decomposition — it splits exactly the boundary the hot loop needs, because that
boundary is load-bearing and the rest is not. A broader reorganisation can
follow if it earns its own measurement.

## Files

- `molpack/src/context/pack_context.rs` — `PairInputs` / `PairOutputs`, field
  moves, accessors for the ~40 call sites that touch `xcart` / `gxcar` today
- `molpack/src/objective.rs` — kernel signatures; then `walk_chain`
- `molpack/src/{initial,movebad,packer}.rs` — call sites that read or write the
  moved fields

## Tasks

1. `PairInputs` / `PairOutputs` with accessors; no behaviour change, no
   traversal change. Measure.
2. `walk_chain` + `scatter_pair_gradient`; collapse the five kernels. Measure.
3. If step 2 still regresses, the cause is (1) not (2): try `#[inline(always)]`
   on the sink boundary, or hand-specialise the three serial kernels while
   keeping the two rayon ones shared.

## Testing

Existing coverage is adequate and was what let the reverted attempt pass:
`src/objective.rs` (`compute_fg_parallel_matches_compute_g_serial_large_system`,
`compute_fg_small_system_parallel_matches_serial`: serial/parallel equivalence of
the pair sums; seed parity is `gencan/gencan_pack.rs::gencan_pack_is_deterministic`),
`src/restraint/geometric/tests/gradient.rs` plus the `*_gradient_matches_finite_difference`
tests in `src/restraint/collective/*.rs` (finite-difference parity per restraint
kind and through the objective), and the five example programs (all five official
Packmol examples pack with a clean `State`). The former integration files
were deleted 2026-09-20. Run the lib tests with `--features cli,ff,rayon`; the
rayon tests are feature-gated and are silently skipped otherwise, which is how
they went stale for four commits.

## Performance protocol (binding)

**Every step is measured on a dedicated node before it is committed.** A shared
login node cannot resolve this work: it reported `compute_f` at +5% when a
dedicated node showed parity, and it could not distinguish 7 µs from 16 µs while
diagnosing the +81%.

Submit `bench_ab.sh` (an exclusive-node A/B of `pair_kernel`, `run_iteration`,
`objective_dispatch`, `restraint_eval`, `pack_end_to_end`, `collective_eval`
between the pre-refactor baseline worktree and HEAD) and read the result before
committing.

The number is recorded, not gated. Architecture comes first: a structural fix
that costs a few percent lands with the cost written down and an optimisation
task behind it. A result like the +81% is different in kind — that is not a cost
to note but evidence the design is wrong, and it sends the design back.

## Out of scope

- Replacing the linked-cell layout with molrs's counting-sorted one. That is the
  larger remaining opportunity — `latomfirst` / `latomnext` chase pointers where
  `LinkCell` walks contiguous slices — but it changes pair emission order, and
  therefore floating-point summation order, so it cannot ride along with a
  change whose acceptance test is "results unchanged".
- A general `PackContext` decomposition.
