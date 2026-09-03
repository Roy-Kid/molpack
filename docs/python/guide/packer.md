# Packer

`GenCanPack` drives the GENCAN-based three-phase optimizer. It is one of
two engine entries — the other, `CbmcGrow`, grows chains instead of
placing rigid bodies; see [Chain growth](growth.md). You pick the
algorithm by picking the entry, and both share the builder names below.
All tuning is through `with_*` builder methods — the `GenCanPack`
constructor takes no arguments.

## Constructor

```python
from molpack import GenCanPack

packer = GenCanPack()
```

All defaults match Packmol's reference behaviour. Override any of
them via the builders below.

## Builder methods

Every builder returns a **new** `GenCanPack`:

```python
packer = (
    GenCanPack()
    .with_tolerance(2.0)            # minimum allowed pairwise distance (Å)
    .with_precision(0.01)           # convergence threshold on fdist and frest
    .with_inner_iterations(20)      # GENCAN inner-loop cap (Packmol `maxit`)
    .with_init_passes(0)            # initial compaction passes (Packmol `nloop0`; 0 = auto)
    .with_init_box_half_size(1000)  # hard bound on init placement (Packmol `sidemax`)
    .with_perturb(0.05, False, True)  # stall heuristic: fraction, random pick, on/off
    .with_seed(42)                  # deterministic RNG
    .with_parallel_eval(False)      # rayon-backed pair-kernel eval (opt-in)
    .with_progress(True)            # LAMMPS-style screen output
)
```

Runs are silent by default; `.with_progress(True)` turns the screen log
on and `.with_progress(False)` turns it back off. Finer control:
`.with_log_level("quiet" | "summary" | "progress" | "verbose")` (an
explicit level wins over `with_progress`) and `.with_log_frequency(n)`
to print every *n*-th loop.

A few more builders cover specific needs:

```python
packer = (
    packer
    .with_periodic_box([0, 0, 0], [30, 30, 30])  # fully-periodic cell (Packmol `pbc`)
    .with_density(0.9)                            # or: size a cubic cell from g/cm³
    .with_avoid_overlap(True)                     # reject init placements onto a fixed molecule
)
```

`with_tolerance`, `with_precision`, `with_seed`, `with_periodic_box`,
`with_density`, `with_parallel_eval`, `with_progress`, `with_handler`,
and `with_global_restraint` are the shared entry builders — they exist
on `CbmcGrow` too. The rest are GENCAN-only.

## Global restraints

Attach a single restraint to every target in a pack:

```python
packer = packer.with_global_restraint(
    InsideBoxRestraint([0, 0, 0], [40, 40, 40])
)
```

Semantically equivalent to calling `.with_restraint(r)` on every
target.

## Handlers

Attach any object implementing some subset of `on_start(ntotat, ntotmol)`,
`on_step(info) -> bool | None`, `on_finish()`:

```python
class MyHandler:
    def on_step(self, info):
        print(f"phase={info.phase} loop={info.loop_idx} fdist={info.fdist:.3f}")
        return None  # or True to request early stop

packer = packer.with_handler(MyHandler())
```

Returning `True` from `on_step` halts the run at the next check. See
the `Handler` Protocol in `molpack`.

## Periodic boundaries

PBC can be declared per-axis on an `InsideBoxRestraint`, or as a
fully-periodic cell directly on the entry via
`.with_periodic_box(min, max)`. See
[Periodic boundaries](periodic-boundaries.md).

## Running

```python
result = packer.run(targets, max_loops=200)
```

- `targets`   — list of `Target` objects (must be non-empty).
- `max_loops` — per-phase outer-iteration budget.

Raises one of the typed `PackError` subclasses on failure
(`NoTargetsError`, `InvalidPBCBoxError`,
`ConflictingPeriodicBoxesError`, …).

`run()` is the entry's only terminal verb and it consumes the entry —
one engine, one run. Calling `run()` twice on the same object raises
`RuntimeError`; build a fresh `GenCanPack` for the next pack.

## PackResult

```python
result.positions   # (N, 3) float64 ndarray — packed coordinates
result.frame       # molrs.Frame — topology-complete packed frame
result.elements    # list[str]  — one entry per atom
result.natoms      # int
result.converged   # bool — True iff both fdist and frest < precision
result.fdist       # float — final distance-violation sum
result.frest       # float — final restraint-violation sum
result.softened    # int — growth-only; always 0 on the GenCanPack path
```

Inspect convergence:

```python
if not result.converged:
    print(f"not converged: fdist={result.fdist:.4f} frest={result.frest:.4f}")
```

`PackResult.frame` is the packed frame. Pass it to a writer of your
choice (e.g. `molrs.io.write_pdb`). molpack does **not** provide
writers.

## Reproducibility

Packing is deterministic for a given `(targets, tolerance, precision,
seed)` tuple. Capture the builder chain and `max_loops` to reproduce
a result later.
