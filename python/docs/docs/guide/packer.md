# Packer

`Molpack` drives the GENCAN-based three-phase optimizer. All tuning
is through `with_*` builder methods — the constructor takes no
arguments.

## Constructor

```python
from molpack import Molpack

packer = Molpack()
```

All defaults match Packmol's reference behaviour. Override any of
them via the builders below.

## Builder methods

Every builder returns a **new** `Molpack`:

```python
packer = (
    Molpack()
    .with_tolerance(2.0)            # minimum allowed pairwise distance (Å)
    .with_precision(0.01)           # convergence threshold on fdist and frest
    .with_inner_iterations(20)      # GENCAN inner-loop cap (Packmol `maxit`)
    .with_init_passes(0)            # initial compaction passes (Packmol `nloop0`; 0 = auto)
    .with_init_box_half_size(1000)  # hard bound on init placement (Packmol `sidemax`)
    .with_perturb_fraction(0.05)    # fraction of atoms re-sampled per stall
    .with_random_perturb(False)     # randomize perturbation target selection
    .with_perturb(True)             # enable the stall-perturbation heuristic
    .with_seed(42)                  # deterministic RNG
    .with_parallel_eval(False)      # rayon-backed pair-kernel eval (opt-in)
    .with_progress(True)            # attach the built-in progress reporter
)
```

Use `.with_progress(False)` to run silently (useful in tests and scripts).

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

PBC is configured on the `InsideBoxRestraint` itself — not the packer.
See [Periodic boundaries](periodic-boundaries.md).

## Running

```python
result = packer.pack(targets, max_loops=200)
```

- `targets`   — list of `Target` objects (must be non-empty).
- `max_loops` — per-phase outer-iteration budget.

Raises one of the typed `PackError` subclasses on failure
(`NoTargetsError`, `InvalidPBCBoxError`,
`ConflictingPeriodicBoxesError`, …).

## PackResult

```python
result.positions   # (N, 3) float64 ndarray — packed coordinates
result.frame       # dict compatible with molrs.Frame
result.elements    # list[str]  — one entry per atom
result.natoms      # int
result.converged   # bool — True iff both fdist and frest < precision
result.fdist       # float — final distance-violation sum
result.frest       # float — final restraint-violation sum
```

Inspect convergence:

```python
if not result.converged:
    print(f"not converged: fdist={result.fdist:.4f} frest={result.frest:.4f}")
```

`result.frame` is the primary Frame-oriented output — pass it to a
writer of your choice (e.g. `molrs.write_pdb`). molpack does **not**
provide writers.

## Reproducibility

Packing is deterministic for a given `(targets, tolerance, precision,
seed)` tuple. Capture the builder chain and `max_loops` to reproduce
a result later.
