# Packer

`Molpack` drives the GENCAN-based three-phase optimizer. All parameters
can be set in the constructor or overridden later via `with_*` builder
methods.

## Constructor

```python
from molpack import Molpack

packer = Molpack(
    tolerance=2.0,         # minimum allowed pairwise distance (Å)
    precision=0.01,        # convergence threshold on fdist and frest
    maxit=20,              # GENCAN inner-loop iteration cap
    nloop0=0,              # phase-0 per-type loops (0 = auto)
    sidemax=1000.0,        # hard bound on global box extent (Å)
    movefrac=0.05,         # fraction of atoms moved per move-bad pass
    movebadrandom=False,   # use random-direction moves
    disable_movebad=False, # skip the move-bad phase entirely
    progress=True,         # attach the built-in progress reporter
)
```

All parameters are keyword-only after `tolerance` and have the defaults
shown above.

## Builder methods

Every builder returns a **new** `Molpack`:

```python
packer = packer.with_tolerance(2.0)
packer = packer.with_precision(0.01)
packer = packer.with_maxit(20)
packer = packer.with_nloop0(0)
packer = packer.with_sidemax(1000.0)
packer = packer.with_movefrac(0.05)
packer = packer.with_movebadrandom(False)
packer = packer.with_disable_movebad(False)
packer = packer.with_progress(True)
```

Use `with_progress(False)` to run silently (useful in tests and scripts).

## Periodic boundaries

Free boundary is the default. To enable PBC:

```python
packer = packer.with_pbc([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
# equivalent shorthand (origin at zero):
packer = packer.with_pbc_box([30.0, 30.0, 30.0])
```

Zero-length axes raise `RuntimeError`. See
[Periodic boundaries](periodic-boundaries.md) for the full semantics.

## Running

```python
result = packer.pack(targets, max_loops=200, seed=None)
```

- `targets`   — list of `Target` objects.
- `max_loops` — per-phase outer-iteration budget.
- `seed`      — optional `int` for reproducible placements; `None`
                uses a non-deterministic seed.

Raises `RuntimeError` on failure (invalid PBC, zero targets, etc.).

## PackResult

```python
result.positions   # (N, 3) float64 ndarray — packed coordinates
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

## Reproducibility

Packing is deterministic for a given `(targets, tolerance, precision, seed)`
tuple. Capture the full `(packer, targets, max_loops, seed)` to reproduce
a result later.
