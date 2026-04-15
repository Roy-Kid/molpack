# Packer

`Packer` drives the GENCAN-based three-phase optimizer. Configure it
via the constructor and the `with_*` builder methods; call `pack()` to
execute.

## Constructor

```python
from molpack import Packer

packer = Packer(tolerance=2.0, precision=0.01)
```

- `tolerance` — minimum allowed pairwise distance (Å). Atoms closer
  than this contribute to `fdist`.
- `precision` — convergence threshold on `fdist` and `frest`.

Both are `float64`. Defaults as shown.

## Builder methods

All builders return a new `Packer`:

```python
packer = (
    Packer(tolerance=2.0, precision=0.01)
    .with_maxit(20)              # GENCAN inner-loop iteration cap
    .with_nloop0(200)            # Phase-0 per-type loops (0 = auto)
    .with_sidemax(1000.0)        # hard bound on global box extent
    .with_movefrac(0.05)         # fraction moved per move-bad pass
    .with_movebadrandom(False)   # use random-direction moves
    .with_disable_movebad(False) # skip the move-bad phase entirely
    .with_progress(True)         # attach the built-in ProgressHandler
)
```

`with_progress(True)` is the default — call `.with_progress(False)` to
run silently.

## Periodic boundaries

Free boundary is the default. To enable PBC:

```python
packer = packer.with_pbc([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
# equivalent shorthand (origin at zero):
packer = packer.with_pbc_box([30.0, 30.0, 30.0])
```

Zero-length axes raise `ValueError`. See
[Periodic boundaries](periodic-boundaries.md) for the full PBC
semantics.

## Running

```python
result = packer.pack(targets, max_loops=200, seed=None)
```

- `targets`   — list of `Target` objects.
- `max_loops` — per-phase outer-iteration budget.
- `seed`      — optional `int` for reproducible placements. `None`
                uses a non-deterministic seed.

Raises `RuntimeError` on packing failures (invalid PBC, too-small box,
…). The raised error message is the Rust `PackError::Display`
rendering.

## PackResult

```python
result.positions   # (N, 3) float64 ndarray
result.elements    # list[str]  — one entry per atom
result.natoms      # int
result.converged   # bool  — True iff both fdist and frest < precision
result.fdist       # float — final distance-violation sum
result.frest       # float — final restraint-violation sum
```

Inspect convergence:

```python
if not result.converged:
    print(f"not converged: fdist={result.fdist:.4f} frest={result.frest:.4f}")
```

## Reproducibility

Packing is deterministic for a given `(targets, tolerance, precision,
seed)` tuple. Every builder method that changes the engine invalidates
prior runs — capture the full `(packer, targets, max_loops, seed)` if
you need to reproduce a result later.
