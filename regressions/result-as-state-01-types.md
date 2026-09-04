# result-as-state-01-types — public frozen outcome is `State`

Public API pin for the 01-types slice. **hard-coded goldens**; no third-party
runtime (no Packmol binary, no live molrs oracle, no subprocess). Provenance:
in-crate fixture, 2026-09-04; tiny pack matches
`python/tests/test_state.py::_make_tiny_pack` (3 copies × 2 atoms) and
`src/gencan/entry.rs::gencan_entry_is_deterministic` (6 copies × 2 atoms).

Owning tests: `state_natoms_reads_frame` (`src/entry/result.rs`),
`python/tests/test_state.py`, `python/tests/test_entry.py::test_public_surface_is_state_not_pack_result`.

## Public names

```rust
use molpack::{GenCanPack, State, Target};

let state: State = GenCanPack::new().run(&targets, 50)?;
let _ = GenCanPack::new().with_restart(&state);
let _ = Target::fixed_from(&state);
```

```python
import molpack

assert "State" in molpack.__all__
assert "PackResult" not in molpack.__all__
assert hasattr(molpack, "State") is True
assert hasattr(molpack, "PackResult") is False
assert hasattr(molpack.GenCanPack(), "with_restart") is True
assert hasattr(molpack.GenCanPack(), "seeded_from") is False
```

## Hard-coded goldens

- Public type name: `State` (Rust and Python). Continuation:
  `GenCanPack::with_restart` / `GenCanPack().with_restart(state)`.
- Python: `type(run_result).__name__ == "State"`.
- `repr(run_result)` starts with `"State("`.
- `type(run_result.intra).__name__ == "IntraResidual"`.
- Tiny pack 3 copies × 2 atoms → `natoms == 6`.
- Deterministic GENCAN fixture (seed 11, tolerance 2.0, 6 copies × 2 atoms)
  → `natoms() == 12`.
