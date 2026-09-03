# molpack conventions

Facts that are not laws (those live in `law.md`). Moved here from CLAUDE.md
on 2026-09-02 so the router stays thin; content is unchanged except for the
harness-layout rows.

## Public surfaces

- `molcrafts-molpack` Rust library (`[lib] name = "molpack"`, edition 2024, Rust 1.91+)
- `molpack` binary (`cli` feature)
- `molcrafts-molpack` Python wheel under `python/`, built with maturin **without** the `io` feature

## Cargo features

| Feature | Pulls in |
|---|---|
| `default` | nothing |
| `io` | `molrs/io` (PDB / XYZ / SDF / LAMMPS readers) |
| `cli` | `clap` + `io` (binary + integration tests) |
| `rayon` | `rayon` + `molrs/rayon` (parallel evaluation) |
| `ff` | `molrs/ff` (MMFF + L-BFGS) — enables the in-loop `Optimizer` bindings |

The Python wheel is built **without** `io` — the wheel relies on the user's
`molrs` Python package for frame loading, then builds targets from the loaded
frame (the PyO3 `target_from_frame` helper) and lowers scripts with
`Script::lower` + `StructurePlan::apply`.

## Coding style

- Immutability — return new values, never mutate in place; builders take `self` and return `Self`
- Files: 200–400 lines typical, 800 max — split when a module grows beyond one concern
- New public types implement `Debug` and (where appropriate) `Clone`
- `cargo fmt` and `cargo clippy -- -D warnings` are mandatory before commit
- Leaf config files (`grow/config.rs` and siblings) import only sibling leaves and molrs types — never `target` / `entry` / `context`

## Tests and gates

- TDD: RED first, then GREEN, then refactor. 80% coverage minimum.
- Unit tests live in-module (`#[cfg(test)]`); integration tests in `tests/`, one file per subsystem (`grow.rs`, `pipeline.rs`, `examples_batch.rs`, …); Python tests in `python/tests/`.
- `stage-pipeline` chain, one file per subsystem, each owning only its own type's contract:
  `tests/topology.rs` (the bond-graph leaf `src/topology.rs`); `tests/context_rigid_view.rs`
  (the rigid placement vector `src/context/rigid_view.rs`); `src/context/pack_state/tests.rs`
  (`PackState` + `evaluate_unscaled`, in-crate rather than in `tests/` because both were
  `pub(crate)` when written and so invisible to an integration test in a separate crate — mounted
  from `src/context/pack_state.rs` via `#[cfg(test)] mod tests;` and still collected by the
  ordinary `--lib` gate); `tests/stage.rs` (the `Stage` seam's own contract — object safety,
  `requires`/`guarantees`, re-entrancy — on fake stages only, never a production stage; what each
  concrete stage declares stays with its owner, `tests/gencan.rs` / `tests/grow.rs`);
  `tests/pipeline.rs` (the multi-stage lifecycle body `src/pipeline/`); `tests/invariant.rs`
  (`Layers` / `Invariant` / `RestraintsSatisfied`, `src/invariant.rs`).
- `cargo test -p molcrafts-molpack --lib --tests` — fast tier, must always be green (`mol_project.build.test`). A single test: append `-- <name filter>` (`build.test_single`).
- `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored` — Packmol regression (five official examples, fixed seed). Requires test data: `bash ../molrs/scripts/fetch-test-data.sh` (one time).
- `uv run --directory python --group dev tox -e py` — Python wheel, isolated and non-editable (`maturin develop` is not the gate).
- `cargo bench --benches` — criterion regression benches (pair kernel, objective dispatch, one `run_iteration` step, restraint eval, tiny end-to-end pack). Each synthesizes its own geometry, so **no** `io` feature is needed; one bench: `cargo bench --bench pack_end_to_end`. Perf history is tracked on canonical pushes by `.github/workflows/bench.yml`.
- CI parity (pre-push hooks in `.pre-commit-config.yaml`, mirrored by `mol_project.ci.local`): Rust tests with `io`, the `cli` test, `--no-default-features` / `rayon` / `ff` checks, bench build, then tox.

## Repo layout

| Path | Purpose |
|---|---|
| `src/` | library + CLI binary (`src/bin/molpack/`) |
| `python/` | PyO3 wheel (`python/src/`) + package (`python/python/molpack/`) + tests |
| `tests/` | Rust integration tests (incl. `examples_batch.rs` regression suite) |
| `benches/` | small criterion regression benches (self-synthesized geometry, no `io`) |
| `examples/` | runnable example programs (need `--features io`) |
| `docs/` | public docs site (Rust guide + `docs/python/` binding docs) — published via Zensical (`zensical.toml`) with the shared `molcrafts` theme |
| `.claude/specs/` | active feature specs, indexed in `INDEX.md`; deleted on close |
| `.claude/notes/` | passive knowledge: `law.md`, `conventions.md`, `architecture.md`, `notes.md` |

Skills and agents come from the `mol` plugin (`molcrafts-harness`); the repo
carries no project-local `.claude/skills/` or `.claude/agents/` (the former
`mpk-*` set was removed on 2026-09-02).

## Sibling layout assumed

```
workspace/
├── molrs/      ← git clone https://github.com/MolCrafts/molrs
└── molpack/    ← this repo
```

The root `Cargo.toml` uses a single path dep on `../molrs/molrs` (the unified
`molcrafts-molrs` crate); `core` is always-on, while `io` and `ff` are pulled in
through the matching molpack features above.

**Build cache:** the committed `.cargo/config.toml` routes every build (root
workspace and `python/`) into `../molrs/target`, shared with the sibling molrs
checkout — molrs compiles once per (rustc, features, profile) across both
repos. `rust-toolchain.toml` matches molrs's so the cache fingerprints one
rustc. CI caches that dir and runs sccache.

**molrs ABI line** (the rule itself is law P10): molpack exchanges `molrs_ffi`
handle capsules with the installed `molcrafts-molrs` wheel; both must embed the
same molrs **major.minor** (minor line = ABI version — see molrs
`docs/interop.md`). Gates: `molpack/version.py` (wheel metadata, at
`import molpack`), `interop::check_abi` (`molrs._ffi_abi_token()`, at extension
init), and the versioned capsule names (`molrs.FrameRef/<line>`) from
`molrs_ffi::abi`.
