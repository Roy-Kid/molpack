# molpack — project conventions

## What this is

A faithful Rust port of Packmol with PyO3 Python bindings and a `.inp`-compatible CLI.
Lives at `/Users/roykid/work/molcrafts/molpack`. Three public surfaces:

- `molcrafts-molpack` Rust library (`[lib] name = "molpack"`)
- `molpack` binary (`cli` feature)
- `molcrafts-molpack` Python wheel under `python/`, built with maturin **without** the `io` feature

## Architecture (lifecycle)

```
.inp script ──parser──▶ Script ──build──▶ Targets ──Molpack::pack──▶ PackResult
                                                       │
                                                       └─ uses: restraint, objective, packer, gencan/
```

Modules:

| Path | Owns |
|---|---|
| `src/script/` | `.inp` parsing (`parser.rs`) + lowering to Targets/plans (`build.rs`) + format readers/writers (`io.rs`, `io` feature) + errors (`error.rs`) |
| `src/target.rs` | `Target` builder (coords/radii + restraints + relaxer + centering/orientation); `Axis` / `Angle` / `Placement` / `CenteringMode` |
| `src/restraint/` | `AtomRestraint` trait + per-atom geometric restraints (`geometric/`) and collective distribution-matching restraints (`collective/`) |
| `src/region.rs` | `Region` predicate trait + `And`/`Or`/`Not` combinators + `RegionRestraint` lift + `Aabb` + `InsideCellRegion` (primitive cell, any lattice) |
| `src/objective.rs` | objective value + gradient over targets (pair terms, PBC, rayon paths); `Objective` impl |
| `src/constraints/` | `Constraints` container + `EvalMode` / `EvalOutput` evaluation entrypoint |
| `src/packer.rs` | `Molpack` builder + `pack()` orchestration + inner-loop driver (move / relax / evaluate); `PackResult` |
| `src/gencan/` | bound-constrained optimizer: GENCAN/SPG driver (`mod.rs`), CG (`cg.rs`), SPG (`spg.rs`) |
| `src/context/` | packing runtime `PackContext` (carries the `SimBox` + molrs `CellGrid` and the forward-neighbour table) + model/state views + reusable work buffers |
| `src/initial.rs`, `src/movebad.rs` | initial placement; `movebad` escape moves |
| `src/relaxer/` | in-loop per-molecule relaxers: `Relaxer` / `TorsionMcRelaxer` (`mod.rs`), `LBFGSRelaxer` (`lbfgs.rs`, `ff` feature) |
| `src/handler.rs` | `Handler` trait + built-in progress / log / XYZ / early-stop handlers |
| `src/assemble.rs`, `src/frame.rs` | build topology-complete `molrs::Frame` from packed coords; frame ↔ coords conversion |
| `src/validation.rs` | post-pack validation (`ViolationMetrics`, `validate_from_targets`) |
| `src/euler.rs`, `src/numerics.rs` | Euler-angle ↔ matrix; precision / numeric-control helpers |
| `src/cases.rs` | example / regression harness (`ExampleCase`, `build_targets`, `io` feature) |
| `src/bin/molpack/` | CLI front-end (`cli` feature) |
| `python/src/` | PyO3 wheel (built **without** `io`) |

## Cargo features

| Feature | Pulls in |
|---|---|
| `default` | nothing |
| `io` | `molrs/io` (PDB / XYZ / SDF / LAMMPS readers) |
| `cli` | `clap` + `io` (binary + integration tests) |
| `rayon` | `rayon` + `molrs/rayon` (parallel evaluation) |
| `ff` | `molrs/ff` (MMFF + L-BFGS) — enables `LBFGSRelaxer` |

The Python wheel is built **without** `io` — the wheel relies on the user's `molrs` Python package
for frame loading, then builds targets from the loaded frame (the PyO3
`target_from_frame` helper) and lowers scripts with `Script::lower` +
`StructurePlan::apply`.

## Hard rules (load-bearing)

- **No "packmol" in any public symbol or module name.** The product is `molpack`. The script format happens to be Packmol-compatible. Doc-comment prose may mention Packmol; identifiers may not.
- **Configuration is Packmol `.inp` only.** Do not invent TOML / YAML configs.
- **molrs path + version pins are managed manually.** Do not automate the check in pre-commit / CI.
- **prek + tox for local gates.** Hooks use prek (pre-commit-compatible config).
  No project `scripts/` test wrappers. Python isolation is
  `uv run --directory python --group dev tox -e py` (tox lives in
  `python/` dependency-group `dev`). Prefer registry-hosted hooks
  (`doublify/pre-commit-rust`, `astral-sh/ruff`).
- **Git workflow:** fork → PR. Never push directly to `MolCrafts/molpack` master. `origin` = Roy-Kid fork, `upstream` = MolCrafts.

## Coding style

- Rust 1.91+, edition 2024
- Immutability — return new values, never mutate in place
- Files: 200–400 lines typical, 800 max — split when modules grow beyond one concern
- New public types implement `Debug` and (where appropriate) `Clone`
- `cargo fmt` and `cargo clippy -- -D warnings` are mandatory before commit

## Tests

- TDD: RED first, then GREEN, then refactor. 80% coverage minimum.
- `cargo test -p molcrafts-molpack --lib --tests` — fast tier, must always be green.
- `cargo test -p molcrafts-molpack --release --test examples_batch -- --ignored` — Packmol regression.
  Requires test data: `bash ../molrs/scripts/fetch-test-data.sh` (one time).
- `cd python && maturin develop --release && pytest` — Python wheel.
- `cargo bench --benches` — criterion regression benches (pair kernel, objective dispatch, one `run_iteration` step, restraint eval, tiny end-to-end pack). Each is small/fast and synthesizes its own geometry, so **no** `io` feature is needed; run one with `cargo bench --bench pack_end_to_end`. Perf history is tracked on canonical pushes by `.github/workflows/bench.yml`.

## Repo layout

| Path | Purpose |
|---|---|
| `src/` | library + CLI binary |
| `python/` | PyO3 wheel + tests |
| `tests/` | Rust integration tests (incl. `examples_batch.rs` regression suite) |
| `benches/` | small criterion regression benches (self-synthesized geometry, no `io`) |
| `examples/` | runnable example programs (need `--features io`) |
| `docs/` | full docs site (Rust guide + `python/` binding docs) — published via Zensical with the shared `molcrafts` theme |
| `.claude/specs/` | feature specs, indexed in `INDEX.md` |
| `.claude/NOTES.md` | evolving decisions; promoted to CLAUDE.md when stable |
| `.claude/skills/` | user-invocable workflows (`/mpk-*`) |
| `.claude/agents/` | single-axis review agents — invoked only by skills |

## Sibling layout assumed

```
workspace/
├── molrs/      ← git clone https://github.com/MolCrafts/molrs
└── molpack/    ← this repo
```

The root `Cargo.toml` uses a single path dep on `../molrs/molrs` (the unified
`molcrafts-molrs` crate); `core` is always-on, while `io` and `ff` are pulled in
through the matching molpack features above.
