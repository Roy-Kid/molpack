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
| `src/script/` | `.inp` parsing + lowering to Targets (`parser.rs`, `build.rs`, `io.rs`) |
| `src/restraint.rs` | `Restraint` trait + concrete restraints (analytic `f` + matching `fg`) |
| `src/objective.rs` | objective + gradient construction over targets |
| `src/packer.rs` | outer loop driver (move / relax / evaluate) |
| `src/gencan/` | bound-constrained optimizer (cg, spg) |
| `src/region.rs`, `src/cell.rs` | neighbor lookup |
| `src/initial.rs`, `src/relaxer.rs`, `src/movebad.rs` | initialization, relaxation, escape moves |
| `src/bin/molpack/` | CLI front-end (cli feature) |
| `python/src/` | PyO3 wheel |

## Cargo features

| Feature | Pulls in |
|---|---|
| `default` | nothing |
| `io` | `molrs_io` (PDB / XYZ / SDF / LAMMPS readers) |
| `cli` | `clap` + `io` (binary + integration tests) |
| `rayon` | `rayon` + `molrs/rayon` (parallel evaluation) |

The Python wheel is built **without** `io` — the wheel relies on the user's `molrs` Python package
for frame loading, then calls `Target::from_frame_parts` / `Script::lower`.

## Hard rules (load-bearing)

- **No "packmol" in any public symbol or module name.** The product is `molpack`. The script format happens to be Packmol-compatible. Doc-comment prose may mention Packmol; identifiers may not.
- **Configuration is Packmol `.inp` only.** Do not invent TOML / YAML configs.
- **molrs path + version pins are managed manually.** Do not automate the check in pre-commit / CI.
- **Pre-commit and CI use first-party tooling.** No bespoke `scripts/check-*.sh` glue. Prefer registry-hosted hooks (`doublify/pre-commit-rust`, `astral-sh/ruff`).
- **Git workflow:** fork → PR. Never push directly to `MolCrafts/molpack` master. `origin` = Roy-Kid fork, `upstream` = MolCrafts.

## Coding style

- Rust 1.85+, edition 2024
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
- `cargo bench --bench pack_end_to_end -- mixture` — perf baseline.

## Repo layout

| Path | Purpose |
|---|---|
| `src/` | library + CLI binary |
| `python/` | PyO3 wheel + tests + docs |
| `tests/` | Rust integration tests (incl. `examples_batch.rs` regression suite) |
| `benches/` | criterion benchmarks |
| `examples/` | runnable example programs (need `--features io`) |
| `docs/` | concepts / architecture / extending — published via Zensical |
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

The root `Cargo.toml` uses path deps on `../molrs/molrs-core` and `../molrs/molrs-io`.
