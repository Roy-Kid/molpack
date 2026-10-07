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
| `io` | `molrs/io` (`molrs::io::{read_frame, write_frame}` for `Script::build` and the CLI) |
| `cli` | `clap` + `io` (the `molpack` binary) |
| `rayon` | `rayon` + `molrs/rayon` (parallel evaluation) |

There is no `ff` feature: the in-loop `Optimizer` seam is always on, and a
caller binding a molrs force-field optimizer (`Lbfgs` over a `Potential`)
enables `ff` on its own molrs dependency. molpack's tests do the same through
a molrs dev-dependency with `ff` on.

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
- **Unit tests only, in-module.** Every test lives in a `#[cfg(test)]` module
  next to the code that owns the behaviour (law § 11); a test body that is too
  large for its file gets a child `tests.rs` / `tests/` module
  (`src/grow/tests/`, `src/pipeline/tests.rs`,
  `src/context/pack_state/tests.rs`, `src/restraint/geometric/tests/`); fixtures
  shared across modules live in `src/testutil.rs` (`cfg(test)`). There
  is **no** `tests/` directory, **no** `benches/` and **no** `regressions/` —
  they were deleted on 2026-09-20 with the end-to-end packing suites, the
  Packmol regression harness and the criterion benches. Reintroducing any of
  them is a design decision, not a convenience.
- **No golden / bitwise-continuity tests.** A test asserts a property the code
  owns (a named refusal, a layering rule, an exact geometric identity), never
  "the numbers this build happened to produce".
- Ownership examples: `src/euler.rs` (rotation algebra), `src/stage.rs` (the
  `Stage` seam, on fake stages only — what a concrete stage declares stays with
  that stage), `src/pipeline/tests.rs` (the multi-stage lifecycle body),
  `src/invariant.rs` (`Layers` / `Invariant` / `RestraintsSatisfied`),
  `src/grow/tests/{internal,field,prior,entry}.rs`, `src/target.rs` +
  `src/script/build.rs` (the four per-atom properties, API side and `.inp`
  side), `src/context/build.rs` (what context construction refuses).
- `cargo test -p molcrafts-molpack --lib --features cli,rayon` — the gate,
  must always be green (`mol_project.build.test`); seconds, not minutes. A
  single test: append `-- <name filter>` (`build.test_single`).
- `cargo test --doc --features cli` — the rustdoc examples, part of the same
  tier. That feature set is the one doctests are written for: a snippet that
  needs `io` or `ff` is a plain `no_run` block, with any setup on hidden `#`
  lines.
- An ignored test is a failing test. No `ignore` doctest fence, no `#[ignore]`:
  the `no-ignored-tests` prek hook (commit stage, mirrored in CI `lint`) rejects
  both. Uncompiled illustration is a `text` fence.
- `uv run --directory python --group dev tox -e py` — Python wheel, isolated and non-editable (`maturin develop` is not the gate).
- CI parity (pre-push hooks in `.pre-commit-config.yaml`, mirrored by
  `mol_project.ci.local`): lib tests, doc tests, `--all-targets` check (bin +
  examples must compile), `--no-default-features` / `rayon` checks, then tox.

## Repo layout

| Path | Purpose |
|---|---|
| `src/` | library + CLI binary (`src/bin/molpack/`) |
| `python/` | PyO3 wheel (`python/src/`) + package (`python/python/molpack/`) + tests |
| `examples/` | runnable example programs (need `--features io`) |
| `docs/` | public docs site (Rust guide + `docs/python/` binding docs) — published via Zensical (`zensical.toml`) with the shared `molcrafts` theme |
| `.claude/specs/` | active feature specs, indexed in `INDEX.md`; deleted on close |
| `.claude/notes/` | passive knowledge: `law.md`, `conventions.md`, `architecture.md`, `notes.md` |

Template geometry is read with `molrs::core::Frame::coords` (Å); the crate-root leaf
`src/template.rs` holds only `coord_rows` (array → `[x, y, z]` rows) and the
rotatable-bond policy. Bond graphs are `molrs::core::Topology`; molpack does not re-export that type.

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
`molcrafts-molrs` crate); `core` is always-on, while `io` is pulled in
through the matching molpack feature above.

**Build cache:** the committed `.cargo/config.toml` routes every build (root
workspace and `python/`) into `../molrs/target`, shared with the sibling molrs
checkout — molrs compiles once per (rustc, features, profile) across both
repos. `rust-toolchain.toml` matches molrs's so the cache fingerprints one
rustc. CI caches that dir and runs sccache.

**molrs ABI line** (the rule itself is law P10): molpack exchanges `molrs_ffi`
handle capsules with the installed `molcrafts-molrs` wheel; both must embed the
same molrs **major.minor** (minor line = ABI version — see molrs
`docs/interop.md`). Gates: `interop::check_abi` (`molrs._ffi_abi_token()`, at
extension init — the one import-time version check) and the versioned capsule
names (`molrs.FrameRef/<line>`, `molrs.RegionRef/<line>`) from `molrs_ffi::abi`.
