# Contributing to molcrafts-molpack

Thanks for your interest in contributing. This document covers how to set up a
development environment, run tests, and get a PR merged.

## Development environment

**Prerequisites:**
- Rust 1.91+ (`rustup update stable`)
- Python 3.12+ with `maturin` and `pytest` for the Python bindings (the wheel's `requires-python` is `>=3.12`)
- [prek](https://github.com/j178/prek) for git hooks (drop-in for pre-commit)
- [uv](https://docs.astral.sh/uv/) (hooks run `uv run --group dev tox …`; tox is in
  `python/` dependency-group `dev`, not a separate global install)
- The [molrs](https://github.com/MolCrafts/molrs) repo checked out as a sibling:

```
workspace/
├── molrs/      ← git clone https://github.com/MolCrafts/molrs
└── molpack/    ← this repo
```

The root `Cargo.toml` uses a path dependency on `../molrs/molrs`. With the
sibling layout above everything resolves automatically.

**Version pins:** path molrs / PyPI `molcrafts-molrs` / `molcrafts-molpy` track
the **0.14.*** minor line (see `Cargo.toml`, `python/pyproject.toml` `[tool.tox]`,
and `MOLRS_GIT_REF` in `.github/workflows/ci.yml`). Patch may differ; only
major.minor must match. Keep the sibling molrs clone on that minor line
(`git checkout v0.14.x` or the matching release branch). The Python wheel
checks this on ``import molpack`` — a molrs/molpy minor mismatch is an
``ImportError``, not a later FFI segfault.

**First-time setup:**

```bash
# Fetch test data used by the regression suite
bash ../molrs/scripts/fetch-test-data.sh

# Git hooks via prek (reads .pre-commit-config.yaml)
prek install

# Verify the Rust build
cargo build -p molcrafts-molpack
```

## Running tests

prek pre-push and CI run the **same** commands (no project `scripts/` wrappers):

```bash
# Rust — same as CI "rust tests" job
cargo test --features io --lib --tests --examples
cargo test --test cli --features cli
cargo check --no-default-features && cargo check --features rayon && cargo check --features ff
cargo build --benches

# Python — tox isolated env (tox itself from python/ dependency-group dev)
uv run --directory python --group dev tox -e py

# Packmol regression (slow; CI master only)
cargo test --release --features io --test examples_batch -- --ignored
```

For a quick local edit loop you may still `maturin develop` in a personal
venv; **do not** rely on that for the gate — pre-push always uses
`uv run --directory python --group dev tox -e py`.

## Code style

- `cargo fmt` / `clippy` / `ruff` / `ty` — commit-stage prek hooks
- pre-push: rust tests + `uv run … tox -e py` (mirrors CI)
- Follow the immutability rule: return new values, never mutate in place
- Keep files under ~400 lines; split at ~200 if the module grows beyond one concern
- New public types must implement `Debug` and, where appropriate, `Clone`

## Adding a new restraint type

1. Add a struct under `src/restraint/` with semantically-named fields — alongside the existing geometric restraints in `src/restraint/geometric/`, or as a new submodule (see `src/restraint/profile/` for a composed example)
2. Implement `Restraint` (both `f` and `fg`; `fg` must match the gradient of `f`)
3. Re-export it from `src/restraint/mod.rs`, then from the crate root in `src/lib.rs`
4. Add a unit test in `tests/restraint.rs`
5. Document it in `docs/concepts.md` under the restraint table

See the `extending` rustdoc chapter (`cargo doc --open`) for detailed tutorials.

## Commit messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>: <short description>

<optional body>
```

Types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`, `ci`
