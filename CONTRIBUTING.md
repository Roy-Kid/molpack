# Contributing to molcrafts-molpack

Thanks for your interest in contributing. This document covers how to set up a
development environment, run tests, and get a PR merged.

## Development environment

**Prerequisites:**
- Rust 1.85+ (`rustup update stable`)
- Python 3.9+ with `maturin` and `pytest` for the Python bindings
- The [molrs](https://github.com/MolCrafts/molrs) repo checked out as a sibling:

```
workspace/
├── molrs/      ← git clone https://github.com/MolCrafts/molrs
└── molpack/    ← this repo
```

The root `Cargo.toml` uses a path dependency on `../molrs/molrs-core`. With the
sibling layout above everything resolves automatically.

**First-time setup:**

```bash
# Fetch test data used by the regression suite
bash ../molrs/scripts/fetch-test-data.sh

# Verify the Rust build
cargo build -p molcrafts-molpack

# Verify the Python bindings
cd python && maturin develop --release && pytest
```

## Running tests

```bash
# Fast: unit + integration (no test data required)
cargo test -p molcrafts-molpack --lib --tests

# Full: includes examples compilation
cargo test -p molcrafts-molpack

# Packmol regression suite (requires test data, slow)
cargo test -p molcrafts-molpack --release --test examples_batch -- --ignored

# Python bindings
cd python && pytest -v
```

## Code style

- `cargo fmt` — enforced by CI, run before committing
- `cargo clippy -- -D warnings` — no warnings allowed
- Follow the immutability rule: return new values, never mutate in place
- Keep files under ~400 lines; split at ~200 if the module grows beyond one concern
- New public types must implement `Debug` and, where appropriate, `Clone`

## Adding a new restraint type

1. Add a struct to `src/restraint.rs` with semantically-named fields
2. Implement `Restraint` (both `f` and `fg`; `fg` must match the gradient of `f`)
3. Export it from `src/lib.rs`
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

## Pull request checklist

All items in the PR template must be checked before requesting review. CI must
be green. For restraint or objective changes, the Packmol regression suite
(`examples_batch`) must also pass.

## Licensing

By contributing you agree that your contributions will be licensed under
the BSD-3-Clause license that covers this project.
