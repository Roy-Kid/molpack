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

**Version pins:** a local build uses the sibling checkout `../molrs`. The
version fields name the **0.16.*** minor line (see `Cargo.toml` and
`python/pyproject.toml`). CI, and the pre-push hooks, build against
`MolCrafts/molrs` at the commit `scripts/partners.py` resolves (molrs's
`dev`; see [Partners](#partners)), in a layout of their own -- not your
sibling. `Cargo.lock`, `python/Cargo.lock` and `python/uv.lock` are committed
and every cargo / uv call in the gates is `--locked`; they record molrs's
package metadata, not a commit, so relock only when molrs's `dev` changes its
version or dependencies (the recipe is in `.github/partners.env`).
The Python wheel checks this on ``import molpack`` — a molrs minor mismatch is an
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

prek pre-push and CI run the **same** commands (spelled out in
`.pre-commit-config.yaml` and `.github/workflows/`):

```bash
# Rust — same as CI `test / rust`. Behaviour lives in `#[cfg(test)]`
# modules next to the code; there is no tests/ directory.
cargo test --locked --lib --features cli
cargo test --locked --doc --features cli
cargo check --locked --all-targets --features cli      # bin + examples must compile
cargo check --locked --no-default-features && cargo check --locked --features rayon

# Python — tox isolated env (tox itself from python/ dependency-group dev)
uv run --locked --python 3.12 --directory python --group dev tox -e py
```

Against your sibling checkouts these are a quick loop; the gate is the same
command in CI's layout: `scripts/partners.py run -- <command>`.

For a quick local edit loop you may still `maturin develop` in a personal
venv; **do not** rely on that for the gate — pre-push always uses
`uv run --directory python --group dev tox -e py`.

## Hooks

`prek install` installs both hook types from `.pre-commit-config.yaml`; every
command in `.github/workflows/{lint,test,docs}.yml` has a hook. **Never `git commit
--no-verify` or `git push --no-verify`, and never merge a red pull request.**

- **pre-commit** (staged files, nothing compiles, in place): file hygiene
  (whitespace, final newline, YAML/TOML, merge markers, line endings), ruff
  format + ruff (`python/`), rustfmt, and the no-ignored-tests guard.
- **pre-push**:
  - the pre-commit hooks again on `--all-files` (CI `lint / hooks` runs them so);
  - `scripts/partners.py check` — every partner in `.github/partners.env`
    resolves, every path dependency (`../molrs`, `../../molrs`, ...)
    lands in a checkout CI makes, and no workflow spells a partner ref of its
    own;
  - the three lock files are current against the resolved molrs;
  - the docs build (`zensical build --clean --strict` from the `doc` group)
    when docs/ or zensical.toml changed;
  - clippy (`--all-targets --all-features -D warnings`), ty, the rust tests
    and tox -- each in CI's sibling layout (`scripts/partners.py run`: a copy
    of this tree next to molrs at its resolved commit), `--locked`,
    on the toolchain `rust-toolchain.toml` pins (1.99.0, the same as molrs)
    and Python 3.12.
- **Dispatch on the MolCrafts cluster:** the compiling gates go through
  `scripts/hook-run.sh`, which hands the command to `$MOLCRAFTS_HOOK_RUNNER`
  when that is set and it is not already inside a Slurm job. The cluster's
  shared `core.hooksPath` sets it to `.build-alloc/hookrun`, which runs the
  command on a compute node (allocation `$USER-hooks`; it fails after 20
  minutes without a node, never passes) and sets `$MOLCRAFTS_PARTNER_CACHE` so
  the partner checkouts and their `target/` stay warm between pushes.
  Everything else runs in place, so a commit never waits for Slurm. Elsewhere
  nothing sets the variables: every hook runs locally, in a temp layout.

## Partners

On `dev`, partners are tracked, not pinned. `.github/partners.env` names
molrs's branch (`MOLRS_REF=dev`), and `scripts/partners.py` resolves it -- for
CI (`partners.py fetch`, into `../molrs`) and for the hooks (`partners.py
run`) alike -- to the first of:

1. molrs's branch named like the one being built (CI: the pushed branch or a
   pull request's head branch; locally: the checked-out branch), looked up
   first on the fork the build comes from (`<owner>/molrs`, where `<owner>`
   owns the pull request's head repository or the repository CI runs in; in a
   git hook, the remote being pushed to), then on MolCrafts/molrs;
2. outside CI only, that branch in your sibling clone `../molrs`, when it has
   one and neither remote does yet;
3. MolCrafts/molrs's `dev`.

A change that needs a molrs change lands as two same-named branches, never by
skipping a gate: create the same branch (say `converge/x`) in both checkouts;
push both to your forks, never to MolCrafts (molpack's gates take molrs's
branch from your fork, or from your sibling before it is pushed); each push
runs the full CI tier on your fork, and molpack's run resolves molrs's
`converge/x` there; only once both forks are
green, open the pull requests into MolCrafts `dev`, land molrs's, then
molpack's (never a red one), and delete the branches.

**Releasing.** A release builds against a fixed molrs: the release commit on
`master` sets `MOLRS_REF` to the molrs release tag of the line `Cargo.toml`
names (`vX.Y.Z`), relocks if that changes molrs's metadata, and is tagged;
`release.yml` resolves that tag. When `master` is merged back into `dev`, keep
`MOLRS_REF=dev` there.

## CI

One workflow per kind of work. Every push of any branch runs `lint`, `test`
and `docs`, on a fork as on MolCrafts. A pull request into `dev` or `master`
runs them again unless it is a pull request inside a fork (that branch was
already built, full tier, by its push). Those decisions (tier, fork or
upstream, the duplicate pull request) are made in one place: every
workflow's first job, `<file> / context`, runs
`MolCrafts/molcrafts-ci/actions/ci-context@master`, and every other job reads
its outputs.

| workflow | feature-branch push to MolCrafts | everything else: `dev`/`master`/`main` on MolCrafts, pull requests, tags, dispatches, any push to a fork | upstream only |
| --- | --- | --- | --- |
| `lint.yml` | `lint / hooks` (commit hooks on every file, partners, lock files), `lint / clippy` (clippy, ty), `lint / workflows` (`check-workflows`) | same | — |
| `test.yml` | fast: `test / rust`, `test / python (ubuntu-latest)` | full: `test / rust`, `test / python` on Linux, macOS and Windows | — |
| `docs.yml` | `docs / build` (zensical `--strict`) | same | Cloudflare Pages deploys the site from MolCrafts |
| `release.yml` | — | dispatch: dry run (gates, builds, `release / crate (dry run)`: `cargo publish --dry-run`; the upload jobs are skipped) | `v*` tag (`publish` from `release / context`): crates.io, PyPI wheels + sdist, GitHub Release |

So a fork branch gets the full tier on its push: push to your fork, wait for
green, then open the pull request into MolCrafts `dev`. Branches pushed to
MolCrafts itself (Dependabot's) get the fast tier, and their pull requests the
full one. The `require-green-ci` (`dev`) and `protect-master` rulesets require
`test / context` and the full tier's jobs. A release tag must be `v` + the `Cargo.toml` version, on
`master`; trusted publishing on crates.io and PyPI names `release.yml` and the
`crates-io` and `pypi` environments.
Shared setup is MolCrafts/molcrafts-ci's `actions/setup-rust`,
`actions/setup-python` and `actions/setup-partners` (`@master`), the same
actions every MolCrafts repository uses.

## Code style

- `cargo fmt` / `ruff` — commit-stage prek hooks; `clippy` / `ty` — pre-push
- pre-push: clippy, ty, rust tests + `uv run … tox -e py` (mirrors CI; see [Hooks](#hooks))
- Follow the immutability rule: return new values, never mutate in place
- Keep files under ~400 lines; split at ~200 if the module grows beyond one concern
- New public types must implement `Debug` and, where appropriate, `Clone`

## Adding a new restraint type

1. Add a struct under `src/restraint/` with semantically-named fields — alongside the existing geometric restraints in `src/restraint/geometric/`, or as a new submodule (see `src/restraint/profile/` for a composed example)
2. Implement `Restraint` (both `f` and `fg`; `fg` must match the gradient of `f`)
3. Re-export it from `src/restraint/mod.rs`, then from the crate root in `src/lib.rs`
4. Add a `#[cfg(test)]` unit test next to the struct — the restraint owns its own tests
5. Document it in `docs/concepts.md` under the restraint table

See the `extending` rustdoc chapter (`cargo doc --open`) for detailed tutorials.

## Commit messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>: <short description>

<optional body>
```

Types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`, `ci`
