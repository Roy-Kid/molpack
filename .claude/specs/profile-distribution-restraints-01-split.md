---
title: Split src/restraint.rs into a src/restraint/ module (no behaviour change)
status: done
created: 2026-06-12
---

# Split src/restraint.rs into a src/restraint/ module (no behaviour change)

## Summary
`src/restraint.rs` has reached 976 LOC, exceeding the 800-line hard cap in CLAUDE.md. This spec is a pure structural refactor that converts the single file into a `src/restraint/` directory: the `Restraint` trait, its doc comments, and the blanket `impl Restraint for Box<dyn Restraint>` move to `src/restraint/mod.rs`, while the 14 concrete geometric restraint types (Packmol kinds 2–15) move verbatim into `src/restraint/geometric.rs`. `mod.rs` re-exports the geometric types with `pub use`, so every external path such as `crate::restraint::InsideSphereRestraint` and the crate-root re-export `molpack::InsideSphereRestraint` resolve exactly as before. There is no behaviour change, no public-API change, and no new restraint. This split unblocks the rest of the profile-distribution-restraints chain, which adds a `ProfileRestraint` family under `src/restraint/profile/` while keeping every file within the 200–400 LOC budget.

## Design
Entities touched: the `restraint` module only. After the refactor the module is a directory with two files.

- `src/restraint/mod.rs` (the thin trait + re-export surface) owns:
  - the module-level `//!` documentation currently at `src/restraint.rs:1–21`,
  - `use molrs::types::F;`,
  - the `pub trait Restraint: Send + Sync + std::fmt::Debug` definition (`src/restraint.rs:46–62`, signatures `fn f(&self, x: &[F;3], scale: F, scale2: F) -> F` and `fn fg(&self, x: &[F;3], scale: F, scale2: F, g: &mut [F;3]) -> F`, plus the defaulted `is_parallel_safe`, `name`, `periodic_box`),
  - the blanket `impl Restraint for Box<dyn Restraint>` (`src/restraint.rs:64–86`),
  - `mod geometric;` and a single `pub use geometric::{ ... }` block listing all 14 concrete types.
- `src/restraint/geometric.rs` owns the 14 concrete structs and their `impl Restraint` blocks, moved verbatim (kinds 2–15: Inside/Outside Cube, Inside/Outside Box, Inside/Outside Sphere, Inside/Outside Ellipsoid, Above/Below Plane, Inside/Outside Cylinder, Above/Below Gaussian). It begins with `use crate::restraint::Restraint;` and `use molrs::types::F;` so the impls resolve the trait through the unchanged `crate::restraint` path.

Ownership / lifecycle: no types are renamed, no fields change, no signatures change. `RegionRestraint<R>` stays in `src/region.rs` (its `impl Restraint` at `region.rs:230–248` already imports the trait via `use crate::restraint::Restraint;` at `region.rs:35`, a path that is identical before and after the split). The five internal consumers — `src/validation.rs:8`, `src/target.rs:7`, `src/region.rs:35`, `src/packer.rs:25`, `src/context/pack_context.rs:7` — all name `crate::restraint::Restraint`, which `mod.rs` continues to export, so none of them change. The crate-root `pub use restraint::{ ... }` block in `src/lib.rs:130–135` is unchanged because `mod.rs` re-exports the same names.

Immutability and the other CLAUDE.md hard rules are preserved trivially: no logic is rewritten, no "packmol" identifier is introduced, and both new files sit within the LOC budget (geometric.rs carries the bulk, under the 800-max; mod.rs is a thin trait + re-export surface).

> **Chain bootstrap note (for the implementer):** this sub-spec creates `src/restraint/mod.rs` and `src/restraint/geometric.rs` only. The first sibling sub-spec to add a `profile/` file (`-02-coordinate`) is responsible for creating `src/restraint/profile/mod.rs` and adding the `pub mod profile;` line to `src/restraint/mod.rs`; later siblings only append `pub mod <name>;` to the existing `profile/mod.rs`.

## Files to create or modify
- `src/restraint/mod.rs` (new) — trait, blanket impl, module docs, `mod geometric;`, `pub use geometric::{...}` re-exports.
- `src/restraint/geometric.rs` (new) — the 14 concrete geometric restraint structs and their `impl Restraint` blocks, moved verbatim.
- `src/restraint.rs` — removed (its contents are partitioned into the two files above; converting the file to a directory module removes the standalone file).

## Tasks
- [x] Run `cargo test -p molcrafts-molpack --lib --tests` to capture the green baseline before any change (RED-baseline: confirm `tests/restraint.rs` and unit tests pass against the current single-file layout).
- [x] Create `src/restraint/` directory and add `src/restraint/geometric.rs` (new) by moving the 14 concrete restraint structs and their `impl Restraint` blocks (`src/restraint.rs:88–976`) verbatim, prefixed with `use crate::restraint::Restraint;` and `use molrs::types::F;`.
- [x] Add `src/restraint/mod.rs` (new) holding the module docs, the `Restraint` trait, the blanket `impl Restraint for Box<dyn Restraint>`, `mod geometric;`, and a `pub use geometric::{...}` block re-exporting all 14 geometric types under their original names.
- [x] Delete the now-empty `src/restraint.rs` so the `restraint` module resolves to the new directory.
- [x] Verify `crate::restraint::<GeometricType>` and `molpack::<GeometricType>` paths are unchanged by running `cargo clippy -p molcrafts-molpack --all-targets -- -D warnings` and confirming zero import churn outside the `restraint` module (only `mod.rs`/`geometric.rs` differ).
- [x] Run `cargo fmt --all` and confirm `src/restraint/mod.rs` and `src/restraint/geometric.rs` each stay within the 800-LOC hard cap.
- [x] Run full check + test suite (`cargo fmt --all -- --check`, `cargo clippy -p molcrafts-molpack --all-targets -- -D warnings`, `cargo test -p molcrafts-molpack --lib --tests`) to prove zero behaviour change with the existing test files unchanged.

## Testing strategy
- Happy path: the existing `tests/restraint.rs` regression suite (including `assert_gradient_opposes_violation` at `tests/restraint.rs:24–43`) is the unchanged net — it imports all 14 types plus `Restraint` from `molpack::` and must pass byte-for-byte identical assertions after the split, proving value (`f`) and gradient (`fg`) behaviour is preserved.
- Edge cases: confirm no import path outside the `restraint` module is touched — the five internal consumers (`validation.rs`, `target.rs`, `region.rs`, `packer.rs`, `context/pack_context.rs`) and the `lib.rs` crate-root re-export must compile without edits, and `RegionRestraint<R>` in `region.rs` must still resolve `crate::restraint::Restraint`.
- Build integrity: `cargo clippy -- -D warnings` clean across `--all-targets` (lib, tests, benches) confirms no orphaned references and no dead re-exports.
- Domain validation: not applicable — this is a structural move with no numerical change; Packmol parity is inherited unchanged because the `impl` bodies are moved verbatim.

## Out of scope
- Any new restraint type, including the `ProfileRestraint` family added by later sub-specs of this chain.
- Any change to a `Restraint` trait method signature or to a concrete restraint's fields, constructors, or numerics.
- Moving or modifying `RegionRestraint<R>` (stays in `src/region.rs`).
- Any edit to `src/objective.rs`, `src/script/`, or any consumer's import statements.
- Splitting `geometric.rs` further into per-kind files (deferred; this spec only lifts the trait surface out so subsequent profile work has room).
