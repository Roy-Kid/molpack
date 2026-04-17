# Deferred renames (second iteration)

Internal names that do not match the spirit of `docs/api_redesign.md` but
are **not touched in this pass** — they live behind the public surface,
or would balloon the scope of the current PR chain. Each row records:
what, where, why deferred.

A second rename iteration (after the public surface lands) will work
through this list.

| Item | Location | Kind | Rationale to defer |
|---|---|---|---|
| `molcrafts-molrs-core = "0.0.8"` / `molcrafts-molrs-io = "0.0.8"` | `Cargo.toml:24,25` | dep version pin | Bumped to 0.0.9 during Stage 1b to unblock `cargo check`; local molrs workspace is already on 0.0.9. Pre-existing drift — confirm with user whether to keep bump. |
| 17 doctests in `docs/*.md` | `docs/getting_started.md`, `docs/extending.md`, `docs/index.md`, plus `README.md` | rustdoc code blocks | Use old 2-arg `InsideBoxRestraint::new`, old `.tolerance()` / `.pbc()` / `.pack(..,seed)`. Scheduled for Stage 6 (docs rewrite). |
| `PackContext` fields `pbc_min`, `pbc_length`, `pbc_periodic` still pub | `src/context/pack_context.rs` | struct field | Kept pub for test / bench access; spec §2 keeps PackContext internals off the supported surface. No rename needed. |
| `MoveBadConfig.{movefrac, movebadrandom, gencan_maxit}` | `src/movebad.rs` | internal struct fields | Internal, referenced by `run_phase` / `run_iteration`. Renames would cascade into several files without user-visible value. |
| `GencanParams.{maxit, maxfc, iprint}` | `src/gencan/mod.rs` | internal struct fields | Internal solver parameters; Packmol-faithful names. Not on the user-facing surface. |
| `initial()` param `sidemax` | `src/initial.rs:302` | fn param | Internal; called only from `packer::pack`. Rename when touching initial.rs next. |
| `initial()` param `nloop0` alias | — | — | Already fixed to `init_passes` internally in packer.rs; param name on `initial()` still `nloop0`. Defer internal. |

## Conventions for adding rows

- **Item** — the exact identifier (e.g. `compute_f`, `sys.fdist`, `PBCBox`).
- **Location** — `file.rs:line` or module path.
- **Kind** — `internal fn` / `struct field` / `internal type` / `var name` / `test helper`.
- **Rationale to defer** — one sentence. Don't omit; if there's no
  rationale, rename in the current stage instead.
