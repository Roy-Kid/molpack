# rust-review-improvements

Status: **DONE** — all tiers implemented (one documented exception: the
`tn_linesearch` half of #8 was deliberately left cohesive; see §7 note).
Each item verified with Packmol `examples_batch` parity + 122 pytest +
clippy (`--features cli` and `--no-default --features rayon`), committed
individually on branch `rust-review-improvements`.
Origin: formalized from upstream branch `claude/rust-project-review-eajU5`
(`docs/rust_review_spec.md`), reconciled against the current working tree.
Scope: Rust core (`src/`) + PyO3 wheel (`python/src/`).

> This is an **umbrella spec** — a tiered backlog of refinements, not one
> feature. Each tier is independently shippable. Implement tier by tier via
> `/mpk-impl rust-review-improvements` (or cherry-pick a single item). Nothing
> here regresses Packmol parity.

## Goal

Reduce cognitive load and close robustness gaps in a high-quality
Fortran→Rust port, **without changing packing results**. Three thrusts:

1. **Robustness** — stop the FFI boundary and `.inp` parser from panicking or
   silently clamping on malformed input; centralize validation.
2. **Maintainability** — collapse duplicated hot kernels, de-Fortran the
   naming (with traceability comments), name magic constants, split the two
   monster functions.
3. **Architecture / extensibility** — decompose the `PackContext` god-object,
   type the phase loop, and unify the restraint surface so adding a restraint
   touches one source of truth instead of three.

Observable outcome: identical numerics and CLI/Python behavior on valid input;
clean typed errors (not panics) on invalid input; more restraint kinds
reachable from `.inp`/Python; smaller, better-named modules.

## Already done (excluded from open scope)

Verified present in the current staged working tree — do **not** re-spec:

- `GencanWorkspace` preallocates `lastgpns` + `bounds_l`/`bounds_u`; `pgencan`
  reuses them via `mem::take` swap (`src/gencan/mod.rs`).
- Relaxer MC loop reuses a hoisted `coords_snapshot` buffer instead of
  `coor[..].to_vec()` per iteration (`src/packer.rs`); `TorsionMcRelaxerRunner`
  hoists the `trial` buffer + `copy_from_slice` instead of `clone()`
  (`src/relaxer.rs`).
- `project_cartesian_gradient` uses `g.fill(0.0)` (`src/objective.rs`).
- `render_packmol_input` → `render_molpack_input`; `cell` module is now
  `pub(crate)`; dead `delta_vector` removed (`src/cases.rs`, `src/lib.rs`,
  `src/cell.rs`) — satisfies the no-"packmol"-in-public-symbols rule.

## Non-goals

- No new packing algorithm, restraint kind, or optimizer behavior.
- No `unsafe`. No change to the `Objective`/`Restraint`/`Relaxer` trait
  *semantics* (signature changes for `Handler` are in scope; numerics are not).
- No TOML/YAML config — `.inp` only.
- Not re-doing the perf items already landed (see above).
- The published `docs/` site is **not** the home for this plan (this spec is).

## Public surface

| Change | Surface | Kind |
|---|---|---|
| 7–8 `.expect()` → `PyValueError` via `?` | `python/src/packer.rs` | behavior on bad input (panic → typed error) |
| `Handler` callbacks return `Result` (or RAII) | `src/handler.rs` trait | **breaking** trait-impl change (pre-1.0) |
| `enum Phase { PerType(usize), AllTypes }` | `src/packer.rs` | internal only |
| `IterationState` struct (run_iteration 19→~4 params) | `src/packer.rs` | internal only |
| Restraint registry (kind→keyword→Rust ctor→Py class) | `src/restraint.rs` / new `restraint/registry.rs` | additive; more `.inp` keywords + Python classes reachable |
| De-Fortran renames | private fns/fields only — keep `// Packmol: <old>` | non-breaking (no public symbol renamed) |
| Drop `Constraints` ZST | remove `pub mod constraints` from `lib.rs` | **breaking** if anything re-exports it — confirm first |
| Reject `idx < 1`, radius `> 0`, `count > 0` | `src/validation.rs` (+ parser/build error paths) | behavior on bad input |

## Module placement

From the `mpk-architect` placement note:

- **`PackContext` split** → keep `PackContext` as composer in
  `context/pack_context.rs`; put `Topology` / `CellGrid` / `EvalState` /
  `PbcParams` each in its own file under `src/context/`, re-exported from
  `context/mod.rs`. (Mandatory — `pack_context.rs` is already 908 lines.)
- **`IterationState`** → new file; promote `src/packer.rs` → `src/packer/mod.rs`
  and add `packer/iteration_state.rs` (packer.rs is already 1069 lines).
- **`Phase` enum** → `src/packer/phase.rs`, co-located with the `0..=ntype` loop.
- **Drop `Constraints` ZST** → inline `EvalMode` match into
  `PackContext::evaluate` in `src/context/`; remove `src/constraints/`.
- **Merge 5 pair kernels** → stays inside `src/objective.rs` (net line
  *reduction*).
- **`evaluate_unscaled` radius RAII guard** → private, in `src/packer/`.
- **Split `pack()`** (`validate_inputs`/`prepare_system`/`run_main_loop`) →
  `src/packer/`; **split `tn_linesearch`** (`tn_extrapolation`/
  `tn_interpolation`) → its own file under `src/gencan/`.
- **Magic constants** → module-level `const` in owning file: movebad `10.0` →
  `src/movebad.rs`; CG `GAMMA`/`THETA` → `src/gencan/cg.rs`; kernel `4.0` /
  `fimp` clamp → `src/objective.rs` / `src/packer/`.
- **Restraint registry** → single source of truth in `src/restraint.rs` (or
  `src/restraint/registry.rs`) — the **lowest** layer. `src/script/parser.rs`
  and `python/src/` *consume* it. **Do not** place it in `script/` (would
  force `python/src/`, which bypasses the parser, to depend upward on `script`).
- **Centralized validation** → extend existing `src/validation.rs` (already
  imports only lower layers `restraint`/`target`/`numerics`; already consumed
  by all three frontends). Do not scatter into `parser.rs`/`python/src/`.
- **`Handler` `Result`/RAII** → `src/handler.rs` trait, consumed by `packer`.

## Numerical contract

- **Invariant: zero numeric change on valid input.** Every tier must keep
  `tests/examples_batch.rs` (the Packmol regression suite) bit-for-bit green.
- Kernel merge (§ pair-kernel) is the highest-risk numeric item: the merged
  parameterized kernel must reproduce `fparc`/`gparc`/`fgparc`/`fgparc_stats`/
  `fparc_stats` exactly. Gate with the `pair_kernel` criterion bench (perf) and
  a gradient finite-difference check (correctness).
- De-Fortran renames + magic-constant extraction must be pure substitutions —
  no value changes (e.g. `4.0`, `GAMMA=1e-4`, `THETA=1e-6`, `fimp ∈
  [-99.99, 99.99]` clamp preserved exactly).

## Test plan

- **Fast tier** (`cargo test -p molcrafts-molpack --lib --tests`) green after
  every item.
- **Parity** (`cargo test --release --test examples_batch -- --ignored`) green
  after each numeric-touching item (kernel merge, constant extraction, renames).
- **Robustness — new unit tests (RED first):**
  - parser rejects `atoms 0` / `idx < 1`; rejects sphere `radius <= 0`.
  - `build.rs` index lowering errors (not clamps) on invalid index.
  - `PyTarget::new` rejects `count == 0` and empty index lists.
  - malformed core result surfaces `PyValueError` (not `PanicException`) across
    `python/tests/test_*.py`.
- **Orchestration unit tests** (§6 of the review): phase transitions, the
  swap-state save/restore contract, handler callback ordering.
- **Bench gate:** `cargo bench --bench pack_end_to_end -- mixture` and the
  `pair_kernel` bench bracket the kernel merge.
- Consider `proptest` for the parser (stretch).

## Doc plan

- `docs/extending.md` — document the restraint registry + which kinds are
  reachable from `.inp`/Python (currently a 14-vs-5 mismatch); add the
  "write your own `Restraint`" example.
- README keyword table — mark restraints not yet wired into parser/Python
  until the registry closes the gap.
- `docs/packmol_parity.md` — cross-reference the `// Packmol:` rename comments.
- rustdoc — implementation docstrings on the 14 restraints (penalty formula +
  params), esp. cylinder/gaussian.
- `docs/architecture.md` — update once `PackContext` is decomposed.
- `CHANGELOG` — note the `Handler`-trait breaking change and the deprecated
  alias removals (`BBox`/`FixedPlacement`/`FromRegion`) scheduled for 0.2.

## Tiers (implementation order)

**Tier 1 — low risk, high value (do first)**
1. ✅ **DONE** FFI: `.expect()` → `PyValueError` in `positions`/`frame`/
   `elements` getters (`python/src/packer.rs`).
2. ✅ **DONE** Index bug: parser rejects `atoms 0` (`idx < 1`); `build.rs`
   `saturating_sub` documented as guarded by the parser invariant.
3. ✅ **DONE (partial)** Validation: sphere `radius > 0` (parser), `count > 0`
   (`python/src/target.rs`), non-empty index list (already present in parser).
   *Validation lives inline at the parse/decode boundary, not centralized in
   `validation.rs` — accepted for this slice (see Risks; tracked as debt).*
4. ✅ **DONE** De-Fortran naming with `// Packmol:` comments.
   - ✅ kernels `fparc/gparc/fgparc` → `pair_f_atom/pair_g_atom/pair_fg_atom`;
     `*_stats` → `*_parallel`.
   - ✅ euler `beta/gama/teta` → `euler_beta/euler_gamma/euler_theta` (only in
     euler-context files; CG/line-search `beta` correctly left alone).
   - ✅ `ilubar/ilugan` (+ `_start`/`_offset`/`_bad`/`_good` compounds) →
     `x_com_offset/x_euler_offset` family.
   - ✅ `fdist/frest` documented (public API — not renamed).
   - ✅ `ibtype/ibmol` → `atom_type_idx/atom_mol_idx`, `comptype` →
     `is_type_active`, `move_flag` → `selective_repack_mode`, with the
     mirror-sync methods (`set_ibtype`→`set_atom_type_idx`,
     `cached_comptype`→`cached_is_type_active`, …) co-renamed.
5. ✅ **DONE** Magic constants: `PENALTY_GRAD_COEFF` (kernel `4.0`, 6 sites),
   `MOVEBAD_FIMP_THRESHOLD` (movebad `10.0`), `FIMP_CLAMP` (`fimp` clamp). CG
   `GAMMA`/`THETA` were already named constants (`gencan/mod.rs:161-162`).
6. ✅ **DONE** `restraint.rs` (952 lines) split into `restraint/{mod,inside,
   outside,plane,cylinder,gaussian}.rs`; added a "write your own restraint"
   doctest + penalty-formula docstrings on the opaque kinds (cylinder/gaussian).

> **Progress (this slice):** Tier-1 #1–#3 landed + regression tests
> (`reject_zero_atom_index`, `reject_inside/outside_sphere_nonpositive_radius`,
> `test_zero_count_target_raises_value_error`). `parse_inside`/`parse_outside`/
> `parse_plane_above`/`parse_plane_below` extracted to
> `src/script/restraint_parse.rs` to hold `parser.rs` under the 800-line cap and
> pre-stage the §4 registry home. Verified: core suite, clippy (`--features
> cli`), Packmol `examples_batch`, and 120 pytest all green. CHANGELOG updated.
> Prereq: bumped molrs pin `0.0.15`→`0.0.16` (separate change) to unblock the
> local build. Remaining: Tier-1 #4–#6, then Tier 2/3.

**Tier 2 — medium risk (bench-gated)**
7. ✅ Centralized the duplicated short-radius penalty formula across the 5 pair
   kernels into `#[inline(always)]` helpers (parity bit-identical; full
   single-kernel merge not needed — the prologue was already shared via
   `AtomHotState`).
8. ✅ Split `pack()` (354→~180) into `validate_inputs` + `prepare_system`.
   **`tn_linesearch` left cohesive** — its extrapolation/interpolation loops
   share ~25 locals; extracting them would need 25-arg helpers, hurting
   readability. (Documented exception.)
9. ✅ `run_iteration` 19 params → `IterationConfig` + `IterationState`.
10. ✅ `evaluate_unscaled` radius swap → `UnscaledRadii` RAII guard.

**Tier 3 — architecture (medium-term)**
11. ✅ Split `PackContext` into `EvalState`/`Topology`/`CellGrid`/`PbcParams`.
12. ✅ Single restraint registry: all 14 kinds reachable from `.inp` + Python
    via one keyword→spec table + one `spec_to_restraint` map.
13. ✅ Dropped the `Constraints` ZST shell; typed the phase with `enum Phase`.
14. ✅ `Handler` I/O callbacks return `Result` (`PackError::HandlerIo`).

## Risks / open questions

- **File-size budget (HIGH, dominant).** Six touched files already exceed the
  800-line cap: `objective.rs` 1305, `packer.rs` 1069, `restraint.rs` 952,
  `pack_context.rs` 908, `parser.rs` 777 (near), `handler.rs` 378. Sequence
  *reducing* changes (kernel merge, `pack()`/`tn_linesearch` splits) before/with
  *adding* ones; promote `packer.rs`→`packer/` and `restraint.rs`→`restraint/`
  so new concerns get new files, not monolith growth. New file >800 = CRITICAL.
- **Layering (MEDIUM, one trap).** The restraint registry must live in
  `restraint.rs`/`region.rs` (lowest layer), never in `script/` — else
  `python/src/` gains an upward `script` dependency. The `Constraints` removal
  must not leave `script`/`objective` importing a deleted path.
- **Feature-gate hygiene (LOW).** Registry must not pull `io`-gated symbols into
  the default-feature path (Python builds without `io`); `validation.rs` must
  stay `io`-free. Re-run `--no-default-features` and
  `--no-default-features --features rayon` after registry + validation work.
- **Naming hygiene (LOW).** Renames use `// Packmol:` in comments (allowed);
  verify no rename introduces a `Packmol`-bearing *public* identifier. Registry
  keyword strings ("inside box", …) are data, not identifiers.
- **Immutability (LOW).** `IterationState` / `PackContext` sub-structs are
  internal mutable eval state, not public builders — fine. If the registry adds
  a public builder, it must use `self -> Self`, not `&mut self -> &mut Self`.
- **Open:** Does anything outside `src/` re-export `constraints`? Confirm before
  dropping the module (Tier 3 #13).
- **Open:** Is the `Handler`-trait breaking change acceptable pre-1.0, or should
  it ship behind a default method to stay source-compatible? (Tier 3 #14.)
- **Open:** `cargo build/clippy/test/bench` require the `../molrs` path dep —
  Tier 2/3 must be implemented where the workspace builds and bench-gated.
