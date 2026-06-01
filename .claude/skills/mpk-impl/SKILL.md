---
name: mpk-impl
description: Implement a spec end-to-end via TDD — RED first, then GREEN, then verify across architecture, numerics, bindings, and docs. Writes code, tests, and docs.
argument-hint: "<spec slug or short description>"
user-invocable: true
---

# mpk-impl — TDD Implementation

Read CLAUDE.md for molpack conventions.

## Procedure

1. **Locate the spec.** If `<arg>` matches a slug in `.claude/specs/INDEX.md`, read that spec. Otherwise treat `<arg>` as an inline mini-spec; surface the gap to the user before proceeding.
2. **Scope-classify:**
   - *small* — ≤2 files in one module
   - *medium* — touches one subsystem (script, restraint, packer, gencan, bindings)
   - *large* — crosses the lifecycle boundary (script → packer → bindings)
3. **Architecture pre-check.** Delegate to `mpk-architect` with the planned file list. If it returns CRITICAL or HIGH, stop and revise the plan with the user.
4. **RED.** Delegate to `mpk-tester` to write failing tests in the right files. Confirm the tests fail with the expected reason.
5. **GREEN.** Implement the minimum code to pass. Respect immutability, file-size budget, naming hygiene (no "packmol" in public identifiers), and feature-gate boundaries.
6. **Numerics check.** For changes in `restraint.rs`, `objective.rs`, `packer.rs`, `gencan/`, `initial.rs`, `relaxer.rs`, `movebad.rs`, or `script/`: delegate to `mpk-numerics`.
7. **Architecture post-check.** Delegate to `mpk-architect` again on the actual diff.
8. **Bindings mirror.** If public Rust API changed, delegate to `mpk-bindings`.
9. **Tests verify.** Re-delegate to `mpk-tester` to run the fast suite plus the Packmol regression where applicable; confirm GREEN.
10. **Docs.** Delegate to `mpk-documenter` to update rustdoc, `docs/`, `python/docs/`, README, CONTRIBUTING per the sync table.
11. **Mark spec DONE** in `.claude/specs/INDEX.md`.

## Output

`impl: {slug} — {N} files, {M} tests added; arch ✓ / numerics ✓ / bindings ✓ / docs ✓`
