---
name: mpk-refactor
description: Restructure code preserving observable behavior — same tests, same numerics, same public surface. Writes Rust source only.
argument-hint: "<what to restructure and why>"
user-invocable: true
---

# mpk-refactor — Behavior-Preserving Refactor

Read CLAUDE.md for molpack conventions.

## Procedure

1. **State the invariant.** Write down (in the user-facing response) what must NOT change: public API names, `.inp` keyword surface, gradient values, packed-output bit-equivalence under fixed seed. Confirm with the user before proceeding.
2. **Architecture pre-check.** Delegate to `mpk-architect` with the planned restructure (split file, move module, change feature gating). Stop on CRITICAL or HIGH.
3. **Pin behavior.** Delegate to `mpk-tester` to ensure the affected behavior is covered. Add coverage if thin **before** refactoring — never after.
4. **Refactor.** Move / split files. Respect the file-size budget (target 200–400 lines, hard cap 800).
5. **Architecture post-check.** Delegate to `mpk-architect` on the diff.
6. **Verify behavior unchanged.** Run the fast suite plus Packmol regression if the refactor touched any numerical file. The regression suite output must remain byte-identical.
7. **Bindings + docs.** If a public symbol moved, delegate to `mpk-bindings` (re-export paths, mirror) and `mpk-documenter` (rustdoc cross-links).

## Output

`refactor: {summary} — {N} files moved/split; behavior pinned by {M} tests; arch ✓ / parity ✓`
