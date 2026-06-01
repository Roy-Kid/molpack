---
name: mpk-fix
description: Minimal-diff bug fix — reproduce, write a regression test that fails, fix, verify. Writes code + tests, never refactors opportunistically.
argument-hint: "<bug description or failing test name>"
user-invocable: true
---

# mpk-fix — Minimal-Diff Bug Fix

Read CLAUDE.md for molpack conventions.

## Procedure

1. **Reproduce.** From `<arg>`, identify the failing scenario. If the user supplied a failing test, run it; otherwise locate or write the smallest reproduction.
2. **Regression test first.** Write the failing test in the right file (use `mpk-tester`'s test-layout knowledge). Confirm RED with the right reason.
3. **Diagnose** in place — read the suspect code path. Do not modify code yet.
4. **Smallest fix.** Edit only what is necessary to flip the test GREEN. Resist the urge to refactor adjacent code, rename variables, or "clean up" while you're there — open a separate `/mpk-refactor` if cleanup is warranted.
5. **Run impacted tier.** If the fix touches a hot-path file, run the relevant criterion bench. If it touches `restraint.rs`, `objective.rs`, `packer.rs`, `gencan/`, `initial.rs`, `relaxer.rs`, or `movebad.rs`, run `examples_batch`.
6. **Changelog.** Update `CHANGELOG.md` if the bug was user-visible.

## Output

`fix: {one-line cause} — {N} lines changed in {file}; regression test {test_name} now green`
