---
name: mpk-test
description: Single-axis test review or test-add task — coverage gaps, missing gradient checks, regression-suite gating. Writes tests when asked.
argument-hint: "<git ref, path, or 'add: <description>'>"
user-invocable: true
---

# mpk-test — Test Coverage / TDD

Read CLAUDE.md for molpack conventions.

## Procedure

1. **Mode select:**
   - If `<arg>` starts with `add:`, treat as a TDD ask. Delegate to `mpk-tester` with the description; the agent writes the failing test in the right file (RED only — no implementation).
   - Otherwise (review mode): resolve scope from `<arg>` (default `git diff --name-only`). Delegate to `mpk-tester` for a coverage audit.
2. Render findings (review mode) or report which test was added and that it is RED (add mode).

## Output

```
| Severity | File:line | Message |
|---|---|---|
verdict: {APPROVE | REQUEST CHANGES | BLOCK}     ← review mode
```

```
added: tests/{file}::{test_name} (RED) — run /mpk-impl or /mpk-fix to make it green     ← add mode
```
