---
name: mpk-numerics
description: Single-axis deep review of gradient correctness, restraint math, optimizer behavior, and Packmol parity. Read-only.
argument-hint: "<git ref or path scope>"
user-invocable: true
---

# mpk-numerics — Numerics Deep Dive

Read CLAUDE.md for molpack conventions.

## Procedure

1. Resolve scope from `<arg>` (default: `git diff --name-only` working tree).
2. Delegate to the `mpk-numerics` agent with the file list. Ask it to confirm whether the Packmol regression suite (`tests/examples_batch.rs`) was rerun for the diff.
3. Render the agent's findings as a severity table.

## Output

```
| Severity | File:line | Message |
|---|---|---|

verdict: {APPROVE | REQUEST CHANGES | BLOCK} — packmol parity: {ran ✓ | not run | N/A}
```
