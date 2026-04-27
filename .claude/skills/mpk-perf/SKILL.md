---
name: mpk-perf
description: Single-axis deep perf review — hot-path allocations, criterion-bench coverage, rayon equivalence, release-profile guards. Read-only.
argument-hint: "<git ref or path scope>"
user-invocable: true
---

# mpk-perf — Performance Deep Dive

Read CLAUDE.md for molpack conventions.

## Procedure

1. Resolve scope from `<arg>` (default: `git diff --name-only` working tree).
2. Delegate to the `mpk-perf` agent with the file list. Ask it to note whether the corresponding criterion benches were re-run, and (for rayon-touching code) whether `tests/parallel_equivalence.rs` was re-run.
3. Render the agent's findings as a severity table.

## Output

```
| Severity | File:line | Message |
|---|---|---|

verdict: {APPROVE | REQUEST CHANGES | BLOCK} — benches: {re-run ✓ | not re-run | N/A} — rayon equivalence: {re-run ✓ | N/A}
```
