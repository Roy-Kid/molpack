---
name: mpk-review
description: Multi-axis review of the current diff — fans out to architect, numerics, perf, bindings, tester, documenter in parallel and aggregates findings.
argument-hint: "<git ref or 'HEAD' for working tree>"
user-invocable: true
---

# mpk-review — Aggregate Multi-Axis Review

Read CLAUDE.md for molpack conventions.

## Procedure

1. **Scope.** From `<arg>` (default working tree), determine the diff range:
   - working tree: `git diff --name-only`
   - committed: `git diff <arg>^..<arg> --name-only`
2. **Fan out in parallel** — issue all six delegate calls in a single message:
   - `mpk-architect` — layering, feature gates, naming, immutability
   - `mpk-numerics` — gradients, Packmol parity, optimizer
   - `mpk-perf` — hot-path allocations, benches, rayon equivalence
   - `mpk-bindings` — Python mirror, naming, error mapping, docs
   - `mpk-tester` — coverage gaps, missing gradient checks, regression-suite gating
   - `mpk-documenter` — doc-sync gaps per the sync table
3. **Aggregate** all `[SEVERITY] file:line — message` findings into one table.
4. **Verdict:**
   - any CRITICAL → **BLOCK**
   - any HIGH → **REQUEST CHANGES**
   - else → **APPROVE**

## Output

```
| Severity | File:line | Agent | Message |
|---|---|---|---|

verdict: {APPROVE | REQUEST CHANGES | BLOCK} — {C} critical, {H} high, {M} medium, {L} low
```
