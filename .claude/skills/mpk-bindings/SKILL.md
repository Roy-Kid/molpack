---
name: mpk-bindings
description: Single-axis deep review of the Python wheel — naming, feature gate, error mapping, API mirror, test/doc/stub coverage. Read-only.
argument-hint: "<git ref or path scope>"
user-invocable: true
---

# mpk-bindings — Python Bindings Deep Dive

Read CLAUDE.md for molpack conventions.

## Procedure

1. Resolve scope from `<arg>` (default: `git diff --name-only` working tree). **Always** include both `src/` and `python/` in the file list passed to the agent — Rust changes can break the Python mirror even when no `python/` file is in the diff.
2. Delegate to the `mpk-bindings` agent with the combined file list.
3. Render the agent's findings as a severity table.

## Output

```
| Severity | File:line | Message |
|---|---|---|

verdict: {APPROVE | REQUEST CHANGES | BLOCK}
```
