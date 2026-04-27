---
name: mpk-arch
description: Single-axis deep architecture review — layering, feature gates, naming hygiene, immutability, file-size budget. Read-only.
argument-hint: "<git ref or path scope>"
user-invocable: true
---

# mpk-arch — Architecture Deep Dive

Read CLAUDE.md for molpack conventions.

## Procedure

1. Resolve scope from `<arg>` (default: `git diff --name-only` working tree).
2. Delegate to the `mpk-architect` agent with the file list.
3. Render the agent's findings as a severity table.

## Output

```
| Severity | File:line | Message |
|---|---|---|

verdict: {APPROVE | REQUEST CHANGES | BLOCK}
```
