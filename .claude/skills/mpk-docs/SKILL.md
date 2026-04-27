---
name: mpk-docs
description: Audit or update docs — rustdoc, docs/, python/docs/, README, CONTRIBUTING, CHANGELOG. Writes docs when asked.
argument-hint: "<git ref, path, or 'update: <topic>'>"
user-invocable: true
---

# mpk-docs — Documentation

Read CLAUDE.md for molpack conventions.

## Procedure

1. **Mode select:**
   - If `<arg>` starts with `update:`, delegate to `mpk-documenter` with the topic; it edits the relevant docs per its sync table.
   - Otherwise (audit mode): resolve scope from `<arg>` (default `git diff --name-only`). Delegate to `mpk-documenter` for a doc-sync audit.
2. Render findings (audit mode) or list edited files (update mode).

## Output

```
| Severity | File:line | Message |
|---|---|---|
verdict: {APPROVE | REQUEST CHANGES | BLOCK}     ← audit mode
```

```
updated: {file}, {file}, ...     ← update mode
```
