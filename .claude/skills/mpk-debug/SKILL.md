---
name: mpk-debug
description: Diagnose a failure or unexpected behavior in molpack. Read-only — never edits code, tests, or docs.
argument-hint: "<symptom — failing test, wrong output, perf regression>"
user-invocable: true
---

# mpk-debug — Diagnose Only

Read CLAUDE.md for molpack conventions.

**This skill never edits files.** Its only output is a diagnosis. If a fix is needed, the user runs `/mpk-fix` next.

## Procedure

1. **Restate the symptom** in one sentence.
2. **Reproduce read-only.** Run the failing test or offending command (`cargo test ...`, `cargo run --example ... --features io`, `pytest ...`). Capture actual output.
3. **Bisect by signal.**
   - *Numerical drift* (Packmol regression failed): inspect commits since the last green baseline; concentrate on `restraint.rs`, `objective.rs`, `gencan/`, `initial.rs`, `relaxer.rs`, `movebad.rs`.
   - *Crash / panic*: read the backtrace; locate the exact `file:line`.
   - *Perf regression*: run the relevant criterion bench; compare against `target/criterion/` baseline if present.
   - *Python test failure*: check both the PyO3 binding under `python/src/` and the underlying Rust API.
4. **Form a hypothesis** with `file:line` citations.
5. **Stop here.** Do not edit. Hand off to `/mpk-fix` if the user agrees with the hypothesis.

## Output

```
diagnosis: {symptom}
hypothesis: {one-paragraph cause, with file:line citations}
suggested next: /mpk-fix "{summary}"
```
