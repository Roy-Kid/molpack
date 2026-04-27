---
name: mpk-spec
description: Capture a feature spec into .claude/specs/{slug}.md and append to the spec INDEX. Writes specs only — does not implement code.
argument-hint: "<feature description>"
user-invocable: true
---

# mpk-spec — Feature Spec

Read CLAUDE.md for molpack conventions.

## Procedure

1. **Distill** `<arg>` into a slug (kebab-case, ≤4 words). If the user gave a long description, restate the intent in one sentence and confirm before continuing.
2. **Delegate** to the `mpk-architect` agent with the question: "Where in the existing module layering should this land, and what compliance risks does it carry?" — receive a placement note.
3. **Draft** `.claude/specs/{slug}.md` with these sections:
   - **Goal** — one paragraph; what observable behavior is added or changed.
   - **Non-goals** — what this spec deliberately does *not* do.
   - **Public surface** — Rust + Python + CLI changes; exact symbol / keyword / type names.
   - **Module placement** — from the architect note.
   - **Numerical contract** — gradient invariants, parity-suite expectations, optimizer impact (leave blank if N/A).
   - **Test plan** — which `tests/*.rs` and / or `python/tests/test_*.py` files; whether `examples_batch` must pass.
   - **Doc plan** — which entries in the documenter's sync table must follow.
   - **Risks / open questions**.
4. **Append** to `.claude/specs/INDEX.md`: `- [{slug}](./{slug}.md) — {one-line goal} — DRAFT`.
5. Show the spec path; ask the user to review before `/mpk-impl {slug}`.

## Output

`spec: .claude/specs/{slug}.md (DRAFT) — review then run /mpk-impl {slug}`
