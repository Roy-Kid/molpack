---
name: mpk-note
description: Capture a project decision into .claude/NOTES.md, detect conflicts with existing notes and CLAUDE.md, sweep stale notes, and promote stable ones into CLAUDE.md.
argument-hint: "<the decision or rule to remember>"
user-invocable: true
---

# mpk-note — Project Memory

Read CLAUDE.md for molpack conventions and `.claude/NOTES.md` for prior notes.

## Procedure

1. **Capture.** Restate `<arg>` as a one-paragraph note: rule + **Why:** + **How to apply:**.
2. **Conflict scan.** Grep CLAUDE.md and existing NOTES.md entries for overlap. If a conflict exists:
   - report it with `file:line` citations
   - ask the user: "supersede the old note, narrow the scope, or drop the new one?"
   Do not write until the user resolves.
3. **Stale sweep.** Re-read NOTES.md. For each entry: does the cited code path still exist? Does the rule still match observed code? Mark stale entries for the user's review (do not delete unilaterally).
4. **Promote.** If a NOTES.md entry has been stable across ≥3 invocations of `/mpk-note`, or the user says "this is a hard rule", move it into CLAUDE.md (under the appropriate section) and DELETE the NOTES.md entry. Never leave the same rule in both files.
5. **Append** the new note to NOTES.md with today's date.

## Output

`note added: {one-line title} — wrote to {NOTES.md | CLAUDE.md} — conflicts: {N} — stale flagged: {N}`
