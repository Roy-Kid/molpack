# .claude/specs/INDEX.md — molpack feature specs

Status legend: **DRAFT** (under review) → **IN PROGRESS** (being implemented) → **DONE**.

Add via `/mpk-spec <feature description>`. Implement via `/mpk-impl <slug>`.

- [gencan-anneal-budget](./gencan-anneal-budget.md) — pack_solvprotein 5.5× slowdown: real root cause was a missing `avoid_overlap` fixed-atom rejection in `initial.rs` (solvent seeded inside the fixed protein); fixed → now faster than packmol. Proposed anneal-budget workaround SUPERSEDED — RESOLVED
