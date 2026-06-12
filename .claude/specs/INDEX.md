# .claude/specs/INDEX.md — molpack feature specs

Status legend: **DRAFT** (under review) → **IN PROGRESS** (being implemented) → **DONE**.

Add via `/mpk-spec <feature description>`. Implement via `/mpk-impl <slug>`.

- [gencan-anneal-budget](./gencan-anneal-budget.md) — pack_solvprotein 5.5× slowdown: real root cause was a missing `avoid_overlap` fixed-atom rejection in `initial.rs` (solvent seeded inside the fixed protein); fixed → now faster than packmol. Proposed anneal-budget workaround SUPERSEDED — RESOLVED

## profile-distribution-restraints (chain — implement in order, start with `-01-split`)

Soft restraints that bias selected sites toward a target spatial distribution ρ*(ξ) via Boltzmann inversion, composing a geometry coordinate (planar/radial/cylindrical/region-distance) × a distribution (Gaussian/erf/tanh/exponential/tabulated). Split into 6 sub-specs by the large-spec rule.

- [profile-distribution-restraints-01-split](./profile-distribution-restraints-01-split.md) — refactor `src/restraint.rs` (976 LOC) into a `src/restraint/` module; no behaviour change [done]
- [profile-distribution-restraints-02-coordinate](./profile-distribution-restraints-02-coordinate.md) — `coordinate.rs`: the 4 reaction coordinates ξ(x) + analytic ∇ξ + r_guard + minimum-image [done]
- [profile-distribution-restraints-03-distribution](./profile-distribution-restraints-03-distribution.md) — `distribution.rs`: closed-form Boltzmann-inversion penalties + density floor + shell-Jacobian (§6.2) [done]
- [profile-distribution-restraints-04-spline](./profile-distribution-restraints-04-spline.md) — `spline.rs`: C¹ monotone cubic + numerical inversion for the tabulated profile [done]
- [profile-distribution-restraints-05-compose](./profile-distribution-restraints-05-compose.md) — `ProfileRestraint`: compose coordinate × distribution and impl `Restraint` (chain rule) [done]
- [profile-distribution-restraints-06-script](./profile-distribution-restraints-06-script.md) — `.inp` `profile` keyword + parser + lowering + end-to-end pack tests (criteria 1–4) [done]
