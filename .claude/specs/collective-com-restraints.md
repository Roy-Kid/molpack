---
title: collective-com-restraints — per-molecule collective coordinates and a reciprocal-space density-fluctuation restraint
status: draft
created: 2026-07-25
---

# collective-com-restraints — packing to collective targets

## Summary

Extend the collective-restraint family from per-atom coordinates to **per-copy
(per-molecule) collective coordinates**, and add a **reciprocal-space
density-fluctuation restraint** on molecule centres.

Three additions, all rigid-body and all force-field-free:

1. a `com` reduction that turns a species' atoms into one site per copy, with
   the gradient scattered back over that copy's atoms;
2. a `StructureFactor` restraint that drives the centre-of-mass structure factor
   `S(q)` towards a target over a shell of reciprocal-lattice vectors;
3. an `axis` geometry giving a per-copy orientation coordinate, so an
   orientation *distribution* can be a packing target.

Together these make the polymer community's "prepacking" step a first-class
constraint inside the packer, and generalise it from *uniform density* to *any
target profile* plus orientation.

> **Partly unblocked (2026-09-08).** Two prerequisites this spec assumed it would
> have to build have landed with `SelfSeparation` (see *Relation to
> `SelfSeparation`* below):
>
> - **The collective seam now carries a context.** `Restraint::f`/`fg` take a
>   `GroupEvaluation { scale, scale2, natoms_per_copy, mic }` instead of the two bare
>   scales. Task 2's wavevector enumeration needs the cell and Task 1's
>   reduction needs the per-copy width; both now arrive at evaluation time
>   rather than being cached at construction — which matters, because the cell
>   is resolved *after* the targets are lowered, so a restraint that stored a
>   box at construction could store the wrong one. The `Mic` on `GroupEvaluation` is
>   the pair loop's own; a reciprocal lattice for `StructureFactor` should be
>   added to `GroupEvaluation` the same way rather than re-derived.
> - **The forward/backward COM math exists**, in
>   `src/restraint/collective/com.rs` (`centroids` / `scatter`, geometric
>   weights, unit-tested). Task 1's `SiteReduction` should **wrap** it, not
>   restate it — `ComWeights::Mass` is the part still to build.
> - **`GroupEvaluation` carries the cell**, not just the minimum image, precisely so a
>   group-level term can *partition* space. `SelfSeparation` bins its centres
>   into a `molrs` `CellGrid` at `d_min` and sweeps the forward stencil, the
>   same primitive the pair loop uses. `StructureFactor` needs the same box for
>   its reciprocal lattice and should take it from there rather than adding a
>   second channel.

## Domain basis

- In a melt, excluded volume is screened (Flory ideality), so single-chain
  statistics are known independently of packing. Auhl, Everaers, Grest, Kremer &
  Plimpton (J. Chem. Phys. **119**, 12718, 2003) exploit exactly this scale
  separation: generate chains as random walks with the correct intramolecular
  statistics ignoring inter-chain overlap, then run a **prepacking Monte Carlo
  that translates and rotates whole chains, with intramolecular structure held
  fixed, to minimise density fluctuations**, and only then switch on excluded
  volume through a soft push-off.
- The prepacking step's degrees of freedom are therefore exactly a packer's
  degrees of freedom (rigid-body centre + orientation), and its objective is a
  collective functional over all molecules — not a pairwise distance. This is
  the argument for putting it in molpack rather than in a chain generator.
- What it buys: independently placed chains leave density fluctuations
  correlated on the scale of Rg. Once excluded volume is switched on, those
  relax only by centre-of-mass motion over the chain's own size, i.e. on the
  disentanglement time τ_d ∝ N^3.4, with total equilibration effort scaling as
  ≈ N^4.9. The collective term is paid once, at packing time.
- Long-wavelength density fluctuations are the small-|q| limit of the
  centre-of-mass structure factor, so `S(q)` on a shell of small reciprocal
  vectors is the natural, grid-free, analytically differentiable form of the
  same objective. For a periodic cell the admissible wavevectors are the
  reciprocal-lattice vectors `q = 2π (h⁻¹)ᵀ n`, `n ∈ ℤ³`.
- Packmol has no counterpart. Its objective terms are all local and hard: an
  atom pair, or an atom against a geometric primitive. `xygauss` is the closest
  thing and it is a hard-coded Gaussian *surface bound* restricted to the xy
  plane, not a distribution target.

## Design

### 1. `com` reduction

`src/restraint/collective/` currently maps a species' **atom** coordinates to a
scalar reaction coordinate ξ through a *geometry* (`plane`, `point`) and matches
the empirical distribution of ξ to a target through the squared 1-D Wasserstein
(sorted-CDF) engine.

Add a reduction layer in front of the geometry:

```rust
pub enum SiteReduction {
    /// One site per atom — today's behaviour, the identity reduction.
    Atoms,
    /// One site per copy at its geometric or mass-weighted centre.
    Com { weights: ComWeights },
}
```

Forward: for copy `c` with atoms `i ∈ c` and weights `w_i` (normalised within
the copy), `R_c = Σ_i w_i r_i`. Backward: `∂L/∂r_i += w_i · ∂L/∂R_c`. The
existing geometries and distributions are unchanged and compose on top, so
`Com + plane + tabulated` immediately gives "the centres of this species follow
this density profile along z" — the polymer form of the existing profile
restraint.

`ComWeights::Uniform` (geometric centre) is the default; `ComWeights::Mass`
requires masses on the target and errors at build time if absent.

### 2. `StructureFactor` restraint

New `src/restraint/collective/structure_factor.rs`. Over the reduced sites
`{R_c}`, `c = 1..N`:

```
A_q = Σ_c cos(q·R_c)          B_q = Σ_c sin(q·R_c)
S(q) = (A_q² + B_q²) / N
E    = Σ_{q ∈ Q} w_q · (S(q) − S*(q))²
```

with the analytic gradient

```
∂S(q)/∂R_c = (2/N) · [ B_q cos(q·R_c) − A_q sin(q·R_c) ] · q
∂E/∂R_c    = Σ_q 2 w_q (S(q) − S*(q)) · ∂S(q)/∂R_c
```

then scattered onto atoms by the `com` reduction. Cost is `O(N · |Q|)` with no
grid and no FFT; `|Q|` is a few hundred at most.

**Wavevector set.** `Q` is generated from the `SimBox` reciprocal lattice
`q = 2π (h⁻¹)ᵀ n` for `n ∈ ℤ³ \ {0}` with `|q| ≤ q_max`, deduplicated by the
`q ↔ −q` symmetry. Only lattice-commensurate wavevectors are admissible under
PBC, so this set is not a modelling choice — it is the complete admissible set
below the cutoff. Non-periodic axes contribute no reciprocal direction and are
excluded from the enumeration.

**Target.** `S*(q)` is a constant (default) or a tabulated curve of `|q|`.
`S* = 0` is the maximally-suppressed limit; an incompressible melt has a small
non-zero plateau, so the constant is user-supplied rather than hard-coded to
zero. A configuration of independently placed centres has `S(q) ≈ 1`, which is
the baseline the figure compares against.

### 3. `axis` geometry

New geometry in `src/restraint/collective/geometry.rs`: the per-copy unit vector
`u_c` defined by two named atom indices within the template
(`u_c = (r_b − r_a)/|r_b − r_a|`), with reaction coordinate `ξ_c = u_c · n̂`
for a user-given reference direction `n̂` (optionally `|u_c · n̂|` when the
molecular vector is apolar). The gradient is the standard normalisation
derivative scattered onto atoms `a` and `b` only.

Composed with the existing distributions this expresses "the chain vectors of
this species follow this orientation distribution" — e.g. a target nematic order
at an interface. Packmol's `constrain_rotations` can only impose a hard angular
bound.

### 4. Surface

Rust builder: `Target::with_collective_restraint` already exists and takes the
new restraint types unchanged. Python: expose `ComReduction`, `StructureFactor`
and `AxisGeometry` alongside the existing collective restraints. Script: a
`structure_factor` keyword and an `orient` keyword mirroring the existing
`profile` keyword.

## Files

- `molpack/src/restraint/collective/reduction.rs` — new, `SiteReduction`
- `molpack/src/restraint/collective/structure_factor.rs` — new
- `molpack/src/restraint/collective/geometry.rs` — `axis` geometry
- `molpack/src/restraint/collective/mod.rs` — wiring and re-exports
- `molpack/src/target.rs` — reduction attached to a collective binding
- `molpack/src/script/…` — `structure_factor` / `orient` keywords
- `molpack-python` bindings
- in-module `#[cfg(test)]` cases in the owning `src/restraint/collective/*.rs` files, plus
  through-the-objective parity in `src/restraint/geometric/tests/gradient.rs` (the former integration files were deleted 2026-09-20)

## Tasks

1. `SiteReduction` + gradient scatter; `Atoms` path provably identical to today.
2. `StructureFactor` restraint: reciprocal-lattice enumeration, energy, analytic
   gradient, finite-difference test.
3. `axis` geometry + gradient.
4. Script + Python surface.
5. End-to-end cases: uniform melt, tabulated lamellar profile, orientation
   target.

## Testing

**Gradients.** Every new term gets a central finite-difference check in
the owning collective module's tests (pattern: `gaussian.rs::plane_gradient_matches_finite_difference`)
and through the objective in `src/restraint/geometric/tests/gradient.rs`
(pattern: `collective_restraint_gradient_matches_finite_difference_through_the_objective`) at the existing tolerance, on both an orthorhombic and a
triclinic cell (the reciprocal lattice of a tilted cell is the case most likely
to be wrong).

**Uniform melt (prepacking reproduction).** Pack the same conformer pool at the
same density and tolerance twice — once with the pairwise objective alone, once
with `StructureFactor` added — and compare `S(q)` on the lowest shells. The
collective run must suppress the small-|q| plateau by at least an order of
magnitude while still satisfying the packing tolerance (`State::frest == 0` and
`State::fdist <= precision` in both runs; the validation report no longer exists). *Identical conformer pool in both arms* — the single-variable
comparison is the whole point, and it is also how the paper figure is built
against Packmol.

**Target profile.** With a tabulated lamellar profile as target, the packed
centre-of-mass profile along the lamellar normal matches the target within a
stated W₂ tolerance, and the pair-distance constraint is still satisfied.

**Orientation.** With a target distribution over `cos θ`, the packed
distribution matches within a stated W₂ tolerance; the resulting order parameter
is reported.

**No harm.** All existing collective-restraint tests pass unchanged; a packing
run with no collective restraint is byte-identical to before.

## Out of scope

- **Internal degrees of freedom.** Molecules stay rigid; per-copy Rg / Ree are
  fixed by the input conformer and are not optimisation variables. Matching an
  intramolecular size distribution is a **conformer-pool** question, solved in
  closed form by sorted quantile assignment (1-D W₂), and belongs to input
  preparation, not to the packer. One sentence in Methods, no code here.
- Chain generation of any kind.
- Push-off / soft-core ramps / MD equilibration — those belong to the MD engine.
- Grid- or FFT-based density functionals; the reciprocal-space form is the
  scope.
- Force fields. Nothing here requires the `ff` feature.
- Any performance claim.

## Relation to `SelfSeparation`

`SelfSeparation` (landed 2026-09-08, `src/restraint/collective/separation.rs`)
is the real-space, pairwise member of the collective family: a lower bound on
the centre-to-centre distance between two copies of one species. It is **not**
this spec's `StructureFactor` in another basis, and the two are not parallel
abstractions for one concept — they measure different quantities and compose:

| | `SelfSeparation` | `StructureFactor` (this spec) |
|---|---|---|
| space | real, pairwise | reciprocal, collective |
| constraint | **local**: no two copies closer than `d_min` | **global**: suppress the small-\|q\| density-fluctuation plateau |
| satisfied when | a lower bound holds | the empirical `S(q)` approaches `S*(q)` |
| says nothing about | how copies are distributed above `d_min` | how close any particular pair may come |
| dependencies | none | `triclinic-cell-downshift` |

A melt run may want both: the separation bound keeps chains off each other
locally, `S(q)` removes the Rg-scale density correlations. Neither subsumes the
other, and neither should be reimplemented in the other's terms.

## Dependencies

`triclinic-cell-downshift` — for the reciprocal lattice of a general cell.
The orthorhombic path can be developed and tested before it lands.
