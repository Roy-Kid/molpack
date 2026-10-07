# Software Engineering Laws

Every rule here outranks scope, minimal-diff, and convenience. There is
no "just this once". CLAUDE.md indexes one line per law.

Adding, changing, or repealing a law is the operator's act via
`/mol:note`. No skill retires a law on its own judgment.

## 0. Purpose

This file is not a style guide and not a pattern catalog.

Laws define **non-negotiable design constraints**. They constrain
architecture, ownership, dependency, state, and change. Patterns and
implementation techniques are subordinate to these laws. A law must be
strong enough to reject a concrete design in review.

If a rule cannot identify a concrete forbidden design, it is guidance,
not law.

A carve-out exists only where this file records it under § VII, naming
the subsystem. An agent never grants itself one.

---

# I. System shape

How the system as a whole should look.

<!-- mol:law:id:conceptual-integrity -->
## 1. Conceptual integrity

**Principle.** One problem should have one coherent conceptual model.

**Intent.** The system uses one set of concepts, terms, and abstractions.
A subsystem does not invent a sibling model of the same idea.

**Never**

- Never create parallel abstractions for the same concept.
- Never introduce aliases that develop independent semantics.
- Never solve local inconvenience by inventing a new conceptual layer.

**Derived guidance.** Prefer extending an existing concept over a sibling
concept. Shared vocabulary is part of architecture.

<!-- mol:law:id:architecture-first -->
## 2. Architecture first

**Principle.** Preserve a simple system shape before optimizing local
convenience.

**Intent.** Local coding convenience must not buy itself by breaking the
whole architecture. Simple, clear shape is the foundation of
maintainability and performance.

**Never**

- Never add a layer only because it may be useful later.
- Never introduce infrastructure for unmeasured performance concerns.
- Never let a local feature dictate global architecture.

**Derived guidance.** Prefer fewer architectural concepts. Prefer
removing indirection over explaining it.

<!-- mol:law:id:earn-complexity -->
## 3. Earn complexity

**Principle.** Every unit of complexity must be justified by demonstrated
pressure.

**Intent.** Complexity is not free. Need, performance, compatibility, or
extension must already exist before an abstraction does. Architecture
first governs shape; this law governs the complexity budget.

**Never**

- Never generalize for hypothetical future requirements.
- Never optimize without evidence.
- Never make something configurable merely because it could vary.
- Never add extensibility without an actual extension point.

---

# II. Boundaries and ownership

How the system is cut.

<!-- mol:law:id:locality-of-change -->
## 4. Locality of change

**Principle.** A local requirement should require a local change.

**Intent.** A good module boundary shows up as change locality, not as
an abstract cohesion score. High cohesion and low coupling are
consequences of this law.

**Never**

- Never require unrelated modules to change in lockstep.
- Never create dependency cycles.
- Never spread one responsibility across multiple owners.
- Never make callers understand unrelated subsystem details.

**Derived guidance.** A unit is green via `$META.build.test_single` on
its mirrored tests with fakes for outbound deps. If proving the unit
requires the full graph, the boundary is wrong — split, inject, or
`/mol:refactor`. Do not compensate with more integration tests.

<!-- mol:law:id:hide-decisions -->
## 5. Hide decisions, expose contracts

**Principle.** Implementation decisions stay behind their owning
boundary.

**Intent.** A module hides **decisions that may change**, not merely
lines in a different file.

**Never**

- Never leak representation details across module boundaries.
- Never expose internal lifecycle or storage decisions as public
  contract.
- Never require callers to reproduce internal policy.

**Derived guidance.** Program against stable contracts. An
implementation detail should be replaceable without rewriting
consumers.

<!-- mol:law:id:dependencies-follow-policy -->
## 6. Dependencies follow policy

**Principle.** Replaceable mechanisms depend on stable policy, never
the reverse.

**Intent.** Core semantics are not defined by UI, binding, framework,
storage, or transport.

    mechanism → policy

not

    policy → mechanism

**Never**

- Never make domain/core depend on UI.
- Never make core depend on a serialization format.
- Never make core depend on Python / Rust / WASM binding concerns.
- Never let a framework define domain semantics.

---

# III. Public surface

How others use the system.

<!-- mol:law:id:primitive-surface -->
## 7. Primitive public surface

**Principle.** Public APIs expose orthogonal primitives; composition
belongs to callers.

**Intent.** The API ships building blocks, not a hidden workflow.

**Never**

- Never provide an all-in-one façade for unrelated operations.
- Never make one public method perform several independently
  meaningful steps.
- Never encode one preferred workflow as the only API.
- Never duplicate primitives with convenience aliases that become
  separate contracts.

**Derived guidance.** High-level workflows may live outside the
primitive core (docs, examples, caller code).

<!-- mol:law:id:explicit-flow -->
## 8. Explicit flow

**Principle.** State transitions, ownership, and required ordering
must be explicit and enforceable.

**Intent.** A user must not enter an illegal state by forgetting a
step. This covers initialization, validation, lifecycle, context, and
state machines.

**Never**

- Never rely on hidden ambient context.
- Never expose `validate()` / `init()` steps callers can forget.
- Never depend on undocumented call ordering.
- Never encode required state in conventions alone.
- Never make illegal states trivially representable when the
  type/model can prevent them.

---

# IV. Truth and state

Who the system believes.

<!-- mol:law:id:one-home -->
## 9. One home per fact

**Principle.** Every authoritative fact has exactly one owner.

**Intent.** Avoid synchronization and drift. **Representation may be
duplicated; authority cannot.** A serialization copy may exist; it
must not become a second mutable truth.

**Never**

- Never maintain two independently mutable representations of the
  same truth.
- Never cache authoritative state without explicit invalidation
  semantics.
- Never copy configuration into another source of truth.
- Never infer and persist information that can be derived cheaply
  from its owner.

---

# V. Evolution

How the system changes without rotting.

<!-- mol:law:id:no-silent-debt -->
## 10. No silent debt

**Principle.** Debt must be explicit, bounded, and owned.

**Intent.** The worst debt is not a hack — it is a hack packaged as
normal architecture. A conscious exception is debt. An invisible
exception becomes architecture.

**Never**

- Never hide an architectural compromise inside an unrelated change.
- Never introduce temporary duplication without marking its removal
  path.
- Never normalize a workaround by silently building on top of it.
- Never leave known invariant violations undocumented.
- Never ignore, skip-mark, or weaken an assert on rot you already
  saw. Fix it if local and stage-allowed; else stop, report
  path:line, route `/mol:debug` / `/mol:refactor` / supersede, and
  name it in the summary.

Outranks "stay in scope" and "minimal diff".

---

# VI. Verification

How we prove the design has not decayed. Separate from architecture
laws.

<!-- mol:law:id:tests-owned-behavior -->
## 11. Tests verify owned behavior

**Principle.** Tests belong to the owner of the behavior they verify.

**Intent.** Tests verify a module's own contract, not the choreography
of the whole system.

**Never**

- Never test implementation details as public behavior.
- Never require unrelated subsystems merely to verify local
  semantics.
- Never use broad integration setup where a unit boundary is
  sufficient.

### Project testing policy: unit-only by default

New behavior must be unit-testable at its ownership boundary
(a `#[cfg(test)]` module next to the code, one module,
`$META.build.test_single`, fakes for outbound deps). molpack has **no**
integration, end-to-end or regression harness at all: a scenario that
can only be checked by running a whole pack is not a test here — it is a
runnable example, or evidence of a missing boundary.

Layout details: `tester` agent.

---

# VII. Exceptions

Any design that violates a law is recorded, not inferred:

    Law violated:
    Reason:
    Evidence:
    Scope:
    Removal condition:
    Owner:

Convenience is not sufficient justification. The exception is itself
an architecture decision. An agent never grants one from the task
text.

---

# VIII. Derived principles

These are consequences or heuristics, not laws. SOLID, YAGNI, DRY,
and the rest do not outrank this file. If someone cites them, first
show which law they serve.

| Heuristic | Comes from |
|---|---|
| YAGNI | Earn complexity |
| High cohesion / low coupling | Locality of change |
| Dependency inversion | Dependencies follow policy |
| Information hiding | Hide decisions, expose contracts |
| DRY (authority only) | One home per fact |
| Deep modules | Hide decisions + Primitive public surface |
| Composition over inheritance | Locality of change (a common means) |
| KISS / fewer boxes | Architecture first + Earn complexity |
| Program to interfaces | Hide decisions, expose contracts |
| Single responsibility / SoC | Locality of change + Primitive public surface |

---

# IX. Project invariants (molpack)

Moved here from CLAUDE.md "Hard rules" and from the 2026-08-28 solver
principles in `notes.md` (user ruling: "这个必须牢记"). Same template,
same protection. Where a rule maps onto a general law, the mapping is
named so reviewers can cite either.

<!-- mol:law:id:no-packmol-identifiers -->
## P1. No "packmol" in public identifiers

**Principle.** The product is `molpack`. Packmol is the compatible
script format and the reference implementation, not the code.

**Intent.** Public naming states what this project is. Prose (doc
comments, docs site) may cite Packmol freely; identifiers may not.

**Never**

- Never name a public symbol, module, type, feature, crate, or Python
  class with "packmol".
- Never let a ported Fortran routine's name become part of the public
  surface (keep it in a comment for traceability if needed).

**Derived guidance.** Conceptual integrity (§ 1): one product name,
one vocabulary.

<!-- mol:law:id:inp-only-config -->
## P2. Configuration is Packmol `.inp` only

**Principle.** The script format is the one configuration surface.

**Intent.** Users bring existing `.inp` files; a second configuration
language would fork the user base and the parser.

**Never**

- Never invent a TOML / YAML / JSON configuration surface.
- Never add a second way to express something the `.inp` grammar
  already expresses.
- Never introduce a growth or engine keyword per-structure; if the
  grammar is ever extended, it is a script-level entry-selection
  keyword (engine-entry-split ruling).

<!-- mol:law:id:molrs-pins-manual -->
## P3. molrs path and version pins are managed manually

**Principle.** The `../molrs/molrs` path dependency and its version
line are edited by a human, deliberately.

**Never**

- Never automate the pin check in pre-commit hooks or CI.
- Never let a hook or script rewrite `Cargo.toml` / `pyproject.toml`
  pins.

<!-- mol:law:id:local-gates-prek-tox -->
## P4. Local gates are prek + tox

**Principle.** Hooks use prek (pre-commit-compatible config); Python
isolation is tox from the `python/` `dev` dependency group.

**Never**

- Never add a project `scripts/` test wrapper.
- Never hand-write a local hook where a registry-hosted one exists
  (`doublify/pre-commit-rust`, `astral-sh/ruff`).
- Never run Python binding tests any way other than
  `uv run --directory python --group dev tox -e py` (non-editable,
  sibling molrs path install + maturin wheel + pytest).

<!-- mol:law:id:fork-pr-workflow -->
## P5. Fork → PR

**Principle.** `origin` is the Roy-Kid fork, `upstream` is
`MolCrafts/molpack`; changes land through pull requests.

**Never**

- Never push directly to `MolCrafts/molpack` master.
- Never force-push or rebase a branch on `upstream`.

<!-- mol:law:id:pure-geometry-solvers -->
## P6. Solvers are pure geometry

**Principle.** Everything on the packing-solver seam consumes only
topology (bond graph) and geometry (coordinates, radii, priors,
restraints). No force field, no chemistry perception.

**Intent.** molpack builds initial conformations; energy belongs to the
user's force field downstream. Conformer statistics come from priors
the caller supplies as geometric data (torsion-state weights, C∞,
persistence length, template values). Elements, bond orders, reactions,
hydrogens, and stereo centers enter only as user data at the boundary
(rigid groups, decoration-atom sets, bond labels, reaction templates).

**Never**

- Never make a solver, stage, term, or prior depend on the `ff`
  feature or on `molrs::ff` / `molrs::optimize`.
- Never derive a prior from a force field inside molpack.
- Never gate core behavior on chemistry perception (bond orders,
  aromaticity, element symbols); `molrs::perceive` may only produce
  data the user passes in.

**Derived guidance.** The in-loop optimizer (`GencanPack::with_optimizer`) is an optional
enhancement and must never become a solver dependency. Dependencies
follow policy (§ 6).

<!-- mol:law:id:solvers-are-peers -->
## P7. Solvers are peers, not layers

**Principle.** GENCAN, growth, and every future packing algorithm are
interchangeable implementations of one seam; they share lifecycle and
infrastructure, never each other's internals.

**Intent.** Shared: lifecycle stages ①②⑤, `PackContext` / `PackState`
(live run), the shared objective, frozen public `State`. Not shared:
drivers.

**Never**

- Never call `pgencan` / `run_phase` / `run_iteration` (or any other
  solver's driver) from another solver or stage. Reusing a numeric
  primitive (`spg::spgls`, `cg::cg_solve`) is not calling a driver.
- Never let a solver self-report `fdist` / `frest`; the final verdict
  comes from the shared objective on the final state — one ruler.
- Never chain a second algorithm silently on non-convergence; chaining
  is explicit user code (`with_restart`, `fixed_from`, a `Pipeline`).

<!-- mol:law:id:user-picks-method -->
## P8. The user picks the method

**Principle.** The caller chooses the algorithm by choosing the entry
or stage; molpack never decides for them.

**Never**

- Never infer the algorithm from the molecule (small molecule → rigid,
  polymer → growth).
- Never silently degrade to another method when the chosen one cannot
  handle a target; return a named error that says what to change.
- Never let a per-target choice be folded away by an `any()` / `first()`
  over targets.

<!-- mol:law:id:aa-and-cg -->
## P9. All-atom and coarse-grained alike

**Principle.** Every algorithm serves both all-atom and coarse-grained
templates.

**Never**

- Never hard-code an all-atom assumption: exclusion depth, angle
  handling, rotatable-bond detection, hydrogen identification are
  per-target data with all-atom defaults, not constants.
- Never key behavior on element symbols inside the core (see P6).

<!-- mol:law:id:molrs-abi-line -->
## P10. molrs ABI line

**Principle.** molpack exchanges `molrs_ffi` handle capsules with the
installed `molcrafts-molrs` wheel; both must embed the same molrs
**major.minor** (the minor line is the ABI version; see molrs
`docs/interop.md`).

**Never**

- Never hard-code a capsule name; take the versioned names
  (`molrs.FrameRef/<line>`) from `molrs_ffi::abi`.
- Never bypass the two gates: `molrs_capsule::check_abi` (`molrs._ffi_abi_token()`
  at extension init — the one import-time version check) and the versioned
  capsule names. Never add a second, metadata-based version check beside it.

<!-- add further project invariants below, one `<!-- mol:law:id:<slug> -->` each,
     using the same Principle / Intent / Never / Derived guidance template. -->
