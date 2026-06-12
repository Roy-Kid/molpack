---
slug: profile-distribution-restraints-01-split
criteria:
  - id: ac-001
    summary: Existing restraint test suite passes unchanged after the split
    type: code
    evaluator_hint: "cargo test -p molcrafts-molpack --lib --tests"
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests` passes with
      tests/restraint.rs unmodified (same imports from molpack:: and same
      assertions, including assert_gradient_opposes_violation), proving
      f/fg behaviour is preserved across all 14 geometric restraints.
    status: verified
  - id: ac-002
    summary: External crate::restraint paths resolve with no import churn outside the module
    type: code
    evaluator_hint: "cargo clippy -p molcrafts-molpack --all-targets -- -D warnings"
    pass_when: |
      `cargo clippy -p molcrafts-molpack --all-targets -- -D warnings`
      is clean AND the only files changed under src/ besides the new
      restraint/ dir are none — every site naming crate::restraint::<Type>
      or molpack::<Type> (lib.rs:130, validation.rs, target.rs, region.rs,
      packer.rs, context/pack_context.rs) compiles without edits.
    status: verified
  - id: ac-003
    summary: Both new module files stay within the 800-LOC hard cap
    type: code
    evaluator_hint: "wc -l src/restraint/mod.rs src/restraint/geometric.rs"
    pass_when: |
      src/restraint/mod.rs and src/restraint/geometric/*.rs each report
      <= 800 lines, and the standalone src/restraint.rs no longer exists.
      (Note: geometric split into geometric/{mod,bounded,surface}.rs to
      honor the 800-LOC cap — a single geometric.rs would have been 891 LOC.)
    status: verified
---

# Acceptance criteria

- **ac-001** — The behaviour-preservation gate. The regression net (`tests/restraint.rs`) is the load-bearing proof: because the `impl Restraint` bodies are moved verbatim, the suite must pass with zero edits to its imports or assertions.
- **ac-002** — The path-stability gate. The whole point of `mod.rs` re-exporting `geometric::*` is that no consumer outside the `restraint` module changes. A clean `-D warnings` clippy across all targets confirms every `crate::restraint::<Type>` and crate-root `molpack::<Type>` still resolves.
- **ac-003** — The LOC-budget gate that motivated the refactor. Both new files must sit under the 800-line CLAUDE.md cap, and the original `src/restraint.rs` must be gone (directory module, not a stray sibling file).
