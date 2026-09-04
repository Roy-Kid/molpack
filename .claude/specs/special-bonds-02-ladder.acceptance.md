---
slug: special-bonds-02-ladder
created: 2026-09-04
criteria:
  - id: ac-001
    summary: Three clocks with three named readers; deadends is gone
    type: code
    pass_when: |
      Chain has deadends_total, deadend_streak, rungs_earned and no field
      named deadends. Successful commit sets deadend_streak = 0 only.
      The retract 1<< expression reads deadend_streak and does not mention
      deadends_total. Softening uses rung_due(total, rungs_earned, soften_after)
      with saturating arithmetic. Force-place uses force_due(..., streak, ...)
      evaluated before this round's streak increment.
    status: verified
    last_checked: 2026-09-04
  - id: ac-002
    summary: Probe::Blocked stays a unit variant; probe stays first-hit
    type: code
    pass_when: |
      Probe::Blocked is still a unit variant. probe still returns Blocked on
      the first hard-core hit and does not call block_kind. tests/grow.rs
      asserts of Probe::Blocked compile unchanged.
    status: verified
    last_checked: 2026-09-04
  - id: ac-003
    summary: block_kind splits self-blocked from inter-chain
    type: runtime
    pass_when: |
      cargo test -p molcrafts-molpack --lib --tests -- field_block_kind_self_vs_inter
      is green with the hard-coded SelfBlocked / InterChain / None / SelfBlocked-wins
      cases in the spec Testing strategy.
    status: verified
    last_checked: 2026-09-04
  - id: ac-004
    summary: Softening ladder fires through intermittent commits
    type: runtime
    pass_when: |
      grow_ladder_fires_after_intermittent_commits is green: soften_after = 2
      yields radscale < 1.0 at least once, every drop is one 0.97 rung.
      In-module retract_depth / force_due / rung_due tests are green, including
      force_due(0.8, 0.8, 0, 50) == false.
    status: verified
    last_checked: 2026-09-04
  - id: ac-005
    summary: Existing ladder pins stay green unedited
    type: runtime
    pass_when: |
      assert_radscale_ladder, grow_softening_needs_repeated_dead_ends,
      grow_dense_strict_core_terminates_unconverged, and
      grow_cg_kremer_grest_c_inf are bit-for-bit the pre-spec source and green.
      KG is run first.
    status: verified
    last_checked: 2026-09-04
  - id: ac-006
    summary: propose returns DeadEnd; commit classifies before rollback
    type: code
    pass_when: |
      score_atoms returns Result<(Vec<PlacedAtom>, F), DeadEnd>. propose
      returns Result<Proposal, DeadEnd> with DeadEnd = Overlap(BlockKind) |
      Restraint. commit returns Result<(), BlockKind> and classifies in the
      Probe::Blocked arm before rollback; it does not expect() after a
      cleared field. Across candidates, SelfBlocked wins.
    status: verified
    last_checked: 2026-09-04
  - id: ac-007
    summary: Auhl 0.8 floor; hybrid is not B1
    type: code
    pass_when: |
      GrowConfig::new still defaults min_hard_scale to 0.8. PackResult has no
      new fields. StageOutcome::new still takes (bool, usize). driver.rs still
      has a single run-local hard_scale: F and the min_hard_scale fold.
    status: verified
    last_checked: 2026-09-04
  - id: ac-008
    summary: Regression example embeds hard-coded ladder goldens
    type: runtime
    pass_when: |
      regressions/special-bonds-02-ladder.md names CbmcGrow, StepInfo.radscale,
      OverlapField::block_kind, embeds 0.97 / 0.8 / SelfBlocked / InterChain
      and the 2026-09-04 defect observation, with no third-party runtime.
    status: verified
    last_checked: 2026-09-04
  - id: ac-009
    summary: Rustdoc states three readers and the new soften_after clock
    type: docs
    pass_when: |
      cargo doc --no-deps reports 0 warnings. with_soften_after rustdoc is
      cumulative not consecutive. hard_scale is documented as dimensionless.
    status: pending
out_of_scope:
  - "PackResult residual (04), Topology sink (03), Target table (05), Python (06)"
  - "grow-axes B1 per-chain recoverable hard_scale"
  - "Unfolding min_hard_scale.min() / serial.first()"
---

# Acceptance — special-bonds-02-ladder

Done means the per-chain softening ladder fires through intermittent commits,
the three clock readers cannot be cross-wired, commit classifies before
rollback, and the four named ladder pins stay unedited.
