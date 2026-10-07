//! What the growth engines refuse by name, and what they declare to the
//! stage chain. No successful growth is asserted here — that is the
//! driver's behaviour, not the engine's.

use super::*;

// ── CbmcGrow — named rejections (no silent degradation) ───────────────────

/// `Grow` on a `from_coords`-built target (template == None) is a named
/// error, not a silent fall-back to rigid-body packing.
#[test]
fn grow_rejects_target_without_template() {
    // The frozen prior surface: all three variants construct, and config /
    // prior are Debug + Clone.
    let priors = [
        TorsionPrior::Uniform,
        TorsionPrior::Template { kappa: 2.0 },
        TorsionPrior::States(vec![(0.0, 0.645), (2.094, 0.1775), (-2.094, 0.1775)]),
    ];
    let cfg = GrowConfig::new(priors[0].clone());
    let _ = format!("{:?}", cfg.clone());

    let coords: [[F; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.5, 0.0, 1.1],
        [4.0, 0.0, 1.1],
    ];
    let target = Target::from_coords(&coords, &[1.0; 4], 2);

    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(42)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("Grow without a template must be rejected");

    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::MissingTemplate,
                ..
            }
        ),
        "expected PackError::Grow(MissingTemplate), got {err:?}"
    );
}

/// `Grow` on a template frame with an atoms block but NO bonds block is a
/// named `NoBonds` error — chemistry comes from the bond graph, and molpack
/// does not guess it.
#[test]
fn grow_rejects_template_without_bonds() {
    let target = Target::new(chain_frame(5, 1.5, false), 2);

    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(42)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("Grow on a bond-less template must be rejected");

    let PackError::Grow { source, .. } = err else {
        panic!("expected PackError::Grow, got {err:?}");
    };
    assert!(
        matches!(source, GrowError::NoBonds),
        "expected GrowError::NoBonds, got {source:?}"
    );
    assert_eq!(
        source.to_string(),
        "the template frame carries no bonds; growth needs the bond graph — pack this \
         target with GencanPack or supply connectivity"
    );
}

/// A 2-atom bondless frame is `NoBonds`, not `TemplateTooSmall` or
/// `Disconnected` — refusal order is `NoBonds` before those variants.
#[test]
fn grow_rejects_two_atom_bondless_as_no_bonds() {
    let err = InternalTree::from_frame(
        &chain_frame(2, 1.53, false),
        &BondDistanceWeights::from_exclusion_depth(3),
    )
    .expect_err("a 2-atom bondless template must be refused");
    assert!(
        matches!(err, GrowError::NoBonds),
        "expected GrowError::NoBonds (not TemplateTooSmall or Disconnected), got {err:?}"
    );
}

/// An isolated atom next to a bonded pair is `Disconnected`.
#[test]
fn grow_rejects_disconnected_isolated_atom() {
    let frame = frame_from_parts(&zigzag_coords(3, 1.53), &[(0, 1)]);
    let err = InternalTree::from_frame(&frame, &BondDistanceWeights::from_exclusion_depth(3))
        .expect_err("an isolated atom must be refused");
    assert!(
        matches!(err, GrowError::Disconnected),
        "expected GrowError::Disconnected, got {err:?}"
    );
}

/// Two disjoint 3-atom chains: every atom has a bond, but `n_components > 1`.
#[test]
fn grow_rejects_two_disjoint_3atom_chains() {
    let bonds = [(0, 1), (1, 2), (3, 4), (4, 5)];
    let frame = frame_from_parts(&zigzag_coords(6, 1.53), &bonds);
    let err = InternalTree::from_frame(&frame, &BondDistanceWeights::from_exclusion_depth(3))
        .expect_err("two disjoint chains must be refused");
    assert!(
        matches!(err, GrowError::Disconnected),
        "expected GrowError::Disconnected, got {err:?}"
    );
}

/// No `atoms` block is `GrowError::NoAtomsBlock`.
#[test]
fn grow_rejects_no_atoms_block() {
    let mut frame = Frame::new();
    let mut bonds = Block::new();
    bonds
        .insert("atomi", Array1::from_vec(vec![0u32, 1]).into_dyn())
        .expect("atomi column");
    bonds
        .insert("atomj", Array1::from_vec(vec![1u32, 2]).into_dyn())
        .expect("atomj column");
    frame.insert("bonds", bonds);

    let err = InternalTree::from_frame(&frame, &BondDistanceWeights::from_exclusion_depth(3))
        .expect_err("no atoms block must be refused");
    assert!(
        matches!(err, GrowError::NoAtomsBlock),
        "expected GrowError::NoAtomsBlock, got {err:?}"
    );
}

/// An atoms block missing `z` is the same `NoAtomsBlock` as a missing block.
#[test]
fn grow_rejects_atoms_block_missing_z() {
    let coords = zigzag_coords(3, 1.53);
    let mut atoms = Block::new();
    for (name, k) in [("x", 0), ("y", 1)] {
        let col: Vec<F> = coords.iter().map(|p| p[k]).collect();
        atoms
            .insert(name, Array1::from_vec(col).into_dyn())
            .expect("coordinate column");
    }
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);
    let mut bonds = Block::new();
    bonds
        .insert("atomi", Array1::from_vec(vec![0u32, 1]).into_dyn())
        .expect("atomi column");
    bonds
        .insert("atomj", Array1::from_vec(vec![1u32, 2]).into_dyn())
        .expect("atomj column");
    frame.insert("bonds", bonds);

    let err = InternalTree::from_frame(&frame, &BondDistanceWeights::from_exclusion_depth(3))
        .expect_err("missing z column must be refused");
    assert!(
        matches!(err, GrowError::NoAtomsBlock),
        "expected GrowError::NoAtomsBlock, got {err:?}"
    );
}

/// Bond row 0 naming atom 99 on a 3-atom frame reports structured `BondOutOfRange`.
#[test]
fn grow_rejects_bond_out_of_range() {
    let frame = frame_from_parts(&zigzag_coords(3, 1.53), &[(0, 99)]);
    let err = InternalTree::from_frame(&frame, &BondDistanceWeights::from_exclusion_depth(3))
        .expect_err("out-of-range bond must be refused");
    match err {
        GrowError::BondOutOfRange { row, atom, n } => {
            assert_eq!(row, 0, "bonds-block row");
            assert_eq!(atom, 99, "out-of-range endpoint");
            assert_eq!(n, 3, "atom count from the frame");
        }
        other => panic!("expected GrowError::BondOutOfRange {{ row, atom, n }}, got {other:?}"),
    }
}

/// `Grow` on a bonded template with fewer than 3 atoms is a named
/// `TemplateTooSmall` error carrying the actual atom count.
#[test]
fn grow_rejects_tiny_template() {
    let target = Target::new(chain_frame(2, 1.5, true), 2);

    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(42)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("Grow on a 2-atom template must be rejected");

    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::TemplateTooSmall(2),
                ..
            }
        ),
        "expected PackError::Grow(TemplateTooSmall(2)), got {err:?}"
    );
}

/// A Grow target with neither a periodic box nor a cell is a named `NoBox`
/// error — growth needs a final volume from atom 0 (spec §4b/§7), and
/// molpack does not guess one.
#[test]
fn grow_requires_box() {
    let target = Target::new(chain_frame(12, 1.53, true), 2);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .run(&[target], 50)
        .expect_err("Grow without a periodic box or cell must be rejected");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::NoBox,
                ..
            }
        ),
        "expected PackError::Grow(NoBox), got {err:?}"
    );
}

/// A Grow target in a non-orthorhombic cell is a named `TriclinicCell` error
/// — the v1 `OverlapField` is orthorhombic-only (spec §10, Out of scope).
#[test]
fn grow_rejects_triclinic() {
    let target = Target::new(chain_frame(12, 1.53, true), 2);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_cell([20.0, 20.0, 20.0], [60.0, 60.0, 60.0], [true; 3])
        .run(&[target], 50)
        .expect_err("Grow in a 60°/60°/60° cell must be rejected");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::TriclinicCell,
                ..
            }
        ),
        "expected PackError::Grow(TriclinicCell), got {err:?}"
    );
}

/// `fixed_at` combined with Grow is a named `FixedTarget` error. The box is
/// valid, so the fixed placement is the only reason to refuse.
#[test]
fn grow_rejects_fixed_target() {
    let target = Target::new(chain_frame(12, 1.53, true), 1).fixed_at([5.0, 5.0, 5.0]);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [26.0; 3], [true; 3])
        .run(&[target], 50)
        .expect_err("fixed_at combined with Grow must be rejected");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::FixedTarget,
                ..
            }
        ),
        "expected PackError::Grow(FixedTarget), got {err:?}"
    );
}

/// Compile-level pin of the frozen `GrowConfig` builder surface (consuming
/// builders, `Target`-style; every knob has a default — values are NOT
/// asserted here, only the method names and signatures).
#[test]
fn grow_config_builder_chain() {
    let cfg = GrowConfig::new(TorsionPrior::Uniform)
        .with_trials(16)
        .with_selectivity(1.5)
        .with_soft_shell(0.8)
        .with_retract(6)
        .with_relax(20, 4)
        .with_soften_after(30)
        .with_min_hard_scale(0.85);
    let dbg = format!("{cfg:?}");
    assert!(!dbg.is_empty(), "GrowConfig must keep a useful Debug impl");
}

/// `with_density` together with an explicit box is a named
/// `DensityConflictsWithBox` error — one source of truth for the volume,
/// never a silent precedence rule (spec §7).
#[test]
fn density_conflicts_with_box() {
    let target = Target::new(chain_frame(5, 1.5, true), 2).with_mass(72.0);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(3)
        .with_tolerance(2.0)
        .with_density(0.05)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("with_density combined with with_periodic_box must be rejected");
    assert!(
        matches!(err, PackError::DensityConflictsWithBox),
        "expected PackError::DensityConflictsWithBox, got {err:?}"
    );
}

/// `with_density` on a target whose element masses cannot be resolved
/// (element "X", no `with_mass` override) is a named
/// `UnknownMass { target }` error — molpack does not guess masses (spec §7,
/// principle 3).
#[test]
fn density_needs_mass() {
    // chain_frame has no element column → Target::new resolves "X" masses.
    let target = Target::new(chain_frame(5, 1.5, true), 2);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(3)
        .with_tolerance(2.0)
        .with_density(0.05)
        .run(&[target], 20)
        .expect_err("density with an unresolvable target mass and no with_mass must be rejected");
    assert!(
        matches!(err, PackError::UnknownMass { target: 0 }),
        "expected PackError::UnknownMass {{ target: 0 }}, got {err:?}"
    );
}

// ── Section: the growth engine ─────────────────────────────────────────────

/// A ring template is refused by name, never silently grown open: a tree
/// decomposition would drop the ring-closing bonds without a word.
#[test]
fn ring_template_is_refused() {
    let coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.75, 1.3, 0.0]];
    let frame = frame_from_parts(&coords, &[(0, 1), (1, 2), (2, 0)]);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(frame.clone(), 2)], 60)
        .expect_err("a ring must be refused");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::RingTemplate,
                ..
            }
        ),
        "expected PackError::Grow(RingTemplate), got {err:?}"
    );
    let msg = format!("{err}");
    assert!(msg.contains("ring"), "named rejection, got: {msg}");
    let tree_err = InternalTree::from_frame(&frame, &BondDistanceWeights::from_exclusion_depth(3))
        .expect_err("a ring must not construct InternalTree");
    assert!(
        matches!(tree_err, GrowError::RingTemplate),
        "expected GrowError::RingTemplate, got {tree_err:?}"
    );
}

/// The seeded chain's named rejections and composability.
#[test]
fn seeded_run_contract() {
    let grown = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(9)
        .with_tolerance(1.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(chain_frame(5, 1.5, true), 2)], 60)
        .expect("grow stage runs");
    assert!(grown.converged);

    // Shape mismatch is a named rejection, not a scrambled pack.
    let err = GencanPack::new()
        .with_restart(&grown)
        .run(&[Target::new(chain_frame(5, 1.5, true), 3)], 10)
        .expect_err("a seed for 2 copies must refuse 3");
    assert!(
        matches!(err, PackError::SeedMismatch { .. }),
        "expected SeedMismatch, got {err:?}"
    );

    // The cell travels with the seed; a second box is the existing
    // mutual-exclusion error, never a silent precedence rule.
    let err = GencanPack::new()
        .with_restart(&grown)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(chain_frame(5, 1.5, true), 2)], 10)
        .expect_err("seed cell + declared box must be rejected");
    let msg = format!("{err}");
    assert!(msg.contains("mutually exclusive"), "got: {msg}");

    // A fixed matrix may be appended after the seeded free targets.
    let dimer = Target::from_coords(&[[0.0; 3], [1.5, 0.0, 0.0]], &[0.5, 0.5], 1)
        .with_centering(crate::CenteringMode::Off)
        .fixed_at([2.0, 2.0, 2.0]);
    let packed = GencanPack::new()
        .with_restart(&grown)
        .with_seed(3)
        .with_tolerance(1.0)
        .run(&[Target::new(chain_frame(5, 1.5, true), 2), dimer], 40)
        .expect("seeded run with an appended fixed matrix runs");
    assert_eq!(packed.natoms(), grown.natoms() + 2);
}

/// Degree > 4 is still a named tetrahedral refusal, not "branched staged".
#[test]
fn lattice_grow_rejects_degree_gt_4() {
    use crate::LatticeGrow;
    let mut coords = vec![[0.0; 3]];
    let mut bonds = Vec::new();
    for k in 0..5u32 {
        coords.push([(k as F + 1.0) * 1.5, 0.0, 0.0]);
        bonds.push((0, k + 1));
    }
    let err = LatticeGrow::new(TorsionPrior::Uniform)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(frame_from_parts(&coords, &bonds), 1)], 60)
        .expect_err("degree 5 must be named");
    let msg = format!("{err}");
    assert!(
        msg.contains("5") || msg.contains("tetrahedral") || msg.contains("degree"),
        "named degree>4, got: {msg}"
    );
    assert!(
        !msg.contains("staged for the lattice branch phase"),
        "{msg}"
    );
}

/// A mesh that does not overlap the cell leaves no diamond site.
#[test]
fn lattice_grow_empty_region_is_named() {
    use crate::LatticeGrow;
    use molrs::core::Polyhedron;
    use molrs::core::TriMesh;
    let cube = Polyhedron::new(TriMesh::from_triangles(&cube_tris([100.0; 3], [101.0; 3])))
        .expect("far cube");
    let err = LatticeGrow::new(TorsionPrior::Uniform)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(
            &[Target::new(chain_frame(6, 1.53, true), 1)
                .with_restraint(RegionRestraint(Arc::new(cube)))],
            4,
        )
        .expect_err("empty Region ∩ lattice must be named");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::LatticeRegionEmpty,
                ..
            }
        ),
        "expected PackError::Grow(LatticeRegionEmpty), got {err:?}"
    );
}

/// Same seed, same targets: bit-identical placements, for a pure grow and
/// for the grow → `GencanPack::with_restart` push-off chain.
#[test]
fn grow_and_push_off_are_deterministic_for_a_seed() {
    let grow = || {
        CbmcGrow::new(TorsionPrior::Uniform)
            .with_seed(9)
            .with_tolerance(1.0)
            .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
            .run(&[Target::new(chain_frame(5, 1.5, true), 2)], 60)
            .expect("grow stage runs")
    };
    let push_off = |seed: &crate::State| {
        GencanPack::new()
            .with_restart(seed)
            .with_seed(3)
            .with_tolerance(1.0)
            .run(&[Target::new(chain_frame(5, 1.5, true), 2)], 40)
            .expect("push-off runs")
    };
    let same = |a: &crate::State, b: &crate::State| {
        assert_eq!(a.natoms(), b.natoms());
        for (pa, pb) in a.positions().iter().zip(b.positions().iter()) {
            for k in 0..3 {
                assert_eq!(pa[k].to_bits(), pb[k].to_bits());
            }
        }
    };
    let (g1, g2) = (grow(), grow());
    same(&g1, &g2);
    same(&push_off(&g1), &push_off(&g2));
}

// ── Section: the seam markers the two growth stages declare ───────────────
//
// Owner-side half of acceptance ac-008 (`stage-pipeline-04-stage`): what
// `GrowStage` and `LatticeStage` declare on the stage seam belongs to their
// owner — `stage::tests` verifies the seam's own behaviour with fakes and
// boots no real algorithm. Both stages build their placements from nothing
// and return with every free molecule placed, so both read
// `Placed::None -> Placed::All`.

/// Two templates keep two tables: GrowStage must not min-fold them into
/// one skip set (law P8). Depth 1 excludes self+1-2; depth 3 excludes
/// self+1-2/1-3/1-4. Residual classes on the same 5-bead stick differ
/// the same way (scored 2.0 Å vs 4.0 Å).
#[test]
fn grow_stage_two_targets_keep_distinct_tables() {
    use crate::grow::driver::GrowStage;
    use molrs::core::SimBox;
    use ndarray::array;

    let coords: Vec<[F; 3]> = (0..5).map(|i| [i as F, 0.0, 0.0]).collect();
    let frame = frame_from_parts(&coords, &chain_bonds(5));
    let shallow = Target::new(frame.clone(), 1)
        .with_special_bonds(BondDistanceWeights::from_exclusion_depth(1));
    let deep =
        Target::new(frame, 1).with_special_bonds(BondDistanceWeights::from_exclusion_depth(3));

    let stage = GrowStage::from_targets(
        &[shallow.clone(), deep.clone()],
        &GrowConfig::new(TorsionPrior::Uniform),
        7,
    )
    .expect("two templates with distinct special-bonds tables must both compile");

    assert_eq!(
        stage.tree(0).exclusions(0),
        &[0u32, 1][..],
        "depth-1 species excludes self + 1-2 only"
    );
    assert_eq!(
        stage.tree(1).exclusions(0),
        &[0u32, 1, 2, 3][..],
        "depth-3 species excludes self + 1-2/1-3/1-4"
    );
    assert_ne!(
        stage.tree(0).exclusions(0),
        stage.tree(1).exclusions(0),
        "GrowStage must not min-fold two tables into one skip set"
    );

    let cell = SimBox::ortho(
        array![100.0, 100.0, 100.0],
        array![0.0, 0.0, 0.0],
        [false; 3],
    )
    .expect("orthorhombic box");
    let r1 = IntraResidual::from_targets(std::slice::from_ref(&shallow), &coords, &cell);
    let r3 = IntraResidual::from_targets(std::slice::from_ref(&deep), &coords, &cell);
    assert_eq!(
        r1.scored.to_bits(),
        (2.0 as F).to_bits(),
        "depth-1 scored 1-3 at 2.0 Å"
    );
    assert_eq!(
        r3.scored.to_bits(),
        (4.0 as F).to_bits(),
        "depth-3 scored 1-5 at 4.0 Å"
    );
}

/// The continuum growth stage requires nothing placed and guarantees
/// everything placed.
#[test]
fn grow_stage_requires_none_guarantees_all() {
    use crate::grow::driver::GrowStage;
    use crate::{Placed, Stage};

    let targets = [Target::new(chain_frame(6, 1.53, true), 2)];
    let stage = GrowStage::from_targets(&targets, &GrowConfig::new(TorsionPrior::Uniform), 7)
        .expect("a bead-chain target with a template builds a growth stage");

    assert_eq!(
        stage.requires().placed,
        Placed::None,
        "growth constructs its own placements, so it requires nothing placed"
    );
    assert_eq!(
        stage.guarantees().placed,
        Placed::All,
        "growth returns with every free molecule placed"
    );
}

/// The diamond-lattice growth stage declares the same two markers.
#[test]
fn lattice_stage_requires_none_guarantees_all() {
    use crate::grow::lattice::LatticeConfig;
    use crate::grow::lattice::LatticeStage;
    use crate::{Placed, Stage};

    let targets = [Target::new(chain_frame(6, 1.53, true), 2)];
    let stage = LatticeStage::from_targets(&targets, &LatticeConfig::new(TorsionPrior::Uniform), 7)
        .expect("a bead-chain target with a template builds a lattice stage");

    assert_eq!(
        stage.requires().placed,
        Placed::None,
        "lattice growth decorates its own placements onto the lattice"
    );
    assert_eq!(
        stage.guarantees().placed,
        Placed::All,
        "lattice growth returns with every free molecule placed"
    );
}
