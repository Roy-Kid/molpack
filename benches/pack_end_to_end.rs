//! Regression alarm for the full `Molpack::pack` pipeline.
//!
//! Drives a small, fully-synthesized two-species mixture through the public
//! `Molpack::pack` entry point — initial placement, the three-phase GENCAN
//! loop, and frame assembly — so a regression anywhere in the end-to-end path
//! shows up. This is the coarse "catastrophic-regression" alarm that
//! complements the per-kernel microbenches in the sibling bench files.
//!
//! **Kept deliberately small and fast** (a dilute ~28-molecule / 124-atom box,
//! few outer loops) so it is a regression guard, not a throughput measurement,
//! and finishes well under a second. Geometry is built in-process via the
//! public `Target` API — no example PDBs, so the bench needs no `io` feature.
//! The larger canonical workloads live in `examples/pack_*` and the
//! `examples_batch` regression test.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use molpack::{F, InsideBoxRestraint, InsideCellRegion, Molpack, RegionRestraint, Target};

const SEED: u64 = 42;
const BOX_SIDE: F = 25.0;
const MAX_LOOPS: usize = 20;

/// A compact 3-atom water-like species.
fn water_target() -> Target {
    let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let radii = [1.5, 1.2, 1.2];
    Target::from_coords(&positions, &radii, 20)
        .with_name("water")
        .with_restraint(InsideBoxRestraint::new([0.0; 3], [BOX_SIDE; 3], [false; 3]))
}

/// A small 8-atom rigid species with deterministic ~5 Å-wide geometry.
fn blob_target() -> Target {
    let mut rng = SmallRng::seed_from_u64(0xA);
    let positions: Vec<[F; 3]> = (0..8)
        .map(|_| {
            [
                (rng.random::<F>() * 2.0 - 1.0) * 2.5,
                (rng.random::<F>() * 2.0 - 1.0) * 2.5,
                (rng.random::<F>() * 2.0 - 1.0) * 2.5,
            ]
        })
        .collect();
    let radii = vec![1.5; positions.len()];
    Target::from_coords(&positions, &radii, 8)
        .with_name("blob")
        .with_restraint(InsideBoxRestraint::new([0.0; 3], [BOX_SIDE; 3], [false; 3]))
}

/// The same two species confined to a hexagonal cell of equal volume.
///
/// Equal volume, not equal edge length: a hexagonal cell of side `s` holds
/// `s³·sin(120°)`, so matching the edge would make the triclinic run denser and
/// the two timings incomparable.
fn hexagonal_targets() -> [Target; 2] {
    let side = BOX_SIDE / (120.0 as F).to_radians().sin().powf(1.0 / 3.0);
    let cell =
        InsideCellRegion::from_lengths_angles([side, side, side], [90.0, 90.0, 120.0], [false; 3])
            .expect("hexagonal cell");

    let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let radii = [1.5, 1.2, 1.2];
    let water = Target::from_coords(&positions, &radii, 20)
        .with_name("water")
        .with_restraint(RegionRestraint(cell.clone()));

    let mut rng = SmallRng::seed_from_u64(0xA);
    let blob_positions: Vec<[F; 3]> = (0..8)
        .map(|_| {
            [
                (rng.random::<F>() * 2.0 - 1.0) * 2.5,
                (rng.random::<F>() * 2.0 - 1.0) * 2.5,
                (rng.random::<F>() * 2.0 - 1.0) * 2.5,
            ]
        })
        .collect();
    let blob = Target::from_coords(&blob_positions, &[1.5; 8], 8)
        .with_name("blob")
        .with_restraint(RegionRestraint(cell));

    [water, blob]
}

fn run(targets: &[Target], parallel: bool) {
    let mut packer = Molpack::new()
        .with_tolerance(2.0)
        .with_precision(0.01)
        .with_parallel_eval(parallel)
        .with_seed(SEED);
    std::hint::black_box(packer.pack(targets, MAX_LOOPS).expect("pack"));
}

fn bench_pack_end_to_end(c: &mut Criterion) {
    let mut group = c.benchmark_group("pack_end_to_end");
    group.sample_size(10);
    group.measurement_time(Duration::from_millis(500));

    let targets = [water_target(), blob_target()];

    let hexagonal = hexagonal_targets();

    group.bench_function("mixture_small", |b| b.iter(|| run(&targets, false)));

    // Same workload in a hexagonal cell of equal volume. Not a comparison of
    // "which is faster" — the triclinic minimum image is strictly more work —
    // but a permanent baseline for a code path that had none, since it could
    // not be expressed before.
    group.bench_function("mixture_small_hexagonal", |b| {
        b.iter(|| run(&hexagonal, false))
    });

    // Serial vs rayon over cells. Forward neighbours are now selected by cell
    // index rather than by a fixed set of 13 offset directions, so per-cell
    // neighbour counts run from 0 to 26 where they used to be a flat 13. Total
    // work is unchanged (13 per cell on average) but its distribution across
    // rayon tasks is not, and this pair is what would show that going wrong.
    group.bench_function("mixture_small_parallel", |b| b.iter(|| run(&targets, true)));
    group.bench_function("mixture_small_hexagonal_parallel", |b| {
        b.iter(|| run(&hexagonal, true))
    });

    group.finish();
}

criterion_group!(benches, bench_pack_end_to_end);
criterion_main!(benches);
use rand::RngExt;
use rand::SeedableRng;
use rand::rngs::SmallRng;
