# Getting Started

This chapter walks through the minimum you need to run a packing job
and the most common customization points. Read the reference docs of
individual types on the crate landing page for full details.

## CLI

The `molpack` binary accepts the same `.inp` script format as the
original Packmol, making it a drop-in replacement for the command-line
workflow. See [Install](install.md) for setup.

### Running a script

```bash
# File argument — paths inside the script resolve relative to its directory
molpack mixture.inp

# Stdin — compatible with the classic Packmol invocation
molpack < mixture.inp
cat mixture.inp | molpack
```

### Script format

The `.inp` format is line-oriented. Lines starting with `#` are comments.
Global keywords come first; each molecule type occupies a
`structure … end structure` block.

```text
# Global settings
tolerance 2.0          # minimum atom–atom distance in Å (default 2.0)
seed 42                # random seed for reproducibility (optional)
filetype pdb           # input format for all structure files (optional)
output packed.pdb      # output file; format inferred from extension
nloop 400              # max outer-loop iterations (default 400)
pbc 40.0 40.0 40.0     # optional: fully-periodic box [0,0,0]..[40,40,40]

# One block per molecule type
structure water.pdb
  number 1000
  inside box 0. 0. 0. 40. 40. 40.
end structure

structure urea.pdb
  number 400
  inside box 0. 0. 0. 40. 40. 40.
end structure
```

**`pbc` keyword** — declares a fully-periodic box (every axis wraps),
matching Packmol's `pbc`. Two forms are accepted:

- `pbc X Y Z` — `min = [0, 0, 0]`, `max = [X, Y, Z]`
- `pbc X0 Y0 Z0  X1 Y1 Z1` — explicit `min`/`max`

When a script sets `pbc`, the packer's cell grid is built directly from
that box. Without either a `pbc` directive or an `inside`-style
restraint the packer falls back to inferring a box from the initial
random placement (half-width defaults to 1000 Å), which drives the cell
grid to tens of millions of cells — so always give the packer a
spatial constraint.

**Unknown keywords are rejected.** The parser returns a
[`ScriptError::UnknownKeyword`](crate::script::ScriptError) rather than
silently dropping them — a silent `pbc` drop used to burn 40+ GB of
RAM before the packer even started.

**Restraint keywords** (all Packmol restraint types are supported):

| Keyword | Scope | Description |
|---|---|---|
| `inside box x0 y0 z0 x1 y1 z1` | molecule / atoms | Axis-aligned box |
| `inside sphere cx cy cz r` | molecule / atoms | Sphere |
| `outside sphere cx cy cz r` | molecule / atoms | Exclusion sphere |
| `over plane nx ny nz d` | molecule / atoms | Half-space (above) |
| `below plane nx ny nz d` | molecule / atoms | Half-space (below) |

Per-atom-subset restraints use an `atoms … end atoms` sub-block:

```text
structure palmitoil.pdb
  number 10
  inside box 0. 0. 0. 40. 40. 14.
  atoms 31 32          # 1-based atom indices
    below plane 0. 0. 1. 2.
  end atoms
  atoms 1 2
    over plane 0. 0. 1. 12.
  end atoms
end structure
```

Fixed-position placement:

```text
structure protein.pdb
  number 1
  center
  fixed 0. 0. 0. 0. 0. 0.   # x y z  euler_x euler_y euler_z
end structure
```

### Extended format support

molpack extends Packmol's PDB/XYZ input support with the additional readers
from molrs-io. The `filetype` keyword or file extension selects the format:

| Format | Read | Write | Extension / `filetype` keyword |
|---|---|---|---|
| PDB | ✓ | ✓ | `.pdb` / `pdb` |
| XYZ | ✓ | ✓ | `.xyz` / `xyz` |
| SDF / MOL | ✓ | — | `.sdf`, `.mol` / `sdf` |
| LAMMPS dump | ✓ | ✓ | `.lammpstrj` / `lammps_dump` |
| LAMMPS data | ✓ | — | `.data` / `lammps_data` |

The output format is always inferred from the extension of the `output` path.

### Running the canonical Packmol examples via CLI

The five canonical workloads each ship with a `.inp` file:

```bash
molpack examples/pack_mixture/mixture.inp
molpack examples/pack_bilayer/bilayer-comment.inp
molpack examples/pack_interface/interface.inp
molpack examples/pack_solvprotein/solvprotein.inp
molpack examples/pack_spherical/spherical.inp
```

---

## Rust library

Add the dependency (see [Install](install.md) for the full toolchain).
Then write a packing job in code.

### First pack: one molecule type in a box

```rust,no_run
use molpack::{InsideBoxRestraint, Molpack, Target};

let water_positions = [
    [0.0, 0.0, 0.0],
    [0.96, 0.0, 0.0],
    [-0.24, 0.93, 0.0],
];
let water_radii = [1.52, 1.20, 1.20];

let target = Target::from_coords(&water_positions, &water_radii, 100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3]));

let result = Molpack::new()
    .with_tolerance(2.0)
    .with_precision(0.01)
    .with_seed(42)
    .pack(&[target], 200)?;

println!("converged = {}, natoms = {}", result.converged, result.natoms());
# Ok::<(), molpack::PackError>(())
```

Arguments to [`Molpack::pack`](crate::Molpack::pack):

- `&[Target]` — one entry per molecule type.
- `max_loops: usize` — outer-iteration budget per phase.

Seed and every other tuning knob live on the builder (`.with_seed(n)`,
`.with_tolerance(t)` etc.). Returns a [`PackResult`](crate::PackResult)
with the output coordinates as a `molrs::Frame` plus convergence
diagnostics.

## Restraint scopes

Restraints can be attached at three granularities. All of them use the
same [`Restraint`](crate::Restraint) trait underneath.

**Per-target, all atoms** — the most common case:

```rust,no_run
# use molpack::{InsideBoxRestraint, Target};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
let target = Target::from_coords(pos, rad, 100)
    .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3]));
```

**Per-target, atom subset** — restrict to specific atoms of every
molecule copy (0-based, Rust convention — if you are porting from a
Packmol `.inp` file, subtract 1):

```rust,no_run
# use molpack::{BelowPlaneRestraint, Target};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
let target = Target::from_coords(pos, rad, 100)
    .with_atom_restraint(&[0, 1], BelowPlaneRestraint::new([0.0, 0.0, 1.0], 5.0));
```

**Global, all targets** — broadcast. Semantically equivalent to
calling `with_restraint` on every target:

```rust,no_run
# use molpack::{InsideSphereRestraint, Molpack, Target};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
# let t_a = Target::from_coords(pos, rad, 100);
# let t_b = Target::from_coords(pos, rad, 100);
let result = Molpack::new()
    .with_global_restraint(InsideSphereRestraint::new([20.0; 3], 30.0))
    .with_seed(42)
    .pack(&[t_a, t_b], 200)?;
# Ok::<(), molpack::PackError>(())
```

The scope-equivalence law is formal: `molpack.with_global_restraint(r)`
is implemented as `for t in targets { t.with_restraint(r.clone()) }`.
There is no separate "global-restraint" storage path.

## Handlers (progress observation)

Attach [`Handler`](crate::Handler) implementations for progress reporting,
trajectory output, or early-stop logic. Built-ins:

```rust,no_run
# use molpack::{EarlyStopHandler, InsideBoxRestraint, Molpack, ProgressHandler, Target, XYZHandler};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
# let target = Target::from_coords(pos, rad, 100).with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3]));
let mut packer = Molpack::new()
    .with_handler(ProgressHandler::new())
    .with_handler(XYZHandler::new("traj.xyz", /* every = */ 10))
    .with_handler(EarlyStopHandler::new(/* threshold = */ 1e-4));
# let _ = packer.with_seed(42).pack(&[target], 200);
```

Write your own — see the [`extending`](crate::extending) module.

## Relaxers (in-loop conformation sampling)

Flexible molecules benefit from torsion-MC relaxation between outer
optimizer calls. Attach a [`Relaxer`](crate::Relaxer) to a target:

```rust,no_run
# use molrs::molgraph::MolGraph;
# use molpack::{InsideSphereRestraint, Target, TorsionMcRelaxer};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
# let graph = MolGraph::new();
let target = Target::from_coords(pos, rad, 1)  // relaxers require count == 1
    .with_restraint(InsideSphereRestraint::new([0.0; 3], 20.0))
    .with_relaxer(
        TorsionMcRelaxer::new(&graph)
            .with_temperature(0.5)
            .with_steps(20)
    );
```

## Free versus periodic boundary

Free boundary is the default. There are two ways to declare PBC:

**On the packer** (fully periodic box, equivalent to the script's
`pbc` keyword):

```rust,no_run
# use molpack::Molpack;
let packer = Molpack::new().with_periodic_box([0.0; 3], [30.0; 3]);
```

**On an `InsideBoxRestraint`** (per-axis control, e.g. for slab
geometries):

```rust,no_run
# use molpack::{InsideBoxRestraint, Target};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
// Fully-periodic (Packmol-style 3D PBC) via a restraint:
let target = Target::from_coords(pos, rad, 100).with_restraint(
    InsideBoxRestraint::new([0.0; 3], [30.0; 3], [true; 3]),
);
```

Slab geometries with only some axes wrapping (e.g. X/Y periodic
membrane, Z confined):

```rust,no_run
# use molpack::{InsideBoxRestraint, Target};
# let (pos, rad) = (&[[0.0; 3]][..], &[1.0][..]);
let target = Target::from_coords(pos, rad, 100).with_restraint(
    InsideBoxRestraint::new([0.0; 3], [30.0; 3], [true, true, false]),
);
```

When a `pack()` call resolves PBC from both a packer-level
`with_periodic_box` **and** a restraint-level declaration, the two
must match exactly (bounds + per-axis flags) or
[`PackError::ConflictingPeriodicBoxes`](crate::PackError) is returned.
Zero-length axes return [`PackError::InvalidPBCBox`](crate::PackError).

## Running the canonical examples

Five Packmol-equivalent workloads are checked in under `examples/`:

```bash
cargo run -p molcrafts-molpack --release --example pack_mixture
cargo run -p molcrafts-molpack --release --example pack_bilayer
cargo run -p molcrafts-molpack --release --example pack_interface
cargo run -p molcrafts-molpack --release --example pack_solvprotein
cargo run -p molcrafts-molpack --release --example pack_spherical
```

Two heavier ad-hoc demonstrations (not in the regression suite):

```bash
cargo run -p molcrafts-molpack --release --example mc_fold_chain
cargo run -p molcrafts-molpack --release --example pack_polymer_vesicle
```

Set `MOLRS_PACK_EXAMPLE_PROGRESS=1` to enable the built-in
`ProgressHandler`. Set `MOLRS_PACK_EXAMPLE_XYZ=1` to dump a trajectory
under `out/`.

## Python bindings

A PyO3 binding ships under `python/`:

```bash
pip install molcrafts-molpack molcrafts-molrs
```

```python
import molrs
from molpack import InsideBoxRestraint, Molpack, Target

frame = molrs.read_pdb("water.pdb")

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
)
result = (
    Molpack()
    .with_tolerance(2.0)
    .with_seed(42)
    .pack([water], max_loops=200)
)
```

The Python API mirrors the Rust builder surface one-for-one. The
binding is I/O-free — use
[`molcrafts-molrs`](https://pypi.org/project/molcrafts-molrs/) (or any
Frame-compatible loader) to read templates, and write `result.frame`
back out with the same library.

The standalone Python documentation site lives under
[`python/docs/`](https://github.com/MolCrafts/molpack/tree/master/python/docs);
runnable examples are in
[`python/examples/`](https://github.com/MolCrafts/molpack/tree/master/python/examples).

## Where to go next

- [Concepts](concepts.md) — the abstractions in one place.
- [Architecture](architecture.md) — system design, data flow, loops, hot path.
- [Extending](extending.md) — write your own `Restraint` / `Region` /
  `Handler` / `Relaxer`.
