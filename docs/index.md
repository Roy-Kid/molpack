# molpack

Molecular packing in pure Rust, with Python bindings. Part of the
[molrs](https://github.com/MolCrafts/molrs) toolkit.

molpack builds initial configurations for molecular dynamics: give it a set
of molecular templates, per-type counts, and geometric restraints (boxes,
spheres, half-spaces, fixed placements), and it produces non-overlapping
coordinates that satisfy all constraints.

## Install

=== "CLI"

    ```bash
    cargo install molcrafts-molpack --features cli
    ```

=== "Rust library"

    ```bash
    cargo add molcrafts-molpack
    ```

=== "Python"

    ```bash
    pip install molcrafts-molpack molcrafts-molrs
    ```

See **[Install](install.md)** for prerequisites, available wheels, and the
build-from-source flow.

## Quick start

=== "CLI"

    The `molpack` binary accepts Packmol's `.inp` script format, so it works
    as a drop-in replacement on the command line. Beyond PDB/XYZ it also
    reads SDF/MOL, LAMMPS dump, and LAMMPS data.

    ```bash
    molpack mixture.inp         # file argument
    molpack < mixture.inp       # stdin
    cat mixture.inp | molpack   # pipe
    ```

=== "Rust"

    ```rust
    use molpack::{InsideBoxRestraint, Molpack, Target};

    let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    let radii = [1.52, 1.20, 1.20];

    let target = Target::from_coords(&positions, &radii, 100)
        .with_name("water")
        .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [false; 3]));

    let frame = Molpack::new().pack(&[target], 200)?;
    ```

=== "Python"

    ```python
    import molrs
    from molpack import InsideBoxRestraint, Molpack, Target

    frame = molrs.read_pdb("water.pdb")

    water = (
        Target(frame, count=100)
        .with_name("water")
        .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
    )
    frame = Molpack().pack([water], max_loops=200)
    ```

## Documentation map

This guide has eight chapters. The four long-form chapters
(getting started, concepts, architecture, extending) are also embedded in
the crate rustdoc as `molpack::getting_started` and friends.

**Get Started**

- **[Install](install.md)** — CLI, Rust crate, and Python binding.
- **[Getting Started](getting_started.md)** — your first pack, end to end.

**Guide**

- **[Concepts](concepts.md)** — targets, restraints, regions, handlers, relaxers.
- **[Examples](examples.md)** — the five canonical workloads.
- **[Packmol parity](packmol_parity.md)** — what is matched, and how it is verified.

**Internals**

- **[Architecture](architecture.md)** — module map, data flow, loops, hot path.
- **[Extending](extending.md)** — write your own `Restraint` / `Region` / `Handler` / `Relaxer`.

This page is the eighth chapter.

## License

BSD-3-Clause. See [LICENSE](https://github.com/MolCrafts/molpack/blob/master/LICENSE).
