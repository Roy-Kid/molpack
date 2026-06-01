<div align="center">

<h1>
  <img src=".github/assets/moko.svg" alt="" height="48" align="absmiddle">
  &nbsp;molpack
</h1>

<p><strong>Packmol-grade molecular packing — Rust core, Python bindings, and a Packmol-compatible CLI.</strong></p>

<p>
  <a href="https://github.com/MolCrafts/molpack/actions/workflows/ci.yml"><img src="https://img.shields.io/github/actions/workflow/status/MolCrafts/molpack/ci.yml?style=flat-square&logo=githubactions&logoColor=white&label=CI" alt="CI"></a>
  <a href="https://crates.io/crates/molcrafts-molpack"><img src="https://img.shields.io/crates/v/molcrafts-molpack?style=flat-square&logo=rust&logoColor=white" alt="Crates.io"></a>
  <a href="https://pypi.org/project/molcrafts-molpack/"><img src="https://img.shields.io/pypi/v/molcrafts-molpack?style=flat-square&logo=pypi&logoColor=white&label=PyPI" alt="PyPI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-BSD--3--Clause-18432B?style=flat-square" alt="License"></a>
</p>

<p>
  <a href="https://molcrafts.github.io/molpack/"><b>Documentation</b></a> &nbsp;&middot;&nbsp;
  <a href="#quick-start"><b>Quick start</b></a> &nbsp;&middot;&nbsp;
  <a href="#molcrafts-ecosystem"><b>Ecosystem</b></a>
</p>

</div>

molpack produces non-overlapping arrangements of many molecule types with copy
counts and geometric restraints, using a faithful Rust port of Packmol's
GENCAN-driven three-phase algorithm. It ships as a Rust library, a Python wheel,
and a CLI that reads Packmol `.inp` scripts.

> **Under active development.** Public APIs may change between minor releases.

## Capabilities

| Module | Capability |
|---|---|
| `packer` | `Molpack` builder + outer-loop driver; runs the three-phase pack and returns a `PackResult` |
| `gencan` | Bound-constrained GENCAN optimizer (conjugate gradient + spectral projected gradient) |
| `target` | `Target` molecule descriptor — copy counts, centering modes, fixed placements, Euler angles |
| `restraint` | `Restraint` trait + 14 concrete soft-penalty types (box, cube, sphere, ellipsoid, cylinder, plane, gaussian; inside/outside/above/below variants) |
| `region` | `Region` signed-distance predicates with `And` / `Or` / `Not` combinators, liftable to restraints |
| `relaxer` | In-loop reference-geometry relaxers; `TorsionMcRelaxer` for flexible-molecule torsion MC |
| `handler` | Progress-callback handlers — `Progress`, `EarlyStop`, `XYZ`, `Null` |
| `validation` | Post-pack constraint checking — `ValidationReport`, `ViolationMetrics` |
| `script` | Packmol `.inp` parser + lowering to targets, with format-aware file I/O (`io` feature) |
| `bin/molpack` | `molpack` CLI — drop-in `.inp` runner accepting file-arg or stdin input (`cli` feature) |

File formats (via the `io` feature): reads PDB, XYZ, SDF/MOL, LAMMPS dump, and
LAMMPS data; writes PDB, XYZ, and LAMMPS dump.

## Install

```bash
# Rust library
cargo add molcrafts-molpack

# CLI (Packmol-compatible .inp runner)
cargo install molcrafts-molpack --features cli

# Python
pip install molcrafts-molpack molcrafts-molrs
```

The Rust library depends on `molcrafts-molrs-core` for shared data structures;
the `io` and `cli` features additionally pull in `molcrafts-molrs-io`. The Python
wheel loads frames through the `molcrafts-molrs` Python package.

## Quick start

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

Rust API, CLI usage, the `.inp` keyword reference, and more worked examples are
in the documentation.

## Documentation

- **Site** — [molcrafts.github.io/molpack](https://molcrafts.github.io/molpack/) — Home, Get Started, Guide, Internals
- **Rust API** — [docs.rs/molcrafts-molpack](https://docs.rs/molcrafts-molpack)
- **Python** — [`python/docs/`](./python/docs/)
- **Examples** — runnable programs in [`examples/`](./examples/) and [`python/examples/`](./python/examples/)

## MolCrafts ecosystem

| Project | Role |
|---------|------|
| [molpy](https://github.com/MolCrafts/molpy)     | Python toolkit — the shared molecular data model & workflow layer |
| [molrs](https://github.com/MolCrafts/molrs)     | Rust core — molecular data structures & compute kernels (native + WASM) |
| **molpack** | Packmol-grade molecular packing (Rust + Python) — this repo |
| [molvis](https://github.com/MolCrafts/molvis)   | WebGL molecular visualization & editing |
| [molexp](https://github.com/MolCrafts/molexp)   | Workflow & experiment-management platform |
| [molnex](https://github.com/MolCrafts/molnex)   | Molecular machine-learning framework |
| [molq](https://github.com/MolCrafts/molq)       | Unified job queue — local / SLURM / PBS / LSF |
| [molcfg](https://github.com/MolCrafts/molcfg)   | Layered configuration library |
| [mollog](https://github.com/MolCrafts/mollog)   | Structured logging, stdlib-compatible |
| [molhub](https://github.com/MolCrafts/molhub)   | Molecular dataset hub |
| [molmcp](https://github.com/MolCrafts/molmcp)   | MCP server for the ecosystem |
| [molrec](https://github.com/MolCrafts/molrec)   | Atomistic record specification |

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Bugs and feature requests via
[GitHub Issues](https://github.com/MolCrafts/molpack/issues).

## License

BSD-3-Clause — see [LICENSE](LICENSE).

## References

- Martínez, L.; Andrade, R.; Birgin, E. G.; Martínez, J. M.
  **PACKMOL: A package for building initial configurations for molecular dynamics simulations.**
  *J. Comput. Chem.* **2009**, *30* (13), 2157–2164.
  https://doi.org/10.1002/jcc.21224

<hr>

<div align="center">
<sub>Crafted with 💚 by <a href="https://github.com/MolCrafts">MolCrafts</a></sub>
</div>
