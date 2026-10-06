# Packmol Script Notes

The CLI reads Packmol-style `.inp` files. For the script language itself, use
the official [Packmol user guide](https://m3g.github.io/packmol/userguide.shtml)
as the reference. molpack intentionally keeps that job shape: global settings,
then one `structure ... end structure` block per molecule type.

```text
tolerance 2.0
filetype pdb
output packed.pdb

structure water.pdb
  number 1000
  inside box 0. 0. 0. 40. 40. 40.
end structure
```

This page only documents molpack-specific behavior and extensions. If a Packmol
keyword is not listed here, treat Packmol's manual as the source of truth and
use [Packmol parity](../packmol_parity.md) for implementation-level
compatibility notes.

## Invocation and paths

- `molpack mixture.inp`: paths inside `mixture.inp` resolve relative to the
  directory that contains `mixture.inp`.
- `molpack < mixture.inp`: paths resolve relative to the current working
  directory.

This matches the two common Packmol invocation styles while making file-argument
mode safe to run from another directory.

## molpack additions

| Keyword | Meaning |
|---|---|
| `avoid_overlap <yes\|no\|true\|false\|1\|0>` | Controls whether initial random placements are rejected when they overlap a fixed molecule. Default: on. Use `avoid_overlap no` only when you explicitly want Packmol's fixed-solute exclusion guard disabled. |
| `filetype <format>` | Any format `molrs::io::read_frame` reads — `sdf`, `mol2`, `gro`, `cif`, `lammps_data`, `lammps_dump`, … — named or by extension. See [Formats](formats.md). |

Output format is inferred from the `output` file name, through
`molrs::io::write_frame` (see [Formats](formats.md)).

## Strict parsing

molpack rejects unknown keywords instead of silently ignoring them. A misspelled
top-level keyword returns `ScriptError::UnknownKeyword`, which prevents broken
inputs from falling through to an expensive or invalid packing run.

## Supported Packmol subset

The parser covers the Packmol syntax used by the bundled examples and
regression suite:

- global `tolerance`, `seed`, `filetype`, `output`, `nloop`, and `pbc`;
- `structure ... end structure` blocks;
- `number`, `center`, `centerofmass`, and `fixed`;
- whole-molecule and `atoms ... end atoms` restraints;
- `inside`/`outside` box, cube, sphere, ellipsoid, and cylinder;
- `over plane`, `above plane`, and `below plane`.

Atom indices in `.inp` files are 1-based, matching Packmol. When translating
the same model to Python or Rust, subtract 1 because those APIs use 0-based
indices.
