# Formats

The CLI reads each molecule template with the molrs reader of its format
(`molrs::io::read_pdb`, `read_xyz`, …) and writes the final packed structure
with the molrs writer of the output's format to the path named by the
script's `output` keyword. molrs names every door after its format and picks
none for the caller, so the choice below — `filetype`, else the file name —
is the `.inp` grammar's, held once in `molpack::script::StructureFormat`;
reading and writing each format is molrs's.

The output format is inferred from the output file name. Input formats can be
inferred from structure-file names or set globally with `filetype`, which
accepts a format name or any of its extensions.

| Format | Read | Write | File names / `filetype` |
|---|---:|---:|---|
| PDB | Yes | Yes | `.pdb`, `.ent`; `pdb` |
| XYZ / extended XYZ | Yes | Yes | `.xyz`, `.extxyz`; `xyz` |
| SDF / MOL | Yes | No | `.sdf`, `.mol`; `sdf` |
| MOL2 | Yes | Yes | `.mol2`; `mol2` |
| GROMACS GRO | Yes | Yes | `.gro`; `gro` |
| CIF | Yes | Yes | `.cif`; `cif` |
| VASP POSCAR | Yes | Yes | `.poscar`, `.vasp`, `POSCAR*`, `CONTCAR*`; `poscar` |
| XSF | Yes | Yes | `.xsf`; `xsf` |
| Gaussian cube | Yes | Yes | `.cube`, `.cub`; `cube` |
| AMBER inpcrd | Yes | No | `.inpcrd`, `.rst7`, `.restrt`, `.crd`; `inpcrd` |
| LAMMPS data | Yes | Yes | `.data`, `.lmp`; `lammps_data` |
| LAMMPS dump | Yes (first snapshot) | Yes | `.lammpstrj`, `.dump`; `lammps_dump` |

## Example

```text
filetype pdb
output packed.xyz

structure water.pdb
  number 100
  inside box 0. 0. 0. 30. 30. 30.
end structure
```

The input template is read as PDB because of `filetype pdb`; the output is
written as XYZ because the output path ends with `.xyz`.
