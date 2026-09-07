# Regression — lattice-branch-saw-02-docs

Public API / caller contract. 2026-09-05. No third-party runtime.

- Example pick: `LatticeGrow` then `GenCanPack.with_restart` at **2.0 Å** (not Auhl, not `CbmcGrow` fallback).
- Trees including branched are accepted; cycles raise `RingTemplate`.
- Forbidden names: `StarGrow`, `pack_star_lattice`, C3/BFS as public contract, linear-only/v1 on `LatticeGrow`, `packmol` in identifiers.
- Tiny CLI: `python python/examples/pack_peo_topo.py star 2 1 0.2 1`
- 4-arm star `#[X4](#[EO]:n):4` is a tree (`n_bonds = n_atoms - 1`); `arm_length=2` → `n_atoms = 77`.
