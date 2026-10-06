# Examples

The Rust examples live under `examples/pack_*` and use the `io` feature to read
the bundled structure files.

```bash
cargo run --release --features io --example pack_mixture
cargo run --release --features io --example pack_bilayer
cargo run --release --features io --example pack_interface
cargo run --release --features io --example pack_spherical
cargo run --release --features io --example pack_solvprotein
```

| Example | Demonstrates |
|---|---|
| `pack_mixture` | Two-component water/urea cube. |
| `pack_bilayer` | Plane restraints and atom-subset orientation. |
| `pack_interface` | Solvents around a fixed molecule. |
| `pack_spherical` | Nested spherical shells. |
| `pack_solvprotein` | Fixed protein solvation with water and ions. |

Optional progress and trajectory output are controlled by environment
variables:

```bash
MOLPACK_EXAMPLE_PROGRESS=1 cargo run --release --features io --example pack_mixture
MOLPACK_EXAMPLE_XYZ=1 cargo run --release --features io --example pack_mixture
```

For Packmol-style input files, use the matching [CLI examples](../cli/examples.md).
