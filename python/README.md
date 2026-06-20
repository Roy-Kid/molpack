# molcrafts-molpack

Python bindings for **molpack** — a faithful Rust port of Packmol for
molecular packing, with a Packmol-compatible `.inp` script format.

The wheel is built without the `io` feature: structure frames are loaded
through your installed `molrs` Python package, and molpack builds packing
targets from the loaded frame.

See the [project repository](https://github.com/MolCrafts/molpack) for full
documentation, the `.inp` script reference, and examples.
