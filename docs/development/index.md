# Development

These pages are for modifying molpack rather than just using it.

molpack is a Rust packing engine with three public surfaces:

- `molcrafts-molpack` Rust library (`lib` name: `molpack`)
- `molpack` CLI binary (`cli` feature)
- `molcrafts-molpack` Python wheel under `python/`

## Read first

- [Architecture](../architecture.md) maps modules, data flow, optimizer loops,
  and the objective-evaluation hot path.
- [Extending](../extending.md) walks through custom `AtomRestraint`,
  `Region`, and `Handler` implementations, plus binding a custom in-loop
  optimizer (feature `ff`).

## Validation commands

```bash
cargo test -p molcrafts-molpack --lib --tests
cargo test -p molcrafts-molpack --release --test examples_batch -- --ignored
cd python
maturin develop --release
pytest
cargo fmt
cargo clippy -- -D warnings
```

The ignored `examples_batch` test runs the canonical Packmol-equivalent
workloads and is intentionally slower than the fast unit/integration tier.
