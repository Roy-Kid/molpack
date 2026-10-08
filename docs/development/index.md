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
  `Region`, and `Callback` implementations, plus binding a custom in-loop
  optimizer.

## Validation commands

```bash
cargo test -p molcrafts-molpack --lib --features cli
cargo test -p molcrafts-molpack --doc --features cli
cd python
maturin develop --release
pytest
cargo fmt
cargo clippy --all-targets --all-features -- -D warnings
```

Behaviour is tested in `#[cfg(test)]` modules next to the code that owns it;
molpack has no `tests/` directory, no benchmark suite and no regression
harness. The whole Rust tier runs in seconds.
