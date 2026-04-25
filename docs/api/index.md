# API Reference

This chapter is a navigation page, not the final contract itself. The user guide
teaches the model and the workflows. The generated Rust `rustdoc` remains the
authoritative definition of the public API.

## What You'll Learn

Use this page when you need to:

- jump from a concept chapter to the exact Rust type or trait;
- find the main entry points of the public API quickly;
- decide whether you need the narrative guide or the raw API contract.

## 1. Primary Entry Points

- [Crate root](https://docs.rs/molcrafts-molpack/latest/molpack/) for the full
  public surface and module tree.
- [Molpack](https://docs.rs/molcrafts-molpack/latest/molpack/struct.Molpack.html)
  for builder-level configuration and `pack()`.
- [Target](https://docs.rs/molcrafts-molpack/latest/molpack/struct.Target.html)
  for molecule-type specification.
- [PackResult](https://docs.rs/molcrafts-molpack/latest/molpack/struct.PackResult.html)
  for result inspection.

## 2. Extension Traits and Advanced Hooks

- [Restraint](https://docs.rs/molcrafts-molpack/latest/molpack/trait.Restraint.html)
  for custom soft geometry.
- [Region](https://docs.rs/molcrafts-molpack/latest/molpack/trait.Region.html)
  and
  [RegionExt](https://docs.rs/molcrafts-molpack/latest/molpack/trait.RegionExt.html)
  for compositional geometry.
- [Handler](https://docs.rs/molcrafts-molpack/latest/molpack/trait.Handler.html)
  for progress observation and controlled early stop.
- [Relaxer](https://docs.rs/molcrafts-molpack/latest/molpack/trait.Relaxer.html)
  and
  [RelaxerRunner](https://docs.rs/molcrafts-molpack/latest/molpack/trait.RelaxerRunner.html)
  for in-loop conformational updates.

## 3. Narrative Chapters Inside `rustdoc`

The crate also publishes long-form Rust chapters as `rustdoc` modules. Those
chapters are the bridge between this site and the raw API:

- [getting_started](https://docs.rs/molcrafts-molpack/latest/molpack/getting_started/index.html)
- [concepts](https://docs.rs/molcrafts-molpack/latest/molpack/concepts/index.html)
- [architecture](https://docs.rs/molcrafts-molpack/latest/molpack/architecture/index.html)
- [extending](https://docs.rs/molcrafts-molpack/latest/molpack/extending/index.html)

## 4. Work Locally When You Are Changing the Crate

Use local `rustdoc` when you are developing on the crate itself:

```bash
cargo doc --open
```

That is especially useful for checking intra-crate links and reading unpublished
changes before they reach docs.rs.

## 5. Python Note

The Python binding intentionally mirrors the Rust mental model, but this API
section remains Rust-first. On this site, Python coverage is merged into the
workflow pages that matter most: [Install](../install.md), [Getting Started](../getting_started.md),
and [Examples](../examples.md).

## Takeaway

Use this page as the bridge between the guide and the exact API contract. When
you need precision, go to `rustdoc`; when you need explanation, stay in the
narrative chapters.
