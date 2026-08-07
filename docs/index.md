---
title: molpack
description: Packmol-grade molecular packing in pure Rust, with Python bindings.
hide:
  - navigation
  - toc
hero:
  kicker: molpack Manual
  title: molpack
  description: Packmol-compatible molecular packing in pure Rust — a CLI for `.inp` scripts, a Python package for notebooks and pipelines, and a native Rust crate.
  install:
    label: Install
    methods:
      - label: CLI
        command: cargo install molcrafts-molpack --features cli
      - label: Rust library
        command: cargo add molcrafts-molpack
      - label: Python
        command: pip install molcrafts-molpack
  badges:
    - img: https://img.shields.io/crates/v/molcrafts-molpack?color=c2410c&label=crates.io
      href: https://crates.io/crates/molcrafts-molpack
      alt: crates.io version
    - img: https://img.shields.io/pypi/v/molcrafts-molpack?color=d97706&label=PyPI
      href: https://pypi.org/project/molcrafts-molpack/
      alt: PyPI version
    - img: https://github.com/MolCrafts/molpack/actions/workflows/ci.yml/badge.svg
      href: https://github.com/MolCrafts/molpack/actions/workflows/ci.yml
      alt: CI status
    - img: https://img.shields.io/badge/license-BSD--3--Clause-blue.svg
      href: https://github.com/MolCrafts/molpack/blob/master/LICENSE
      alt: License BSD-3-Clause
  actions:
    - label: Get started
      href: getting_started/
      style: primary
    - label: Packmol script
      href: cli/
    - label: Python API
      href: python/
    - label: Rust API
      href: rust/
---

<h1 class="molcrafts-sr-only">molpack</h1>

<div class="molcrafts-manual-home" markdown>

<section class="molcrafts-manual-section molcrafts-manual-section--compact" markdown>

<div class="molcrafts-manual-section__header" markdown>

<span class="molcrafts-manual-eyebrow">Start here</span>

## Find the right page

This page is an index into the manual. Every surface lowers to the same model:
molecule templates, copy counts, geometric restraints, then one packing run.

</div>

<nav class="molcrafts-manual-index" aria-label="Manual entry points">
  <a href="install/">
    <span>01</span>
    <strong>Install</strong>
    <em>CLI binary, crates.io crate, or PyPI wheel — pick a surface and verify it loads.</em>
  </a>
  <a href="getting_started/">
    <span>02</span>
    <strong>Quickstart</strong>
    <em>Pack 100 waters in a 40&nbsp;Å cube and read convergence diagnostics.</em>
  </a>
  <a href="cli/">
    <span>03</span>
    <strong>Packmol script</strong>
    <em>Run Packmol-style `.inp` jobs with structure, number, inside, fixed, and pbc.</em>
  </a>
  <a href="python/">
    <span>04</span>
    <strong>Python API</strong>
    <em>Build targets from frames, pack in notebooks, keep results in memory.</em>
  </a>
  <a href="rust/">
    <span>05</span>
    <strong>Rust API</strong>
    <em>Native Target and Molpack builders inside your crate.</em>
  </a>
  <a href="rust/handlers-relaxers/">
    <span>06</span>
    <strong>Handlers</strong>
    <em>Observe steps, dump trajectories, or stop a run early.</em>
  </a>
</nav>

</section>

<section id="workflows" class="molcrafts-manual-section" markdown>

<div class="molcrafts-manual-section__header" markdown>

<span class="molcrafts-manual-eyebrow">Workflows</span>

## Same packing model, four entry points

All entry points lower to the same target / count / restraint model. Pick the
surface that matches how the rest of your workflow is written.

</div>

<div class="molcrafts-workflow-list" markdown>

<article markdown>

<div class="molcrafts-workflow-list__meta">Path 01 · Script</div>

### Run a `.inp` job

Use the CLI for reproducible packing jobs that already live as Packmol-style
input files. See the [CLI guide](cli/).

```text
structure water.pdb
  number 1000
  inside box 0. 0. 0. 40. 40. 40.
end structure
```

</article>

<article markdown>

<div class="molcrafts-workflow-list__meta">Path 02 · Python</div>

### Pack from a notebook or pipeline

Load or build a frame, create immutable targets, and pass the packed frame to
the writer or analysis code you already use. See the [Python API](python/).

```python
water = Target(frame, 100).with_restraint(
    InsideBoxRestraint([0, 0, 0], [40, 40, 40])
)
packed = Molpack().with_seed(42).pack([water])
```

</article>

<article markdown>

<div class="molcrafts-workflow-list__meta">Path 03 · Rust</div>

### Embed the engine in a crate

Use the native builder API for applications, services, and tests. See the
[Rust API](rust/).

```rust
let frame = Molpack::new()
    .with_seed(42)
    .pack(&[water], 200)?;
```

</article>

<article markdown>

<div class="molcrafts-workflow-list__meta">Path 04 · Handler</div>

### Observe or stop a run

Handlers receive structured events from the packing loop. See
[Handlers and relaxers](rust/handlers-relaxers/).

```rust
impl Handler for WatchFdist {
    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        eprintln!("fdist={}", info.fdist);
    }
}
```

</article>

</div>

</section>

<section id="capabilities" class="molcrafts-manual-section" markdown>

<div class="molcrafts-manual-section__header" markdown>

<span class="molcrafts-manual-eyebrow">Features</span>

## Core capabilities

Dense summary of what the engine and each surface cover.

</div>

<dl class="molcrafts-feature-matrix">
  <div>
    <dt>Targets</dt>
    <dd>Molecule templates, copy counts, fixed placements, centering, and rotation bounds.</dd>
  </div>
  <div>
    <dt>Restraints</dt>
    <dd>Geometric boxes and spheres, planes, and collective distribution-matching terms.</dd>
  </div>
  <div>
    <dt>Packer</dt>
    <dd>GENCAN outer loop, pair-distance objective, optional rayon evaluation, seed control.</dd>
  </div>
  <div>
    <dt>Scripts</dt>
    <dd>Packmol-compatible `.inp` keywords, path resolution, and extra filetypes (SDF, LAMMPS).</dd>
  </div>
  <div>
    <dt>Python</dt>
    <dd>Frame-in / frame-out packing for notebooks and pipelines via molrs.</dd>
  </div>
  <div>
    <dt>Rust</dt>
    <dd>Native builders, handlers, relaxers, and feature-gated IO / FF / rayon.</dd>
  </div>
</dl>

</section>

<section class="molcrafts-manual-section" markdown>

<div class="molcrafts-manual-section__header" markdown>

<span class="molcrafts-manual-eyebrow">Map</span>

## Documentation map

A compact mirror of the navigation tree for returning users.

</div>

<div class="molcrafts-doc-map">
  <section>
    <h3>Tutorial</h3>
    <p>Install, quickstart, concepts, examples, and Packmol parity notes.</p>
  </section>
  <section>
    <h3>Packmol Script</h3>
    <p>`.inp` overview, script notes, formats, and CLI examples.</p>
  </section>
  <section>
    <h3>Python</h3>
    <p>Installation, targets, restraints, packer options, PBC, and API reference.</p>
  </section>
  <section>
    <h3>Rust</h3>
    <p>Builders, restraints, PBC, handlers, relaxers, and example programs.</p>
  </section>
  <section>
    <h3>Development</h3>
    <p>Contributing, architecture, and extending the engine.</p>
  </section>
</div>

</section>

<section class="molcrafts-manual-section molcrafts-manual-section--stack" markdown>

<div class="molcrafts-manual-section__header" markdown>

<span class="molcrafts-manual-eyebrow">Surfaces</span>

## Jump to a surface

</div>

<div class="molcrafts-manual-grid molcrafts-manual-grid--cols-3">
  <a href="cli/">
    <strong>Packmol script</strong>
    <em>Reproducible `.inp` jobs with familiar keywords.</em>
  </a>
  <a href="python/">
    <strong>Python</strong>
    <em>Notebooks, pipelines, and frame-based packing.</em>
  </a>
  <a href="rust/">
    <strong>Rust</strong>
    <em>Native builders, handlers, and relaxers.</em>
  </a>
</div>

</section>

</div>
