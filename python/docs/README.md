# molpack Python docs

[Zensical](https://zensical.org/)-powered documentation site for the
`molpack` Python bindings.

## Build locally

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install zensical

cd python/docs
zensical serve   # live preview on http://localhost:8000
zensical build   # static site under site/
```

## Layout

```
python/docs/
├─ zensical.toml          # site config + navigation
├─ README.md              # this file
└─ docs/
   ├─ index.md
   ├─ installation.md
   ├─ getting-started.md
   ├─ examples.md
   ├─ api-reference.md
   └─ guide/
      ├─ targets.md
      ├─ restraints.md
      ├─ packer.md
      └─ periodic-boundaries.md
```

The Rust crate's own API docs live under `/docs/*.md` at repo root
(loaded into rustdoc via `include_str!`). These Python docs are
independent — they describe only the PyO3 binding surface.
