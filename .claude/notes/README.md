# .claude/notes/ — passive project knowledge

Internal context that agents read and that outlives any single feature.
Nothing here is public documentation (that is `docs/`) and nothing here
is an active spec (that is `.claude/specs/`).

| File | What it is | Who writes it |
|---|---|---|
| `law.md` | The rulebook: eleven engineering laws plus molpack's project invariants (§ IX). Outranks every other file here. CLAUDE.md indexes one line per law. | operator via `/mol:note`; `/mol:bootstrap` appends missing default ids, never rewrites |
| `conventions.md` | Facts agents need but that are not laws: Cargo features, coding style, test tiers and gates, repo layout, the molrs sibling checkout, build cache, ABI gates. | `/mol:note`, `/mol:compact` |
| `architecture.md` | The blueprint: modules, public surface, style, layer roles. `librarian` reads it during `/mol:spec`; `architect` enforces against it. | `/mol:map` (managed block); custom annotations outside the markers |
| `notes.md` | Evolving decisions with their reasons and history. Rules that become absolute move to `law.md`. | `/mol:note` |

Navigation: start at `CLAUDE.md` (router) → `law.md` (what must never happen)
→ `conventions.md` (how things are built and tested) → `architecture.md`
(where things are) → `notes.md` (why things are the way they are).
