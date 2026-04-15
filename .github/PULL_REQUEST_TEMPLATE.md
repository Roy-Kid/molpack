## Summary

<!-- What does this PR do? One paragraph is enough. -->

## Motivation

<!-- Link to the issue this fixes, or explain why the change is needed. -->

Fixes #

## Changes

<!-- Bullet list of the concrete changes: new types, API surface changes, behaviour differences. -->

-

## Test plan

<!-- How did you verify this? Check all that apply. -->

- [ ] `cargo test -p molcrafts-molpack` passes
- [ ] `cargo clippy -- -D warnings` clean
- [ ] `cargo fmt --check` clean
- [ ] New tests added for new behaviour
- [ ] `cargo test --release --test examples_batch -- --ignored` passes (required for restraint / objective changes)
- [ ] Python tests pass (`cd python && pytest`)

## Breaking changes

<!-- Does this change any public API? If yes, describe what callers need to update. -->

None / <!-- describe -->

## Notes for reviewer

<!-- Anything that needs special attention or context. -->
