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

- [ ] `cargo test -p molcrafts-molpack --lib --features cli,ff` passes
- [ ] `cargo test -p molcrafts-molpack --doc --features cli,ff` passes
- [ ] `cargo clippy --all-targets --all-features -- -D warnings` clean
- [ ] `cargo fmt --check` clean
- [ ] New behaviour has a unit test in the module that owns it
- [ ] Python tests pass (`cd python && pytest`)

## Breaking changes

<!-- Does this change any public API? If yes, describe what callers need to update. -->

None / <!-- describe -->

## Notes for reviewer

<!-- Anything that needs special attention or context. -->
