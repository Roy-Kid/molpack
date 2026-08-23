//! ABI-line pairing across the two independently compiled units — and the
//! reference that keeps the path-only `molrs-ffi` dev-dependency actually
//! linked into molpack's test binaries (an unreferenced `--extern` is never
//! loaded by rustc, which would silently drop these binaries out of the
//! shared-dylib circle).
//!
//! The comparison is not configuration self-certification: `molrs::VERSION`
//! is baked in by molpack's own molrs dependency, while
//! `molrs_ffi::abi::abi_line()` is baked into the separately compiled
//! `libmolrs_ffi`. Equality is what makes a capsule crossing between them
//! legal (minor line == ABI version, see molrs `docs/interop.md`), so it is
//! asserted rather than assumed — and never hard-coded.

#[test]
fn abi_line_matches_the_molrs_minor_line() {
    let expected = molrs::VERSION
        .rsplit_once('.')
        .map(|(major_minor, _patch)| major_minor)
        .expect("molrs::VERSION is major.minor.patch");

    assert_eq!(molrs_ffi::abi::abi_line(), expected);
}
