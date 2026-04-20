//! Integration tests for the `molpack` CLI binary.
//!
//! These tests compile and invoke the binary directly so they require the
//! `cli` feature.  Run with:
//!
//! ```bash
//! cargo test --test cli --features cli
//! ```

use std::path::PathBuf;
use std::process::Command;

fn bin_path() -> PathBuf {
    // `cargo test` puts the test binary next to the CLI binary in the same
    // target directory profile sub-dir.  CARGO_BIN_EXE_molpack is injected
    // by cargo when the test crate declares an [[integration-test]] against a
    // binary that shares the same manifest.
    PathBuf::from(env!("CARGO_BIN_EXE_molpack"))
}

fn example_dir(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join(name)
}

// ── parser error cases ───────────────────────────────────────────────────────

#[test]
fn exits_nonzero_on_missing_output_keyword() {
    let out = Command::new(bin_path())
        .args(["--"])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .expect("spawn")
        .wait_with_output()
        .expect("wait");
    assert!(!out.status.success());
}

#[test]
fn exits_nonzero_on_nonexistent_input_file() {
    let out = Command::new(bin_path())
        .arg("does_not_exist.inp")
        .output()
        .expect("run");
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("cannot read"),
        "expected 'cannot read' in stderr, got: {stderr}"
    );
}

#[test]
fn exits_nonzero_on_parse_error() {
    // An .inp with an unclosed structure block.
    let bad_inp = "tolerance 2.0\noutput x.pdb\nstructure mol.pdb\n  number 1\n";
    let out = Command::new(bin_path())
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .and_then(|mut child| {
            use std::io::Write;
            child
                .stdin
                .as_mut()
                .unwrap()
                .write_all(bad_inp.as_bytes())?;
            child.wait_with_output()
        })
        .expect("run");
    assert!(!out.status.success());
}

// ── smoke tests (run actual packing on canonical examples) ──────────────────

/// Pack the mixture example and verify the output file is created.
#[test]
fn smoke_pack_mixture() {
    let dir = example_dir("pack_mixture");
    let out_path = dir.join("_ci_mixture.pdb");
    let _ = std::fs::remove_file(&out_path);

    // Write a temp .inp that redirects output so we don't touch the committed file.
    let inp = format!(
        "tolerance 2.0\nseed 1234567\nfiletype pdb\noutput {}\n\n\
         structure water.pdb\n  number 50\n  inside box 0. 0. 0. 20. 20. 20.\nend structure\n",
        out_path.display()
    );

    let out = Command::new(bin_path())
        .current_dir(&dir)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .and_then(|mut child| {
            use std::io::Write;
            child.stdin.as_mut().unwrap().write_all(inp.as_bytes())?;
            child.wait_with_output()
        })
        .expect("run molpack");

    assert!(
        out.status.success(),
        "molpack exited {:?}\nstdout: {}\nstderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        out_path.exists(),
        "output file not created: {}",
        out_path.display()
    );

    // Clean up.
    let _ = std::fs::remove_file(&out_path);
}

/// Verify file-argument mode resolves paths relative to the .inp directory.
#[test]
fn file_arg_resolves_paths_from_inp_dir() {
    let dir = example_dir("pack_mixture");
    let out_path = dir.join("_ci_mixture_filearg.pdb");
    let _ = std::fs::remove_file(&out_path);

    let inp = format!(
        "tolerance 2.0\nseed 42\nfiletype pdb\noutput {}\n\n\
         structure water.pdb\n  number 30\n  inside box 0. 0. 0. 15. 15. 15.\nend structure\n",
        out_path.display()
    );

    // Write a temp .inp file next to the PDB files.
    let inp_path = dir.join("_ci_test.inp");
    std::fs::write(&inp_path, &inp).expect("write inp");

    let out = Command::new(bin_path())
        .arg(&inp_path)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .output()
        .expect("run molpack");

    let _ = std::fs::remove_file(&inp_path);
    let _ = std::fs::remove_file(&out_path);

    assert!(
        out.status.success(),
        "molpack exited {:?}\nstdout: {}\nstderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}
