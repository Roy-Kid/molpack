//! Error type for the script loader.

use std::fmt;
use std::path::PathBuf;

use crate::error::PackError;

/// Errors produced by the script module — from parsing, file loading, or
/// the downstream `pack()` call invoked by a script-driven run.
#[derive(Debug)]
pub enum ScriptError {
    /// Parse error with line number and diagnostic message.
    Parse { line: usize, message: String },
    /// An unrecognised keyword appeared where a known one was expected.
    /// Unlike a generic parse error, this carries the offending token and
    /// the context block so error messages can suggest fixes.
    UnknownKeyword {
        line: usize,
        keyword: String,
        context: &'static str,
    },
    /// Script is missing a required `output` keyword.
    MissingOutput,
    /// Script contains zero `structure` blocks.
    NoStructures,
    /// A `structure` or `atoms` block was not closed with `end …`.
    UnclosedBlock(&'static str),
    /// Reading or writing a molecular frame failed.
    Io { path: PathBuf, message: String },
    /// The packer itself returned an error.
    Pack(PackError),
}

impl fmt::Display for ScriptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse { line, message } => write!(f, "line {line}: {message}"),
            Self::UnknownKeyword {
                line,
                keyword,
                context,
            } => write!(
                f,
                "line {line}: unknown keyword `{keyword}` in {context}. \
                 Silently dropping it risks wrong semantics (e.g. \
                 `pbc` that was ignored blew the cell grid); if this \
                 keyword should be accepted, add it to the parser",
            ),
            Self::MissingOutput => write!(f, "missing required `output` keyword"),
            Self::NoStructures => write!(f, "no `structure` blocks found in script"),
            Self::UnclosedBlock(kind) => write!(f, "unexpected EOF: unclosed `{kind}` block"),
            Self::Io { path, message } => write!(f, "{}: {message}", path.display()),
            Self::Pack(err) => write!(f, "packing failed: {err}"),
        }
    }
}

impl std::error::Error for ScriptError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Pack(err) => Some(err),
            _ => None,
        }
    }
}

impl From<PackError> for ScriptError {
    fn from(err: PackError) -> Self {
        Self::Pack(err)
    }
}
