//! Script loader: parse molpack's `.inp` input format and turn it into
//! a configured [`Molpack`](crate::Molpack) plus a list of
//! [`Target`](crate::Target)s.
//!
//! Typical use from a frontend (CLI, PyO3 binding, embedding host):
//!
//! ```no_run
//! use std::path::Path;
//! use molpack::script;
//!
//! let src = std::fs::read_to_string("mixture.inp")?;
//! let script = script::parse(&src)?;
//! let built = script.build(Path::new("."))?;
//!
//! let result = built.packer.pack(&built.targets, built.nloop)?;
//! script::write_frame(&built.output, &result.frame)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! The module intentionally keeps parsing, lowering, and I/O separate
//! so embedders can intercept any stage — e.g. mutate the parsed
//! [`Script`] before `build`, attach a custom
//! [`Handler`](crate::Handler) to the packer, or route output through
//! a different writer.

mod build;
mod error;
mod io;
mod parser;

pub use build::BuildResult;
pub use error::ScriptError;
pub use io::{read_frame, write_frame};
pub use parser::{AtomGroup, PbcSpec, RestraintSpec, Script, Structure, parse};
