//! Script loader: parse molpack's `.inp` input format and turn it into
//! a configured [`Molpack`](crate::Molpack) plus a list of
//! [`Target`](crate::Target)s.
//!
//! Two front-end shapes are supported:
//!
//! - **Native (feature `io`)** — [`Script::build`] reads template files
//!   via molrs-io and returns a ready-to-run [`BuildResult`]:
//!
//!   ```ignore
//!   use std::path::Path;
//!   use molpack::script;
//!
//!   let src = std::fs::read_to_string("mixture.inp")?;
//!   let script = script::parse(&src)?;
//!   let built = script.build(Path::new("."))?;
//!
//!   let result = built.packer.pack(&built.targets, built.nloop)?;
//!   script::write_frame(&built.output, &result.frame)?;
//!   # Ok::<(), Box<dyn std::error::Error>>(())
//!   ```
//!
//! - **Embedding hosts (any feature set)** — [`Script::lower`] returns
//!   a [`ScriptPlan`] with file paths resolved but unread. The caller
//!   loads each [`StructurePlan::filepath`] with its own frame loader,
//!   builds a [`Target`], and stamps restraints via
//!   [`StructurePlan::apply`]. This is what the PyO3 wheel uses, so it
//!   does not have to statically link molrs-io.
//!
//! Parsing, lowering, and I/O are kept separate so embedders can
//! intercept any stage — e.g. mutate the parsed [`Script`] before
//! `lower`, attach a custom [`Handler`](crate::Handler) to the packer,
//! or route output through a different writer.

mod build;
mod error;
#[cfg(feature = "io")]
mod io;
mod parser;

#[cfg(feature = "io")]
pub use build::BuildResult;
pub use build::{ScriptPlan, StructurePlan};
pub use error::ScriptError;
#[cfg(feature = "io")]
pub use io::{read_frame, write_frame};
pub use parser::{AtomGroup, PbcSpec, RestraintSpec, Script, Structure, parse};
