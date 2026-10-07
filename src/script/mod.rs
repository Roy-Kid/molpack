//! Script loader: parse molpack's `.inp` input format and turn it into
//! a configured [`GencanPack`](crate::GencanPack) plus a list of
//! [`Target`](crate::Target)s.
//!
//! Two front-end shapes are supported:
//!
//! - **Native (feature `io`)** — `Script::build` reads each template with
//!   the molrs reader of its [`StructureFormat`] (the script's `filetype`,
//!   else the file name) and returns a ready-to-run `ScriptJob`; the
//!   packed frame goes out through `StructureFormat::write`. These names
//!   are compiled only when the `io` feature is on, so they are written in
//!   plain code font here rather than as cross-references:
//!
//!   ```no_run
//!   use std::path::Path;
//!   use molpack::{PackEngine, script};
//!   use molpack::script::StructureFormat;
//!
//!   let src = std::fs::read_to_string("mixture.inp")?;
//!   let script = script::parse(&src)?;
//!   let job = script.build(Path::new("."))?;
//!
//!   let state = job.packer.run(&job.targets, job.nloop)?;
//!   StructureFormat::resolve(&job.output, None)?.write(&job.output, &state.frame)?;
//!   # Ok::<(), Box<dyn std::error::Error>>(())
//!   ```
//!
//! - **Embedding hosts (any feature set)** — [`Script::lower`] returns
//!   a [`ScriptPlan`] with file paths resolved but unread. The caller
//!   loads each [`StructurePlan::filepath`] with its own frame loader
//!   (resolving its format with [`StructureFormat::resolve`]),
//!   builds a [`Target`](crate::Target), and stamps restraints via
//!   [`StructurePlan::apply`]. This is what the PyO3 wheel uses, so it
//!   does not have to statically link molrs-io.
//!
//! Parsing, lowering, and I/O are kept separate so embedders can
//! intercept any stage — e.g. mutate the parsed [`Script`] before
//! `lower`, attach a custom [`Callback`](crate::Callback) to the packer,
//! or route output through a different writer.

mod build;
mod error;
mod parser;
mod structure_format;

#[cfg(feature = "io")]
pub use build::ScriptJob;
pub use build::{ScriptPlan, StructurePlan};
pub use error::ScriptError;
pub use parser::{AtomGroup, PbcSpec, RestraintSpec, Script, Structure, parse};
pub use structure_format::StructureFormat;
