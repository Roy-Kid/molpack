//! Context layer for packmol-aligned packing runtime.

pub(crate) mod build;
pub mod model;
pub mod pack_context;
// Crate-private for the duration of the stage-pipeline chain: `PackState`
// and `Placed` become part of the public seam with the stage signature, not
// before (`src/lib.rs` gains nothing here).
pub(crate) mod pack_state;
pub mod rigid_view;
pub mod state;
pub mod work_buffers;

pub use model::ModelData;
pub use pack_context::{ATOM_FLAG_FIXED, ATOM_FLAG_SHORT, AtomProps, NONE_IDX, PackContext};
pub use rigid_view::RigidView;
pub use state::{RuntimeState, RuntimeStateMut};
pub use work_buffers::WorkBuffers;
