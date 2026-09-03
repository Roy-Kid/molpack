//! Context layer for packmol-aligned packing runtime.

pub(crate) mod build;
pub mod model;
pub mod pack_context;
// Public since the stage seam: `Stage::run` takes `&mut PackState` and the
// two declaration methods speak in `Placed`, so both types are part of the
// published contract (`src/lib.rs` re-exports them). The module's free
// `evaluate_unscaled` stays crate-private — it is a shared primitive, not a
// promise.
pub mod pack_state;
pub mod rigid_view;
pub mod state;
pub mod work_buffers;

pub use model::ModelData;
pub use pack_context::{ATOM_FLAG_FIXED, ATOM_FLAG_SHORT, AtomProps, NONE_IDX, PackContext};
pub use pack_state::{PackState, Placed};
pub use rigid_view::RigidView;
pub use state::{RuntimeState, RuntimeStateMut};
pub use work_buffers::WorkBuffers;
