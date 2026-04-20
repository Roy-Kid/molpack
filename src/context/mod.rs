//! Context layer for packmol-aligned packing runtime.

pub mod model;
pub mod pack_context;
pub mod state;
pub mod work_buffers;

pub use model::ModelData;
pub use pack_context::{ATOM_FLAG_FIXED, ATOM_FLAG_SHORT, AtomProps, NONE_IDX, PackContext};
pub use state::{RuntimeState, RuntimeStateMut};
pub use work_buffers::WorkBuffers;
