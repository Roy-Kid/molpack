//! Context layer for packmol-aligned packing runtime.
//!
//! The run's types — [`PackContext`], [`PackState`], [`Placed`],
//! [`RigidView`] — are published at the crate root; this
//! module publishes only the layout details a custom objective or handler
//! reads off a `PackContext` (per-atom property flags, scratch buffers, the
//! geometry cache key). The leaves are private: one path per item.

pub(crate) mod build;
mod geometry;
pub(crate) mod grid;
pub(crate) mod pack_context;
// `Stage::run` takes `&mut PackState` and the two declaration methods speak in
// `Placed`, so both types are part of the published contract (the crate root
// re-exports them). The module's free `evaluate_unscaled` stays crate-private
// — it is a shared primitive, not a promise.
pub(crate) mod pack_state;
pub(crate) mod rigid_view;
mod work_buffers;

pub use geometry::GeometryKey;
pub use pack_context::{ATOM_FLAG_FIXED, ATOM_FLAG_SHORT, AtomProps, NONE_IDX};
pub(crate) use pack_context::{DEFAULT_SCALE2, PackContext};
pub(crate) use pack_state::{PackState, Placed};
pub(crate) use rigid_view::RigidView;
pub use work_buffers::WorkBuffers;
