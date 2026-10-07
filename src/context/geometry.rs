//! Identity of the packing geometry for the evaluation cache.
//!
//! The packing context builds a key and the work buffers store it. This
//! module imports neither, so those two do not import each other for the key.

use molrs::op::F;

/// Identity of the packing geometry — the cell partition plus the lattice it
/// partitions. Compared by the evaluation cache to decide whether a previous
/// cell assignment is still valid.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct GeometryKey {
    pub celldim: [u32; 3],
    pub pbc: [bool; 3],
    pub h: [F; 9],
    pub origin: [F; 3],
}
