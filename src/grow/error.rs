//! Why a target cannot be grown.
//!
//! This leaf names growth's refusals. It imports neither [`PackError`](crate::PackError)
//! nor a stage, so the crate error can wrap it and the stages can return it
//! without a cycle. Lattice-only variants live here for the same reason:
//! one `GrowError` is what every growth entry maps, and the lattice stage is
//! their only constructor.

use std::fmt;

use molrs::op::F;

/// Why a target cannot be grown.
///
/// Growth consumes the template's *chemistry* (its bond graph), so a target
/// that carries none is refused with a named error — never silently degraded
/// to rigid-body packing. The caller chooses the method per target.
///
/// Template-graph refusals stop at the first match, in this order:
/// `NoAtomsBlock → MissingBondEndpoint → BondOutOfRange → NoBonds →
/// TemplateTooSmall → Disconnected → RingTemplate`.
#[derive(Debug, Clone)]
pub enum GrowError {
    /// The target was built without a template frame
    /// ([`Target::from_coords`][crate::Target::from_coords]), so there is no
    /// bond graph to grow from.
    MissingTemplate,
    /// The template has fewer than 3 atoms; growth needs a rigid seed of 3.
    TemplateTooSmall(usize),
    /// The template frame has no readable `atoms` block (`x` / `y` / `z` in Å).
    NoAtomsBlock,
    /// A non-empty bonds block has no `atomi` or `atomj` column.
    /// `column` is the missing endpoint name.
    MissingBondEndpoint {
        /// The absent endpoint column (`atomi` or `atomj`).
        column: &'static str,
    },
    /// A bond row names an atom outside the template. The fields are
    /// [`molrs::core::TopologyError::EndpointOutOfRange`].
    BondOutOfRange {
        /// 0-based row of the bonds block.
        row: usize,
        /// The out-of-range endpoint.
        atom: usize,
        /// Atom count of the template (`atoms` rows).
        n: usize,
    },
    /// The template frame carries no bonds (missing/empty graph).
    NoBonds,
    /// The template's bond graph does not connect all atoms.
    Disconnected,
    /// No placed reference atom could be found while decomposing atom `.0`.
    NoReference(usize),
    /// Rotatable-bond perception failed.
    Perceive(String),
    /// A Grow target needs a box: neither a periodic box nor a cell (nor a
    /// density, once `with_density` lands) was declared.
    NoBox,
    /// The declared cell is not orthorhombic; the v1 overlap field only
    /// supports orthorhombic boxes (`triclinic-cell-downshift` lifts this).
    TriclinicCell,
    /// `fixed_at` combined with a growth entry — a fixed placement is by
    /// definition not grown.
    FixedTarget,
    /// The template's bond graph contains a cycle. Growth decomposes the
    /// template into a *tree* of internal coordinates; a ring bond would be
    /// silently dropped and the ring grown open — refused instead (lattice
    /// ring closure is its own future spec).
    RingTemplate,
    /// A special-bonds weight is neither 0 nor 1. Growth compiles a binary
    /// skip table; `index` is the 0-based slot (slot 0 is 1-2, slot 2 is 1-4).
    NonBinarySpecialBond {
        /// 0-based table slot.
        index: usize,
        /// The stored weight at `index`.
        weight: F,
    },
    /// The template's heavy-atom backbone does not fit the diamond lattice
    /// (the message names the offense). Tetrahedral heavy degree `1..=4` is
    /// accepted (linear is the `d = 2` degeneracy); degree `> 4`, a
    /// detached heavy, or a non-sp³ interior bond is named here.
    NonTetrahedralTemplate(String),
    /// The attached region contains no usable diamond site for this chain
    /// (empty Region ∩ lattice, or the tree cannot embed in the allowed
    /// subgraph).
    LatticeRegionEmpty,
}

impl fmt::Display for GrowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GrowError::MissingTemplate => write!(
                f,
                "the target has no template frame (built from bare coordinates); growth needs \
                 a bond graph — load the species from a file or build a frame with bonds, or \
                 pack this target with GenCanPack"
            ),
            GrowError::TemplateTooSmall(n) => write!(
                f,
                "the template has {n} atom(s); growth needs at least 3 — pack this \
                 target with GenCanPack"
            ),
            GrowError::NoAtomsBlock => {
                write!(f, "the template frame has no readable atoms block")
            }
            GrowError::MissingBondEndpoint { column } => {
                write!(f, "the template's bonds block has no '{column}' column")
            }
            GrowError::BondOutOfRange { row, atom, n } => write!(
                f,
                "bond row {row} references atom {atom} outside the template (natoms = {n})"
            ),
            GrowError::NoBonds => write!(
                f,
                "the template frame carries no bonds; growth needs the bond graph — pack \
                 this target with GenCanPack or supply connectivity"
            ),
            GrowError::Disconnected => {
                write!(f, "the template's bond graph does not connect all atoms")
            }
            GrowError::NoReference(i) => write!(
                f,
                "no placed reference atom found while decomposing atom {i}"
            ),
            GrowError::Perceive(msg) => write!(f, "rotatable-bond perception failed: {msg}"),
            GrowError::NoBox => write!(
                f,
                "growth needs a box: declare with_periodic_box / with_cell (or a density) \
                 — the box is at its final volume from the first atom"
            ),
            GrowError::TriclinicCell => write!(
                f,
                "growth currently supports orthorhombic boxes only; declare an \
                 orthorhombic cell or periodic box"
            ),
            GrowError::RingTemplate => write!(
                f,
                "the template contains a ring: growth decomposes the bond graph \
                 into a tree, and a ring bond would be silently dropped; ring \
                 templates are refused"
            ),
            GrowError::NonBinarySpecialBond { index, weight } => {
                let pair = index + 2;
                write!(
                    f,
                    "special-bonds 1-{pair} weight is {weight}; growth compiles a \
                     binary skip table (0 or 1). For all-atom explicit hydrogen, \
                     use Target::with_atom_radius rather than a fractional weight"
                )
            }
            GrowError::NonTetrahedralTemplate(msg) => write!(
                f,
                "the template's backbone does not fit the diamond lattice: {msg}"
            ),
            GrowError::FixedTarget => write!(
                f,
                "a fixed target cannot be grown: drop fixed_at or pack this \
                 target with GenCanPack"
            ),
            GrowError::LatticeRegionEmpty => write!(
                f,
                "LatticeGrow: the attached region contains no usable diamond \
                 site for this chain; enlarge the mesh or reduce the template"
            ),
        }
    }
}

impl std::error::Error for GrowError {}
