//! Target builder for molecular packing.

use std::sync::Arc;

use crate::frame::frame_to_coords_and_elements;
use crate::relaxer::Relaxer;
use crate::restraint::Restraint;
use molrs::types::F;

/// Cartesian axis selector used in `Target::with_rotation_bound` and
/// other API surfaces that need to name an axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

/// Angular quantity stored internally as radians.
///
/// Constructors make the unit explicit at the call site:
/// `Angle::from_degrees(30.0)` vs `Angle::from_radians(FRAC_PI_6)`.
/// Implements `Copy` — pass by value, no `&`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angle(F);

impl Angle {
    /// Zero rotation.
    pub const ZERO: Self = Self(0.0);

    pub const fn from_radians(rad: F) -> Self {
        Self(rad)
    }

    pub fn from_degrees(deg: F) -> Self {
        Self(deg * (std::f64::consts::PI as F) / 180.0)
    }

    pub const fn radians(self) -> F {
        self.0
    }

    pub fn degrees(self) -> F {
        self.0 * 180.0 / (std::f64::consts::PI as F)
    }
}

/// Centering behavior for structure coordinates.
///
/// Packmol semantics:
/// - `Auto`: free molecules are centered; fixed molecules are not centered.
/// - `Center`: force centering.
/// - `Off`: keep input coordinates unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CenteringMode {
    #[default]
    Auto,
    Center,
    Off,
}

/// Fixed-molecule placement: translation + Euler orientation.
#[derive(Debug, Clone)]
pub struct Placement {
    /// Translation vector `[x, y, z]`.
    pub position: [F; 3],
    /// Euler rotations around x / y / z in the `eulerfixed` convention,
    /// stored as [`Angle`] triples.
    pub orientation: [Angle; 3],
}

/// Deprecated legacy name for [`Placement`]. Will be removed in the
/// next release.
#[deprecated(
    since = "0.2.0",
    note = "Renamed to `Placement` — the `Fixed` prefix is redundant now that the constructor name carries the semantic."
)]
pub type FixedPlacement = Placement;

/// Describes one type of molecule to be packed.
#[derive(Debug, Clone)]
pub struct Target {
    /// Input coordinates as provided by the source structure.
    pub input_coords: Vec<[F; 3]>,
    /// Flat list of atom positions — the centered reference coordinates.
    /// Shape: natoms × 3, stored as Vec<[F; 3]>.
    pub ref_coords: Vec<[F; 3]>,
    /// Van der Waals radii per atom.
    pub radii: Vec<F>,
    /// Element symbols per atom (e.g. `"C"`, `"O"`). Defaults to `"X"` if unknown.
    pub elements: Vec<String>,
    /// Number of copies to pack.
    pub count: usize,
    /// Optional name for logging.
    pub name: Option<String>,
    /// Restraints applied to every atom of every molecule copy.
    pub molecule_restraints: Vec<Arc<dyn Restraint>>,
    /// Per-atom-subset restraints: `(atom_indices_0_based, restraint)`.
    /// Each entry holds the 0-based atom indices (converted from Packmol's
    /// 1-based convention at registration time) and the restraint applied to them.
    pub atom_restraints: Vec<(Vec<usize>, Arc<dyn Restraint>)>,
    /// Optional structure-level limit for the perturbation heuristic
    /// (Packmol's `maxmove`).
    pub perturb_budget: Option<usize>,
    /// Centering policy.
    pub centering: CenteringMode,
    /// Rotation bounds in Euler variable order
    /// `[beta(y), gama(z), teta(x)]` as `(center, half_width)` [`Angle`] pairs.
    pub rotation_bound: [Option<(Angle, Angle)>; 3],
    /// If `Some`, this molecule is fixed (one copy, placed at the given location).
    pub fixed_at: Option<Placement>,
    /// Per-target in-loop relaxers (e.g. torsion MC). Called in order each iteration.
    pub relaxers: Vec<Box<dyn Relaxer>>,
}

impl Target {
    /// Create a new target from a `molrs::Frame` (read from PDB/XYZ) and a copy count.
    ///
    /// Positions are extracted from the `"atoms"` block (`"x"`, `"y"`, `"z"` columns)
    /// and automatically centered at the geometric center.
    /// VdW radii and element symbols are looked up from the `"element"` column.
    pub fn new(frame: molrs::Frame, count: usize) -> Self {
        let (positions, radii, elements) = frame_to_coords_and_elements(&frame);
        Self::from_parts(&positions, &radii, elements, count)
    }

    /// Create a new target directly from coordinate arrays.
    ///
    /// Useful for testing or when coordinates are already available.
    /// Stores both raw input coordinates and a geometrically centered reference copy.
    /// Effective usage follows [`CenteringMode::Auto`] unless overridden.
    pub fn from_coords(frame_positions: &[[F; 3]], radii: &[F], count: usize) -> Self {
        let n = frame_positions.len();
        Self::from_parts(frame_positions, radii, vec!["X".to_string(); n], count)
    }

    fn from_parts(
        frame_positions: &[[F; 3]],
        radii: &[F],
        elements: Vec<String>,
        count: usize,
    ) -> Self {
        assert_eq!(
            frame_positions.len(),
            radii.len(),
            "positions and radii must have the same length"
        );
        let input_coords = frame_positions.to_vec();
        let ref_coords = centered_coords(frame_positions);
        Self {
            input_coords,
            ref_coords,
            radii: radii.to_vec(),
            elements,
            count,
            name: None,
            molecule_restraints: Vec::new(),
            atom_restraints: Vec::new(),
            perturb_budget: None,
            centering: CenteringMode::Auto,
            rotation_bound: [None, None, None],
            fixed_at: None,
            relaxers: Vec::new(),
        }
    }

    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Attach a restraint applied to every atom of every molecule copy.
    pub fn with_restraint(mut self, r: impl Restraint + 'static) -> Self {
        self.molecule_restraints.push(Arc::new(r));
        self
    }

    /// Attach a restraint for selected atoms of every molecule copy.
    ///
    /// # Atom indexing
    ///
    /// Indices are **0-based**, matching Rust convention: atom `0` is
    /// the first atom in the PDB/XYZ file. For example, `&[0, 1, 2]`
    /// selects the first three atoms. If you are porting from a Packmol
    /// `.inp` file (which uses 1-based indices), subtract 1 at the
    /// call site.
    pub fn with_atom_restraint(mut self, indices: &[usize], r: impl Restraint + 'static) -> Self {
        self.atom_restraints.push((indices.to_vec(), Arc::new(r)));
        self
    }

    /// Attach an in-loop relaxer for this target.
    ///
    /// Multiple relaxers can be attached (called in order).
    /// Relaxers require `count == 1` because all copies share reference coords.
    ///
    /// Mirrors [`with_restraint`](Self::with_restraint) — a per-target builder method.
    pub fn with_relaxer(mut self, relaxer: impl Relaxer + 'static) -> Self {
        assert!(
            self.count <= 1,
            "relaxers require count == 1 (all copies share ref coords)"
        );
        self.relaxers.push(Box::new(relaxer));
        self
    }

    /// Structure-level budget for the perturbation heuristic
    /// (Packmol's `maxmove`). Defaults to `count` when unset.
    pub fn with_perturb_budget(mut self, n: usize) -> Self {
        self.perturb_budget = Some(n);
        self
    }

    /// Set the centering policy.
    ///
    /// - [`CenteringMode::Auto`] (default): free molecules centered,
    ///   fixed molecules kept in place.
    /// - [`CenteringMode::Center`]: always center.
    /// - [`CenteringMode::Off`]: keep input coordinates unchanged.
    pub fn with_centering(mut self, mode: CenteringMode) -> Self {
        self.centering = mode;
        self
    }

    /// Rotation bound on a single Euler axis, analogous to Packmol's
    /// `constrain_rotation <axis> <center> <delta>`. Arguments are
    /// [`Angle`] values — `Angle::from_degrees(30.0)` or
    /// `Angle::from_radians(FRAC_PI_6)`.
    pub fn with_rotation_bound(mut self, axis: Axis, center: Angle, half_width: Angle) -> Self {
        let idx = match axis {
            // Internal index order follows Packmol's Euler variable order
            // `[beta(y), gama(z), teta(x)]`.
            Axis::Y => 0,
            Axis::Z => 1,
            Axis::X => 2,
        };
        self.rotation_bound[idx] = Some((center, half_width));
        self
    }

    /// Fix this molecule at a specific position with zero rotation.
    ///
    /// Forces `count` to 1 — a fixed molecule is by definition a single
    /// copy. Pair with [`with_orientation`][Self::with_orientation] if
    /// a non-zero Euler orientation is needed.
    pub fn fixed_at(mut self, position: [F; 3]) -> Self {
        assert!(
            self.count <= 1,
            "fixed_at() requires count <= 1, got count = {}. \
             A fixed target is a single placed copy.",
            self.count
        );
        self.fixed_at = Some(Placement {
            position,
            orientation: [Angle::ZERO; 3],
        });
        self.count = 1;
        self
    }

    /// Set the Euler orientation of a previously-fixed target. Must be
    /// called after [`fixed_at`][Self::fixed_at]; panics otherwise.
    pub fn with_orientation(mut self, orientation: [Angle; 3]) -> Self {
        let placement = self.fixed_at.as_mut().expect(
            "with_orientation() requires a prior .fixed_at(pos) call — \
             orientation is only meaningful on fixed targets",
        );
        placement.orientation = orientation;
        self
    }

    pub fn natoms(&self) -> usize {
        self.ref_coords.len()
    }
}

fn centered_coords(coords: &[[F; 3]]) -> Vec<[F; 3]> {
    let (cx, cy, cz) = geometric_center(coords);
    coords
        .iter()
        .map(|p| [p[0] - cx, p[1] - cy, p[2] - cz])
        .collect()
}

fn geometric_center(coords: &[[F; 3]]) -> (F, F, F) {
    if coords.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let n = coords.len() as F;
    let cx = coords.iter().map(|p| p[0]).sum::<F>() / n;
    let cy = coords.iter().map(|p| p[1]).sum::<F>() / n;
    let cz = coords.iter().map(|p| p[2]).sum::<F>() / n;
    (cx, cy, cz)
}
