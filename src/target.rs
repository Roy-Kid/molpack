//! Target builder for molecular packing.

use std::sync::Arc;

use crate::frame::frame_to_coords_and_elements;
use crate::restraint::{AtomRestraint, Restraint};
use molrs::BondDistanceWeights;
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

/// Describes one type of molecule to be packed.
#[derive(Debug, Clone)]
pub struct Target {
    /// Input coordinates as provided by the source structure.
    pub input_coords: Vec<[F; 3]>,
    /// Flat list of atom positions — the centered reference coordinates.
    /// Shape: natoms × 3, stored as Vec<[F; 3]>.
    pub ref_coords: Vec<[F; 3]>,
    /// Van der Waals radii per atom, as read from the source structure.
    ///
    /// These are *reference* radii — they are **not** what the packer separates
    /// atoms by. Packing radii come from [`resolved_radii`](Self::resolved_radii);
    /// see [`with_radius`](Self::with_radius) for why.
    pub radii: Vec<F>,
    /// Per-atom packing-radius overrides, one entry per atom. `None` means
    /// "use the packer's default" (`tolerance / 2`).
    ///
    /// This is a **per-type template**: every copy of this target is packed
    /// with the same values, matching Packmol, which broadcasts the first
    /// copy's per-atom values to the rest (`app/packmol.f90` lines 503-515).
    /// The same holds for the three fields below.
    pub atom_radii: Vec<Option<F>>,
    /// Per-atom overlap-penalty weights. `None` means the Packmol default of
    /// `1.0`. The pair term is scaled by `fscale_i * fscale_j`.
    pub atom_fscale: Vec<Option<F>>,
    /// Per-atom short-radius overrides for Packmol's optional second,
    /// shorter-range penalty. `None` means the packer's global short radius.
    pub atom_short_radii: Vec<Option<F>>,
    /// Per-atom weights for that second penalty. `None` means the packer's
    /// global short-radius scale.
    pub atom_short_radius_scale: Vec<Option<F>>,
    /// Element symbols per atom (e.g. `"C"`, `"O"`). Defaults to `"X"` if unknown.
    pub elements: Vec<String>,
    /// Number of copies to pack.
    pub count: usize,
    /// Optional name for logging.
    pub name: Option<String>,
    /// Restraints applied to every atom of every molecule copy.
    pub molecule_restraints: Vec<Arc<dyn AtomRestraint>>,
    /// Per-atom-subset restraints: `(atom_indices_0_based, restraint)`.
    /// Each entry holds the 0-based atom indices (converted from Packmol's
    /// 1-based convention at registration time) and the restraint applied to them.
    pub atom_restraints: Vec<(Vec<usize>, Arc<dyn AtomRestraint>)>,
    /// Group-level restraints evaluated over **all copies** of this type at
    /// once (e.g. distribution matching). Unlike `molecule_restraints`, these
    /// couple the copies through their joint coordinate, so they cannot be
    /// expressed as a per-atom [`AtomRestraint`]. See [`Restraint`].
    pub collective_restraints: Vec<Arc<dyn Restraint>>,
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
    /// Source frame this target was built from, retained so the packer can
    /// replay its full topology (bonds/angles/…) and per-atom metadata onto
    /// the packed coordinates. `None` for targets built from bare coordinates
    /// ([`Target::from_coords`]), whose result frame is coordinates-only.
    pub template: Option<molrs::Frame>,
    /// Intramolecular skip table. This is not `ForceField::special_bonds`.
    ///
    /// Default is [`BondDistanceWeights::from_exclusion_depth`]`(3)`
    /// (`[0, 0, 0, 1]`: 1-2/1-3/1-4 exempt, 1-5+ scored). Set via
    /// [`with_special_bonds`](Self::with_special_bonds).
    pub special_bonds: BondDistanceWeights,
    /// Per-copy total mass override in amu, for the entries'
    /// `with_density` when element symbols cannot provide it (coarse-grained beads, `from_coords`
    /// targets whose elements are `"X"`).
    pub mass: Option<F>,
}

impl Target {
    /// Create a new target from a `molrs::Frame` (read from PDB/XYZ) and a copy count.
    ///
    /// Positions are extracted from the `"atoms"` block (`"x"`, `"y"`, `"z"` columns)
    /// and automatically centered at the geometric center.
    /// VdW radii and element symbols are looked up from the `"element"` column.
    pub fn new(frame: molrs::Frame, count: usize) -> Self {
        let (positions, radii, elements) = frame_to_coords_and_elements(&frame);
        let mut target = Self::from_parts(&positions, &radii, elements, count);
        target.template = Some(frame);
        target
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
            atom_radii: vec![None; radii.len()],
            atom_fscale: vec![None; radii.len()],
            atom_short_radii: vec![None; radii.len()],
            atom_short_radius_scale: vec![None; radii.len()],
            elements,
            count,
            name: None,
            molecule_restraints: Vec::new(),
            atom_restraints: Vec::new(),
            collective_restraints: Vec::new(),
            perturb_budget: None,
            centering: CenteringMode::Auto,
            rotation_bound: [None, None, None],
            fixed_at: None,
            template: None,
            special_bonds: BondDistanceWeights::from_exclusion_depth(3),
            mass: None,
        }
    }

    /// Override the per-copy total mass (amu) used by
    /// [`PackEngine::with_density`](crate::PackEngine::with_density). Needed when
    /// element symbols cannot resolve a mass (CG beads, bare-coordinate
    /// targets).
    pub fn with_mass(mut self, amu: F) -> Self {
        self.mass = Some(amu);
        self
    }

    /// This is not `ForceField::special_bonds`.
    ///
    /// Stores `table` as this target's intramolecular skip weights. Default is
    /// [`BondDistanceWeights::from_exclusion_depth`]`(3)` (`[0, 0, 0, 1]`).
    /// All-atom templates with explicit hydrogen keep that depth-3 table and
    /// shrink hydrogen via [`with_atom_radius`](Self::with_atom_radius); do
    /// not replace it with `[0, 0, 0, 0, 0, 1]`.
    ///
    /// Fractional weights are stored here and refused later when growth
    /// compiles the skip set ([`crate::grow::GrowError::NonBinarySpecialBond`]).
    pub fn with_special_bonds(mut self, table: BondDistanceWeights) -> Self {
        self.special_bonds = table;
        self
    }

    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Set the packing radius for **every atom** of this target.
    ///
    /// Packmol's structure-level `radius` keyword. The packer separates two
    /// atoms by the sum of their radii, so this is how a species is made
    /// bulkier or slimmer than the global `tolerance / 2` default.
    ///
    /// Van der Waals radii read from the source structure are deliberately not
    /// used for this: Packmol packs on a single tolerance so the minimum
    /// separation is a property of the run, not of whichever force field
    /// labelled the input. Opting in per species is what this method is for.
    ///
    /// Order matters — a later call overwrites earlier per-atom values, which
    /// is how the two Packmol passes compose:
    /// ```
    /// # use molpack::Target;
    /// let t = Target::from_coords(&[[0.0; 3], [1.0, 0.0, 0.0]], &[1.5; 2], 1)
    ///     .with_radius(3.0)             // every atom
    ///     .with_atom_radius(&[1], 6.0); // then one of them
    /// assert_eq!(t.resolved_radii(2.0), vec![3.0, 6.0]);
    /// ```
    ///
    /// # Panics
    /// If `radius` is not positive.
    pub fn with_radius(mut self, radius: F) -> Self {
        set_all(&mut self.atom_radii, radius, "packing radius");
        self
    }

    /// Set the packing radius for selected atoms of this target.
    ///
    /// Packmol's `radius` inside an `atoms ... end atoms` block. Indices are
    /// **0-based**, matching [`with_atom_restraint`](Self::with_atom_restraint);
    /// a Packmol `.inp` uses 1-based indices, so subtract one when porting.
    ///
    /// Applies to the same atom of every copy — per-atom values are a per-type
    /// template, not a per-copy one.
    ///
    /// # Panics
    /// If `radius` is not positive, or an index is out of range.
    pub fn with_atom_radius(mut self, indices: &[usize], radius: F) -> Self {
        set_at(&mut self.atom_radii, indices, radius, "packing radius");
        self
    }

    /// Weight this target's atoms in the overlap penalty.
    ///
    /// Packmol's structure-level `fscale`. The pair term is multiplied by
    /// `fscale_i * fscale_j`, so a value below 1 makes a species *softer* —
    /// penalised less for the same overlap — without changing the distance it
    /// is asked to keep. Default `1.0`.
    ///
    /// # Panics
    /// If `fscale` is not positive.
    pub fn with_fscale(mut self, fscale: F) -> Self {
        set_all(&mut self.atom_fscale, fscale, "fscale");
        self
    }

    /// Weight selected atoms in the overlap penalty. Packmol's `fscale` inside
    /// an `atoms ... end atoms` block; indices are **0-based**.
    ///
    /// # Panics
    /// If `fscale` is not positive, or an index is out of range.
    pub fn with_atom_fscale(mut self, indices: &[usize], fscale: F) -> Self {
        set_at(&mut self.atom_fscale, indices, fscale, "fscale");
        self
    }

    /// Give this target's atoms a second, shorter penalty radius.
    ///
    /// Packmol's structure-level `short_radius`. The main radius still governs
    /// the ordinary overlap term; inside this smaller radius an additional,
    /// steeper penalty applies, which lets a pair approach past the main
    /// radius while still being stopped hard. Setting it opts the atoms into
    /// the short-radius term.
    ///
    /// Must be **smaller** than the atom's packing radius — `pack()` rejects
    /// the run otherwise, as Packmol does.
    ///
    /// # Panics
    /// If `short_radius` is not positive.
    pub fn with_short_radius(mut self, short_radius: F) -> Self {
        set_all(&mut self.atom_short_radii, short_radius, "short radius");
        self
    }

    /// Per-atom counterpart of [`with_short_radius`](Self::with_short_radius);
    /// indices are **0-based**.
    ///
    /// # Panics
    /// If `short_radius` is not positive, or an index is out of range.
    pub fn with_atom_short_radius(mut self, indices: &[usize], short_radius: F) -> Self {
        set_at(
            &mut self.atom_short_radii,
            indices,
            short_radius,
            "short radius",
        );
        self
    }

    /// Weight the short-radius penalty for this target's atoms.
    ///
    /// Packmol's structure-level `short_radius_scale`. Like
    /// [`with_short_radius`](Self::with_short_radius), setting it opts the
    /// atoms into the short-radius term even on its own.
    ///
    /// # Panics
    /// If `scale` is not positive.
    pub fn with_short_radius_scale(mut self, scale: F) -> Self {
        set_all(
            &mut self.atom_short_radius_scale,
            scale,
            "short radius scale",
        );
        self
    }

    /// Per-atom counterpart of
    /// [`with_short_radius_scale`](Self::with_short_radius_scale); indices are
    /// **0-based**.
    ///
    /// # Panics
    /// If `scale` is not positive, or an index is out of range.
    pub fn with_atom_short_radius_scale(mut self, indices: &[usize], scale: F) -> Self {
        set_at(
            &mut self.atom_short_radius_scale,
            indices,
            scale,
            "short radius scale",
        );
        self
    }

    /// The packing radius of each atom, resolving unset entries to `default`
    /// (the packer passes `tolerance / 2`).
    ///
    /// One entry per atom — the template every copy of this target is packed
    /// with.
    pub fn resolved_radii(&self, default: F) -> Vec<F> {
        resolve(&self.atom_radii, default)
    }

    /// The overlap-penalty weight of each atom; unset entries resolve to `1.0`.
    pub fn resolved_fscale(&self) -> Vec<F> {
        resolve(&self.atom_fscale, 1.0)
    }

    /// The short-radius of each atom, resolving unset entries to `default`
    /// (the packer's global short radius).
    pub fn resolved_short_radii(&self, default: F) -> Vec<F> {
        resolve(&self.atom_short_radii, default)
    }

    /// The short-radius penalty weight of each atom, resolving unset entries
    /// to `default` (the packer's global short-radius scale).
    pub fn resolved_short_radius_scale(&self, default: F) -> Vec<F> {
        resolve(&self.atom_short_radius_scale, default)
    }

    /// Which atoms this target opts into the short-radius penalty.
    ///
    /// True wherever either short-radius field was set — Packmol raises the
    /// same flag from both keywords (`app/packmol.f90` lines 473 and 495).
    /// The packer ORs this with its global short-tolerance switch.
    pub fn uses_short_radius(&self) -> Vec<bool> {
        self.atom_short_radii
            .iter()
            .zip(&self.atom_short_radius_scale)
            .map(|(r, s)| r.is_some() || s.is_some())
            .collect()
    }

    /// Attach a restraint applied to every atom of every molecule copy.
    pub fn with_restraint(mut self, r: impl AtomRestraint + 'static) -> Self {
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
    pub fn with_atom_restraint(
        mut self,
        indices: &[usize],
        r: impl AtomRestraint + 'static,
    ) -> Self {
        self.atom_restraints.push((indices.to_vec(), Arc::new(r)));
        self
    }

    /// Attach a group-level restraint evaluated over all copies of this type at
    /// once (e.g. distribution matching). The restraint sees every copy's
    /// coordinate jointly and returns a coupled gradient.
    ///
    /// Here `Restraint` is the **group/collective** trait
    /// ([`crate::restraint::Restraint`]) — it sees every copy's coordinate at
    /// once, not the per-atom [`AtomRestraint`].
    pub fn with_collective_restraint(mut self, r: impl Restraint + 'static) -> Self {
        self.collective_restraints.push(Arc::new(r));
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

    /// One fixed obstacle target holding a previous pack's entire output,
    /// coordinates kept verbatim (`CenteringMode::Off` + identity placement).
    ///
    /// The named chaining primitive of engine-entry-split: grow first, then
    /// pack the next stage around the grown matrix held fixed —
    /// `GenCanPack::new().run(&[Target::fixed_from(&grown), solvent], …)`.
    pub fn fixed_from(result: &crate::entry::PackResult) -> Self {
        Self::new(result.frame.clone(), 1)
            .with_centering(CenteringMode::Off)
            .fixed_at([0.0; 3])
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

// ── per-atom override helpers ───────────────────────────────────────────────
//
// The four per-atom properties (radius, fscale, short radius, short-radius
// scale) share one shape: a structure-level set that covers every atom and an
// atom-level set that overrides a selection, with `None` deferring to a
// default the packer supplies. These keep that logic in one place.

fn check_positive(value: F, what: &str) {
    assert!(
        value > 0.0 && !value.is_nan(),
        "{what} must be positive, got {value}"
    );
}

fn set_all(slot: &mut [Option<F>], value: F, what: &str) {
    check_positive(value, what);
    slot.fill(Some(value));
}

fn set_at(slot: &mut [Option<F>], indices: &[usize], value: F, what: &str) {
    check_positive(value, what);
    let n = slot.len();
    for &i in indices {
        assert!(
            i < n,
            "atom index {i} is out of range for a {n}-atom structure"
        );
        slot[i] = Some(value);
    }
}

fn resolve(slot: &[Option<F>], default: F) -> Vec<F> {
    slot.iter().map(|v| v.unwrap_or(default)).collect()
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
