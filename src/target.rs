//! Target builder for molecular packing.

use std::sync::Arc;

use crate::restraint::{AtomRestraint, Restraint};
use molrs::core::BondDistanceWeights;
use molrs::op::F;

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
    pub template: Option<molrs::core::Frame>,
    /// Intramolecular skip table. This is not `ForceField::special_bonds`.
    ///
    /// Default is [`BondDistanceWeights::from_exclusion_depth`]`(3)`
    /// (`[0, 0, 0, 1]`: 1-2/1-3/1-4 exempt, 1-5+ scored). Set via
    /// [`with_special_bonds`](Self::with_special_bonds).
    pub special_bonds: BondDistanceWeights,
    /// Atoms the lattice walk treats as hydrogens: decorations placed off
    /// their backbone neighbour, never on a lattice site. `None` is the
    /// all-atom default — atoms whose element symbol is `H`. Set via
    /// [`with_hydrogens`](Self::with_hydrogens).
    pub hydrogens: Option<Vec<usize>>,
    /// Per-copy total mass override in amu, for the entries'
    /// `with_density` when element symbols cannot provide it (coarse-grained beads, `from_coords`
    /// targets whose elements are `"X"`).
    pub mass: Option<F>,
}

impl Target {
    /// Create a new target from a `molrs::core::Frame` (read from PDB/XYZ) and a copy count.
    ///
    /// Positions are extracted from the `"atoms"` block (`"x"`, `"y"`, `"z"` columns)
    /// and automatically centered at the geometric center.
    /// VdW radii and element symbols are looked up from the `"element"` column.
    ///
    /// # Panics
    /// Panics if the frame has no `"atoms"` block with `"x"` / `"y"` / `"z"`
    /// float columns.
    pub fn new(frame: molrs::core::Frame, count: usize) -> Self {
        let positions = crate::template::coord_rows(
            &frame
                .coords()
                .expect("target frame needs an 'atoms' block with x / y / z columns"),
        );
        let (radii, elements) = radii_and_elements(&frame, positions.len());
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
            hydrogens: None,
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

    /// Name the atoms growth treats as hydrogens (0-based), replacing the
    /// element-symbol default. A coarse-grained model with no hydrogens passes
    /// `&[]`; a united-atom or unusual naming scheme passes its own list.
    ///
    /// # Panics
    /// If an index is out of range.
    pub fn with_hydrogens(mut self, indices: &[usize]) -> Self {
        let n = self.natoms();
        if let Some(&bad) = indices.iter().find(|&&i| i >= n) {
            panic!("with_hydrogens: atom index {bad} out of range for {n} atoms");
        }
        self.hydrogens = Some(indices.to_vec());
        self
    }

    /// Per-atom hydrogen flag: [`hydrogens`](Self::hydrogens) when set, else
    /// element symbol `H`.
    pub(crate) fn hydrogen_mask(&self) -> Vec<bool> {
        match &self.hydrogens {
            Some(indices) => {
                let mut mask = vec![false; self.natoms()];
                for &i in indices {
                    mask[i] = true;
                }
                mask
            }
            None => self
                .elements
                .iter()
                .map(|e| e.eq_ignore_ascii_case("H"))
                .collect(),
        }
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
    /// ([`crate::Restraint`]) — it sees every copy's coordinate at
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
    /// `GenCanPack::new().run(&[Target::fixed_from(&grown.frame), solvent], …)`.
    pub fn fixed_from(frame: &molrs::core::Frame) -> Self {
        Self::new(frame.clone(), 1)
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

/// `coords` shifted so their geometric center sits at the origin.
pub(crate) fn centered_coords(coords: &[[F; 3]]) -> Vec<[F; 3]> {
    let (cx, cy, cz) = geometric_center(coords);
    coords
        .iter()
        .map(|p| [p[0] - cx, p[1] - cy, p[2] - cz])
        .collect()
}

/// The template's geometric center, summed axis by axis in atom order.
///
/// Kept local rather than `molrs::op::centroid` on purpose
/// (module-responsibility ruling 10): this is Packmol's template centering,
/// and its exact summation order and seed value fix the bits of every
/// centered template — and so of every packed coordinate the Packmol-parity
/// goldens pin. Changing it is a change to Packmol parity, not a refactor.
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

/// VdW radius and trimmed symbol per atom, from the `atoms` block's
/// `element` column. An unknown symbol — or no column at all — is `"X"` at
/// 1.5 Å.
fn radii_and_elements(frame: &molrs::core::Frame, n: usize) -> (Vec<F>, Vec<String>) {
    let column = frame
        .get("atoms")
        .and_then(|atoms| atoms.get(molrs::core::keys::ELEMENT))
        .and_then(molrs::core::Column::as_string);
    let Some(symbols) = column else {
        return (vec![1.5; n], vec!["X".to_string(); n]);
    };
    symbols
        .iter()
        .map(|sym| {
            let sym = sym.trim();
            let radius = molrs::core::Element::by_symbol(sym)
                .map(|e| e.vdw_radius() as F)
                .unwrap_or(1.5);
            (radius, sym.to_string())
        })
        .unzip()
}

#[cfg(test)]
mod tests {
    //! Tests for Target builder: construction, natoms/count, fixed_at,
    //! centering modes, restraint attachment, and hook validation.

    use crate::{GenCanPack, PackEngine, RegionRestraint, Target};
    use molrs::core::BondDistanceWeights;
    use molrs::op::F;
    use std::sync::Arc;

    use crate::testutil::inside_box;
    use molrs::core::Sphere;
    use ndarray::array;

    // ── molrs regions lifted to "stay inside" (the one geometric restraint) ─────

    fn inside_sphere(center: [F; 3], radius: F) -> RegionRestraint {
        RegionRestraint(Arc::new(Sphere::new(
            array![center[0], center[1], center[2]],
            radius,
        )))
    }

    // ── helpers ────────────────────────────────────────────────────────────────

    fn water_positions() -> Vec<[F; 3]> {
        vec![
            [0.0, 0.0, 0.0],    // O
            [0.96, 0.0, 0.0],   // H
            [-0.24, 0.93, 0.0], // H
        ]
    }

    fn water_radii() -> Vec<F> {
        vec![1.52, 1.20, 1.20]
    }

    // ── construction ───────────────────────────────────────────────────────────

    #[test]
    fn from_coords_basic() {
        let t = Target::from_coords(&water_positions(), &water_radii(), 5);
        assert_eq!(t.natoms(), 3);
        assert_eq!(t.count, 5);
        assert!(t.name.is_none());
        assert!(t.fixed_at.is_none());
    }

    #[test]
    fn with_name() {
        let t = Target::from_coords(&water_positions(), &water_radii(), 1).with_name("water");
        assert_eq!(t.name, Some("water".to_string()));
    }

    #[test]
    fn ref_coords_are_centered() {
        let coords = vec![[10.0, 20.0, 30.0], [12.0, 20.0, 30.0]];
        let t = Target::from_coords(&coords, &[1.0, 1.0], 1);
        // Center should be at (11, 20, 30), so ref_coords are [-1, 0, 0] and [1, 0, 0]
        let cx: F = t.ref_coords.iter().map(|p| p[0]).sum::<F>() / t.natoms() as F;
        let cy: F = t.ref_coords.iter().map(|p| p[1]).sum::<F>() / t.natoms() as F;
        let cz: F = t.ref_coords.iter().map(|p| p[2]).sum::<F>() / t.natoms() as F;
        assert!(cx.abs() < 1e-6, "ref_coords x center should be 0");
        assert!(cy.abs() < 1e-6, "ref_coords y center should be 0");
        assert!(cz.abs() < 1e-6, "ref_coords z center should be 0");
    }

    #[test]
    fn input_coords_preserved() {
        let coords = vec![[10.0, 20.0, 30.0]];
        let t = Target::from_coords(&coords, &[1.0], 1);
        assert!((t.input_coords[0][0] - 10.0).abs() < 1e-6);
        assert!((t.input_coords[0][1] - 20.0).abs() < 1e-6);
    }

    #[test]
    fn new_uses_geometric_center_even_when_elements_are_known() {
        use molrs::core::Block;
        use ndarray::Array1;

        let mut atoms = Block::new();
        atoms
            .insert("x", Array1::from_vec(vec![0.0, 1.0]).into_dyn())
            .expect("insert x");
        atoms
            .insert("y", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
            .expect("insert y");
        atoms
            .insert("z", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
            .expect("insert z");
        atoms
            .insert(
                "element",
                Array1::from_vec(vec!["O".to_string(), "H".to_string()]).into_dyn(),
            )
            .expect("insert element");

        let mut frame = molrs::core::Frame::new();
        frame.insert("atoms", atoms);

        let t = Target::new(frame, 1);
        let arithmetic_center = (t.ref_coords[0][0] + t.ref_coords[1][0]) / 2.0;
        assert!(
            arithmetic_center.abs() < 1e-6,
            "geometry center should be zero"
        );
    }

    // ── restraints ─────────────────────────────────────────────────────────────

    #[test]
    fn with_restraint() {
        let t = Target::from_coords(&water_positions(), &water_radii(), 5)
            .with_restraint(inside_box([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]));
        assert_eq!(t.molecule_restraints.len(), 1);
    }

    #[test]
    fn with_restraint_chained() {
        let t = Target::from_coords(&water_positions(), &water_radii(), 5)
            .with_restraint(inside_box([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]))
            .with_restraint(inside_sphere([10.0, 10.0, 10.0], 50.0));
        assert_eq!(t.molecule_restraints.len(), 2);
    }

    #[test]
    fn with_atom_restraint() {
        // Indices are now 0-based (matching Rust convention) — no internal
        // conversion happens. Caller subtracts 1 when porting from Packmol
        // `.inp` files.
        let t = Target::from_coords(&water_positions(), &water_radii(), 5)
            .with_atom_restraint(&[0, 1], inside_sphere([0.0, 0.0, 0.0], 5.0));
        assert_eq!(t.atom_restraints.len(), 1);
        assert_eq!(t.atom_restraints[0].0, vec![0, 1]);
    }

    // ── fixed placement ────────────────────────────────────────────────────────

    #[test]
    fn fixed_at_sets_count_to_1() {
        // Explicit count=1 to satisfy the assertion in fixed_at().
        let t =
            Target::from_coords(&water_positions(), &water_radii(), 1).fixed_at([0.0, 0.0, 0.0]);
        assert_eq!(t.count, 1);
        assert!(t.fixed_at.is_some());
        let fp = t.fixed_at.unwrap();
        assert!((fp.position[0]).abs() < 1e-6);
        assert!((fp.orientation[0].radians()).abs() < 1e-6);
    }

    #[test]
    fn fixed_at_with_orientation() {
        use crate::Angle;
        let t = Target::from_coords(&water_positions(), &water_radii(), 1)
            .fixed_at([1.0, 2.0, 3.0])
            .with_orientation([
                Angle::from_radians(0.1),
                Angle::from_radians(0.2),
                Angle::from_radians(0.3),
            ]);
        assert_eq!(t.count, 1);
        let fp = t.fixed_at.unwrap();
        assert!((fp.position[0] - 1.0).abs() < 1e-6);
        assert!((fp.orientation[2].radians() - 0.3).abs() < 1e-6);
    }

    #[test]
    fn fixed_target_auto_centering_disabled() {
        // When fixed_at is used with Auto centering (default), the fixed molecule
        // should NOT be centered — its input coords are used directly.
        let free = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1)
            .with_restraint(inside_box([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0]));
        let fixed = Target::from_coords(&[[10.0, 0.0, 0.0], [12.0, 0.0, 0.0]], &[1.0, 1.0], 1)
            .fixed_at([0.0, 0.0, 0.0]);

        let result = GenCanPack::new()
            .with_seed(1)
            .run(&[free, fixed], 5)
            .expect("pack should succeed");

        // Fixed atoms follow free atoms in output.
        assert!((result.positions()[1][0] - 10.0).abs() < 1e-6);
        assert!((result.positions()[2][0] - 12.0).abs() < 1e-6);
    }

    #[test]
    fn fixed_target_centered() {
        let free = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1)
            .with_restraint(inside_box([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0]));
        let fixed = Target::from_coords(&[[10.0, 0.0, 0.0], [12.0, 0.0, 0.0]], &[1.0, 1.0], 1)
            .with_centering(crate::CenteringMode::Center)
            .fixed_at([0.0, 0.0, 0.0]);

        let result = GenCanPack::new()
            .with_seed(1)
            .run(&[free, fixed], 5)
            .expect("pack should succeed");

        // COM of [10,12] = 11. After centering, ref_coords = [-1, +1].
        // Placed at origin → positions = [-1, +1].
        assert!((result.positions()[1][0] + 1.0).abs() < 1e-6);
        assert!((result.positions()[2][0] - 1.0).abs() < 1e-6);
    }

    // ── centering modes ────────────────────────────────────────────────────────

    #[test]
    fn with_centering_center() {
        use crate::CenteringMode;
        let t = Target::from_coords(&water_positions(), &water_radii(), 1)
            .with_centering(CenteringMode::Center);
        assert_eq!(t.centering, CenteringMode::Center);
    }

    #[test]
    fn with_centering_off() {
        use crate::CenteringMode;
        let t = Target::from_coords(&water_positions(), &water_radii(), 1)
            .with_centering(CenteringMode::Off);
        assert_eq!(t.centering, CenteringMode::Off);
    }

    // ── rotation constraints ───────────────────────────────────────────────────

    #[test]
    fn with_rotation_bound() {
        use crate::{Angle, Axis};
        let t = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1)
            .with_rotation_bound(Axis::X, Angle::from_degrees(0.0), Angle::from_degrees(10.0))
            .with_rotation_bound(Axis::Y, Angle::from_degrees(90.0), Angle::from_degrees(5.0))
            .with_rotation_bound(
                Axis::Z,
                Angle::from_degrees(180.0),
                Angle::from_degrees(15.0),
            );
        // Euler variable order: [beta(Y) = 0, gama(Z) = 1, teta(X) = 2]
        assert!(t.rotation_bound[0].is_some()); // beta (Y)
        assert!(t.rotation_bound[1].is_some()); // gama (Z)
        assert!(t.rotation_bound[2].is_some()); // teta (X)
    }

    // ── perturb budget ─────────────────────────────────────────────────────────

    #[test]
    fn with_perturb_budget() {
        let t = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1).with_perturb_budget(5);
        assert_eq!(t.perturb_budget, Some(5));
    }

    // ── default element ────────────────────────────────────────────────────────

    #[test]
    fn default_elements_are_x() {
        let t = Target::from_coords(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &[1.0, 1.0], 1);
        assert_eq!(t.elements, vec!["X", "X"]);
    }

    // ── special bonds ──────────────────────────────────────────────────────────

    /// Hydrogens default to element symbol `H` (case-insensitive); a stated
    /// list replaces the rule outright.
    #[test]
    fn hydrogens_default_to_element_and_can_be_stated() {
        let mut t =
            Target::from_coords(&[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], &[1.0; 3], 1);
        t.elements = vec!["C".into(), "h".into(), "X".into()];
        assert_eq!(t.hydrogen_mask(), [false, true, false]);
        assert_eq!(t.clone().with_hydrogens(&[]).hydrogen_mask(), [false; 3]);
        assert_eq!(t.with_hydrogens(&[2]).hydrogen_mask(), [false, false, true]);
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn hydrogens_out_of_range_panics() {
        let _ = Target::from_coords(&[[0.0; 3]], &[1.0], 1).with_hydrogens(&[1]);
    }

    /// Default skip table is depth 3: Cassandra `[0, 0, 0, 1]` (1-2/1-3/1-4
    /// exempt, 1-5+ scored). Written once in `from_parts`, so both constructors
    /// carry it.
    #[test]
    fn target_default_special_bonds_is_depth_3() {
        use molrs::core::Block;
        use ndarray::Array1;

        let from_coords = Target::from_coords(&water_positions(), &water_radii(), 5);
        assert_eq!(from_coords.special_bonds.as_slice(), &[0.0, 0.0, 0.0, 1.0]);

        let mut atoms = Block::new();
        atoms
            .insert("x", Array1::from_vec(vec![0.0, 1.0]).into_dyn())
            .expect("insert x");
        atoms
            .insert("y", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
            .expect("insert y");
        atoms
            .insert("z", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
            .expect("insert z");
        let mut frame = molrs::core::Frame::new();
        frame.insert("atoms", atoms);

        let from_frame = Target::new(frame, 1);
        assert_eq!(from_frame.special_bonds.as_slice(), &[0.0, 0.0, 0.0, 1.0]);
    }

    /// Depth 2 is `[0, 0, 1]`. The builder returns `Self` (chainable), not
    /// `Result` — fractional tables are stored here and refused later at
    /// `InternalTree::from_frame`.
    #[test]
    fn target_with_special_bonds_stores_depth_2() {
        let t = Target::from_coords(&water_positions(), &water_radii(), 1)
            .with_special_bonds(BondDistanceWeights::from_exclusion_depth(2))
            .with_name("kg");
        assert_eq!(t.special_bonds.as_slice(), &[0.0, 0.0, 1.0]);
        assert_eq!(t.name.as_deref(), Some("kg"));
    }

    #[test]
    fn target_with_special_bonds_stores_fractional() {
        let table = BondDistanceWeights::new(vec![0.0, 0.0, 0.5, 1.0])
            .expect("fractional 1-4 is a legal BondDistanceWeights table");
        let t =
            Target::from_coords(&water_positions(), &water_radii(), 1).with_special_bonds(table);
        assert_eq!(t.special_bonds.as_slice(), &[0.0, 0.0, 0.5, 1.0]);
    }

    // ── panics ─────────────────────────────────────────────────────────────────

    #[test]
    #[should_panic(expected = "positions and radii must have the same length")]
    fn mismatched_coords_and_radii_panics() {
        Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0, 2.0], 1);
    }
}

#[cfg(test)]
mod atom_property_tests {
    //! Per-atom packing properties: `radius`, `fscale`, `short_radius`,
    //! `short_radius_scale`.
    //!
    //! All four follow the same Packmol scheme (`app/packmol.f90` lines 281-515):
    //! a global default, a structure-level keyword covering every atom of every
    //! copy, and an atom-specific keyword inside an `atoms ... end atoms` block
    //! that overrides the selected atoms. Per-atom values are a **per-type
    //! template** — every copy of a type gets the same ones.

    use crate::{GenCanPack, PackEngine, Target};
    use molrs::op::F;

    // ── molrs regions lifted to "stay inside" (the one geometric restraint) ─────

    fn coords(n: usize) -> Vec<[F; 3]> {
        (0..n).map(|i| [i as F * 3.0, 0.0, 0.0]).collect()
    }

    fn target(n: usize, count: usize) -> Target {
        Target::from_coords(&coords(n), &vec![1.5; n], count)
    }

    // ── layer 1: the global default ────────────────────────────────────────────

    #[test]
    fn radius_defaults_to_the_packer_default() {
        let t = target(3, 2);
        assert_eq!(t.resolved_radii(2.0), vec![2.0; 3]);
    }

    /// The frame's van der Waals radii are *not* the packing radii — Packmol packs
    /// on `tolerance / 2` unless told otherwise.
    #[test]
    fn frame_vdw_radii_do_not_become_packing_radii() {
        let t = Target::from_coords(&coords(3), &[9.0, 9.0, 9.0], 1);
        assert_eq!(t.resolved_radii(2.0), vec![2.0; 3]);
    }

    // ── layer 2: structure-level ───────────────────────────────────────────────

    #[test]
    fn structure_radius_covers_every_atom() {
        let t = target(3, 2).with_radius(3.5);
        assert_eq!(t.resolved_radii(2.0), vec![3.5; 3]);
    }

    // ── layer 3: atom-specific ─────────────────────────────────────────────────

    #[test]
    fn atom_radius_applies_only_to_the_selected_atoms() {
        let t = target(3, 1).with_atom_radius(&[1], 5.0);
        assert_eq!(t.resolved_radii(2.0), vec![2.0, 5.0, 2.0]);
    }

    /// Packmol applies the structure-level pass first and lets the atom-specific
    /// pass overwrite it.
    #[test]
    fn atom_radius_overrides_the_structure_radius() {
        let t = target(3, 1).with_radius(3.0).with_atom_radius(&[2], 6.0);
        assert_eq!(t.resolved_radii(2.0), vec![3.0, 3.0, 6.0]);
    }

    /// A structure-level radius set *after* an atom-specific one still wins for
    /// every atom — the builder applies calls in order, so the caller controls
    /// precedence explicitly rather than by a hidden rule.
    #[test]
    fn a_later_structure_radius_replaces_earlier_atom_radii() {
        let t = target(3, 1).with_atom_radius(&[0], 9.0).with_radius(3.0);
        assert_eq!(t.resolved_radii(2.0), vec![3.0; 3]);
    }

    #[test]
    fn repeated_atom_radius_calls_accumulate() {
        let t = target(4, 1)
            .with_atom_radius(&[0, 1], 5.0)
            .with_atom_radius(&[3], 7.0);
        assert_eq!(t.resolved_radii(2.0), vec![5.0, 5.0, 2.0, 7.0]);
    }

    #[test]
    #[should_panic(expected = "atom index")]
    fn atom_radius_rejects_an_out_of_range_index() {
        let _ = target(3, 1).with_atom_radius(&[3], 5.0);
    }

    #[test]
    #[should_panic(expected = "positive")]
    fn radius_rejects_a_non_positive_value() {
        let _ = target(3, 1).with_radius(0.0);
    }

    // ── the packer honours them ────────────────────────────────────────────────

    // ── Packmol `.inp` parity ──────────────────────────────────────────────────

    // ── fscale ─────────────────────────────────────────────────────────────────
    //
    // Packmol's per-atom weight on the overlap penalty: the pair term is scaled by
    // `fscale_i * fscale_j` (`fparc.f90`), so a species can be made "softer" than
    // the rest without changing how far apart it is asked to sit.

    #[test]
    fn fscale_defaults_to_one() {
        assert_eq!(target(3, 1).resolved_fscale(), vec![1.0; 3]);
    }

    #[test]
    fn structure_fscale_covers_every_atom() {
        let t = target(3, 1).with_fscale(0.25);
        assert_eq!(t.resolved_fscale(), vec![0.25; 3]);
    }

    #[test]
    fn atom_fscale_overrides_the_structure_value() {
        let t = target(3, 1).with_fscale(0.5).with_atom_fscale(&[1], 2.0);
        assert_eq!(t.resolved_fscale(), vec![0.5, 2.0, 0.5]);
    }

    #[test]
    #[should_panic(expected = "positive")]
    fn fscale_rejects_a_non_positive_value() {
        let _ = target(3, 1).with_fscale(0.0);
    }

    // ── short radius ───────────────────────────────────────────────────────────
    //
    // Packmol's optional second, shorter-range penalty (`use_short_tol` /
    // `short_tol_dist` / `short_tol_scale`, overridable per structure and per
    // atom). It lets a pair approach past the main radius while still being
    // stopped hard at a smaller one.

    #[test]
    fn short_radius_is_off_until_asked_for() {
        let t = target(3, 1);
        assert!(!t.uses_short_radius().iter().any(|&b| b));
    }

    #[test]
    fn structure_short_radius_enables_it_for_every_atom() {
        let t = target(3, 1).with_short_radius(0.5);
        assert_eq!(t.resolved_short_radii(1.0), vec![0.5; 3]);
        assert_eq!(t.uses_short_radius(), vec![true; 3]);
    }

    #[test]
    fn atom_short_radius_enables_only_the_selected_atoms() {
        let t = target(3, 1).with_atom_short_radius(&[2], 0.5);
        assert_eq!(t.uses_short_radius(), vec![false, false, true]);
    }

    /// Setting only the scale still opts the atom in, as Packmol does.
    #[test]
    fn short_radius_scale_alone_enables_the_penalty() {
        let t = target(2, 1).with_atom_short_radius_scale(&[0], 5.0);
        assert_eq!(t.uses_short_radius(), vec![true, false]);
        assert_eq!(t.resolved_short_radius_scale(3.0), vec![5.0, 3.0]);
    }

    #[test]
    #[should_panic(expected = "smaller than the tolerance")]
    fn global_short_tolerance_must_be_below_the_tolerance() {
        let _ = GenCanPack::new()
            .with_tolerance(2.0)
            .with_short_tolerance(4.0, 3.0);
    }

    // ── the other three keywords, at both script levels ────────────────────────
}
