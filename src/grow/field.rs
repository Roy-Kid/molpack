//! The overlap field a growing chain is scored against.
//!
//! Holds every atom placed so far in a uniform cell list over the periodic
//! box, and answers one question: *may this atom go here, and how roomy is it?*
//!
//! Two radii per pair. Inside the **hard core**
//! (`(radius_i + radius_j) * hard_scale`, radii in Å, `hard_scale`
//! dimensionless with `1.0` = the same contact distance the GENCAN entry
//! enforces) a placement is refused outright, so a completed structure
//! satisfies molpack's own overlap criterion by construction rather than by
//! convergence. Between the scaled hard core and the **soft shell** (an extra
//! width in Å, charged only for a different-molecule neighbour) the
//! placement is allowed but charged, which is what biases growth towards
//! the roomy directions instead of merely the legal ones.
//!
//! Insertion and removal are both O(1) amortised — the growth driver retracts
//! and regrows constantly, so removal cannot be a rebuild.

use molrs::op::types::F;
use molrs::spatial::neighbors::CellGrid;
use molrs::spatial::{Mic, SimBox};

/// A placed-atom field over an orthorhombic, optionally periodic box.
///
/// Atoms live in a flat *slot* space so one field can carry several species at
/// once: slot `s` belongs to molecule `mol_of[s]` and is that molecule's
/// template atom `atom_of[s]`.
#[derive(Debug, Clone)]
pub struct OverlapField {
    /// The box this field partitions — the single home of origin, lengths and
    /// periodicity.
    bx: SimBox,
    /// Cell indexing and the 3x3x3 stencil, both molrs's: wrap on a periodic
    /// axis, clamp on a free one, dedup on a grid too small for 27 distinct
    /// neighbours. molpack owns the *occupancy* (the linked list and the empty
    /// -cell bookkeeping below), never a second lattice partition.
    grid: CellGrid,
    /// Minimum image, captured once from `bx` — a read-many `Copy` value, as
    /// [`SimBox::mic`] intends.
    mic: Mic,
    /// Scalar restatements of `bx`'s origin and edge lengths for the per-atom
    /// wrap; `bx` stays the authority, these are derived at construction and
    /// never written again.
    origin: [F; 3],
    length: [F; 3],
    /// Per cell: first occupied slot, or `-1`.
    head: Vec<i32>,
    /// Per slot: next slot in the same cell, or `-1`.
    next: Vec<i32>,
    /// Per slot: the cell it sits in, or `-1` when not placed.
    cell_of: Vec<i32>,
    pos: Vec<[F; 3]>,
    radius: Vec<F>,
    mol_of: Vec<u32>,
    atom_of: Vec<u32>,
    /// Per cell: number of placed slots.
    cell_count: Vec<u32>,
    /// Flat indices of currently empty cells, unordered (swap-remove).
    empty_cells: Vec<u32>,
    /// Per cell: its position in `empty_cells`, or `-1` when occupied.
    empty_pos: Vec<i32>,
}

/// Outcome of probing a candidate position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Probe {
    /// Refused: some non-excluded neighbour is inside the hard core.
    Blocked,
    /// Allowed, with a dimensionless crowding penalty (`0.0` = clear).
    Room(F),
}

/// Why a candidate sat inside another atom's hard core.
///
/// [`OverlapField::probe`] still returns the unit variant [`Probe::Blocked`]
/// on the first such neighbour and does not produce this classification —
/// call [`OverlapField::block_kind`] afterwards.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockKind {
    /// A same-molecule, non-excluded neighbour sits inside the hard core.
    SelfBlocked,
    /// A different-molecule neighbour sits inside the hard core.
    InterChain,
}

impl OverlapField {
    /// Build a field over `origin .. origin + length`.
    ///
    /// `cutoff` must be at least the largest hard core plus the soft shell:
    /// the cell edge is chosen so the 27-cell stencil covers it.
    pub fn new(
        origin: [F; 3],
        length: [F; 3],
        pbc: [bool; 3],
        radius: Vec<F>,
        mol_of: Vec<u32>,
        atom_of: Vec<u32>,
        cutoff: F,
    ) -> Self {
        let slots = radius.len();
        debug_assert_eq!(slots, mol_of.len());
        debug_assert_eq!(slots, atom_of.len());
        let bx = SimBox::ortho(
            ndarray::array![length[0], length[1], length[2]],
            ndarray::array![origin[0], origin[1], origin[2]],
            pbc,
        )
        .expect("the growth box is validated by the entry before a field is built");
        let mic = bx.mic();
        let grid = CellGrid::for_cutoff(&bx, cutoff);
        let ncells = grid.n_cells();
        Self {
            bx,
            grid,
            mic,
            origin,
            length,
            head: vec![-1; ncells],
            next: vec![-1; slots],
            cell_of: vec![-1; slots],
            pos: vec![[0.0; 3]; slots],
            radius,
            mol_of,
            atom_of,
            cell_count: vec![0; ncells],
            empty_cells: (0..ncells as u32).collect(),
            empty_pos: (0..ncells as i32).collect(),
        }
    }

    /// How many slots are placed (test bookkeeping check).
    #[cfg(test)]
    pub fn n_placed(&self) -> usize {
        self.cell_of.iter().filter(|&&c| c >= 0).count()
    }

    /// Position of a placed slot, wrapped into the box.
    #[cfg(test)]
    pub fn position(&self, slot: usize) -> [F; 3] {
        self.pos[slot]
    }

    pub fn is_placed(&self, slot: usize) -> bool {
        self.cell_of[slot] >= 0
    }

    /// Wrap a point into the box along every periodic axis
    /// ([`SimBox::wrap_row`]: single point, no allocation on this hot path).
    #[inline]
    pub fn wrap(&self, p: [F; 3]) -> [F; 3] {
        self.bx.wrap_row(p)
    }

    /// The cell holding `p`, then its stencil neighbours written into `buf`.
    ///
    /// Returned as `(icell, n_neighbours)`: a one-sided query reads `icell`
    /// first and then `buf[..n]`. Both come from [`CellGrid`], which clamps on
    /// a free axis rather than aliasing to the opposite face (an out-of-box
    /// atom filed on the far side would be invisible to every probe near its
    /// true position, silently voiding the constructive guarantee) and dedups
    /// a grid with fewer than three cells on an axis.
    #[inline]
    fn cells_around(&self, p: [F; 3], buf: &mut [usize; 27]) -> (usize, usize) {
        let icell = self.grid.cell_of(&self.bx, p);
        (icell, self.grid.stencil_all(icell, buf))
    }

    /// Insert `slot` at `p` (wrapped on insertion).
    pub fn insert(&mut self, slot: usize, p: [F; 3]) {
        debug_assert!(self.cell_of[slot] < 0, "slot already placed");
        let p = self.wrap(p);
        let cell = self.grid.cell_of(&self.bx, p);
        self.pos[slot] = p;
        self.next[slot] = self.head[cell];
        self.head[cell] = slot as i32;
        self.cell_of[slot] = cell as i32;
        self.cell_count[cell] += 1;
        if self.cell_count[cell] == 1 {
            let i = self.empty_pos[cell] as usize;
            let last = *self.empty_cells.last().expect("cell was in the empty list");
            self.empty_cells[i] = last;
            self.empty_pos[last as usize] = i as i32;
            self.empty_cells.pop();
            self.empty_pos[cell] = -1;
        }
    }

    /// Remove a previously inserted slot. A no-op on an empty slot.
    pub fn remove(&mut self, slot: usize) {
        let cell = self.cell_of[slot];
        if cell < 0 {
            return;
        }
        let cell = cell as usize;
        let mut cur = self.head[cell];
        let mut prev: i32 = -1;
        while cur >= 0 {
            if cur as usize == slot {
                if prev < 0 {
                    self.head[cell] = self.next[slot];
                } else {
                    self.next[prev as usize] = self.next[slot];
                }
                break;
            }
            prev = cur;
            cur = self.next[cur as usize];
        }
        self.next[slot] = -1;
        self.cell_of[slot] = -1;
        self.cell_count[cell] -= 1;
        if self.cell_count[cell] == 0 {
            self.empty_pos[cell] = self.empty_cells.len() as i32;
            self.empty_cells.push(cell as u32);
        }
    }

    /// Number of cells currently holding no atom.
    #[cfg(test)]
    pub(crate) fn n_empty_cells(&self) -> usize {
        self.empty_cells.len()
    }

    /// A point inside a uniformly chosen *empty* cell (Mezei-style cavity
    /// seeding), or `None` when every cell is occupied.
    ///
    /// Randomness is the caller's: `u[0]` picks the cell, `u[1..4]` the
    /// point inside it, so the RNG stream stays in the growth driver and
    /// consumption is fixed per call. An empty cell can still refuse a
    /// candidate (a neighbour cell's atom may reach into it) — the caller's
    /// probe stays the arbiter; this only steers trials toward voids.
    pub fn empty_cell_point(&self, u: [F; 4]) -> Option<[F; 3]> {
        if self.empty_cells.is_empty() {
            return None;
        }
        let pick = ((u[0] * self.empty_cells.len() as F) as usize).min(self.empty_cells.len() - 1);
        let c = self.grid.unflat(self.empty_cells[pick] as usize);
        let celldim = self.grid.celldim();
        Some(std::array::from_fn(|k| {
            self.origin[k] + (c[k] as F + u[k + 1]) * self.length[k] / celldim[k] as F
        }))
    }

    /// Score a candidate position for `slot`.
    ///
    /// `p` is a lab-frame position in Å. `excluded` lists same-molecule
    /// template atom indices held at template geometry (sorted). `hard_scale`
    /// is dimensionless: `1.0` is full declared contact (`radius_i + radius_j`);
    /// the growth driver walks it down by
    /// [`GrowConfig::SOFTEN_RUNG`](crate::grow::config::GrowConfig::SOFTEN_RUNG) per softening rung to
    /// `min_hard_scale` (default 0.8, Auhl's 0.8 × σ floor, where σ is the
    /// excluded-volume / bead diameter) only when a placement would otherwise
    /// dead-end, and reports how often it had to. For a different-molecule
    /// pair, `soft_shell` is an extra width in Å added to that *scaled*
    /// contact: `contact < d < contact + soft_shell` is allowed but charged.
    /// Same-molecule non-excluded pairs get the hard core only (no soft
    /// charge) — intramolecular statistics belong to the priors.
    /// Accumulation of that dimensionless crowding penalty stops once it
    /// passes `penalty_cap`, since the configurational-bias (Rosenbluth)
    /// weight is already negligible there.
    pub fn probe(
        &self,
        slot: usize,
        p: [F; 3],
        excluded: &[u32],
        hard_scale: F,
        soft_shell: F,
        penalty_cap: F,
    ) -> Probe {
        let p = self.wrap(p);
        let mol = self.mol_of[slot];
        let mut penalty = 0.0;
        let mut buf = [0usize; 27];
        let (icell, n_near) = self.cells_around(p, &mut buf);
        for &cell in std::iter::once(&icell).chain(buf[..n_near].iter()) {
            let mut cur = self.head[cell];
            while cur >= 0 {
                let q = cur as usize;
                cur = self.next[q];
                // Skip the probing slot's own (old) position, like `nearest`
                // does — a relax pass may score a new candidate before
                // removing the old placement.
                if self.skip_neighbour(slot, q, excluded) {
                    continue;
                }
                let d2 = self.dist2(p, self.pos[q]);
                let contact = self.hard_contact(slot, q, hard_scale);
                if d2 < contact * contact {
                    return Probe::Blocked;
                }
                // The soft shell is an *inter*-molecular crowding bias.
                // Same-molecule non-excluded pairs get the hard core only (no
                // self-threading), never the soft charge: intrachain
                // conformational statistics belong to the priors, and charging
                // them here would double-count sterics and systematically
                // stretch the chains.
                if self.mol_of[q] == mol {
                    continue;
                }
                let shell = contact + soft_shell;
                let shell2 = shell * shell;
                if d2 < shell2 {
                    let gap = (shell2 - d2) / shell2;
                    penalty += gap * gap;
                    if penalty > penalty_cap {
                        // Clamp so equally-hopeless candidates compare equal
                        // instead of by list-traversal order.
                        return Probe::Room(penalty_cap);
                    }
                }
            }
        }
        Probe::Room(penalty)
    }

    /// Classify a hard-core refusal at `p` for `slot`.
    ///
    /// Post-failure query: call after [`Self::probe`] returns [`Probe::Blocked`].
    /// Same skip and scaled-core test as `probe`'s hard-core arm. `p` is a
    /// lab-frame position in Å. `hard_scale` is dimensionless (`1.0` = full
    /// declared contact `radius_i + radius_j`). The soft shell is ignored.
    ///
    /// Same-molecule hits win: the first [`BlockKind::SelfBlocked`] returns
    /// immediately; [`BlockKind::InterChain`] is remembered until the 27-cell
    /// neighbourhood is exhausted. `probe` does not call this and stays
    /// first-hit (returns `Blocked` on the first core neighbour, without a
    /// kind).
    ///
    /// Returns `None` when no non-excluded neighbour sits inside the scaled
    /// core (empty field, excluded-only neighbours, or every hit outside
    /// contact).
    pub fn block_kind(
        &self,
        slot: usize,
        p: [F; 3],
        excluded: &[u32],
        hard_scale: F,
    ) -> Option<BlockKind> {
        let p = self.wrap(p);
        let mut inter = None;
        let mut buf = [0usize; 27];
        let (icell, n_near) = self.cells_around(p, &mut buf);
        for &cell in std::iter::once(&icell).chain(buf[..n_near].iter()) {
            let mut cur = self.head[cell];
            while cur >= 0 {
                let q = cur as usize;
                cur = self.next[q];
                match self.hard_core_kind(slot, q, p, excluded, hard_scale) {
                    Some(BlockKind::SelfBlocked) => return Some(BlockKind::SelfBlocked),
                    Some(BlockKind::InterChain) => inter = Some(BlockKind::InterChain),
                    None => {}
                }
            }
        }
        inter
    }

    /// Smallest non-excluded distance from `p` to any placed atom.
    pub fn nearest(&self, slot: usize, p: [F; 3], excluded: &[u32]) -> F {
        let p = self.wrap(p);
        let mut best = F::INFINITY;
        let mut buf = [0usize; 27];
        let (icell, n_near) = self.cells_around(p, &mut buf);
        for &cell in std::iter::once(&icell).chain(buf[..n_near].iter()) {
            let mut cur = self.head[cell];
            while cur >= 0 {
                let q = cur as usize;
                cur = self.next[q];
                if self.skip_neighbour(slot, q, excluded) {
                    continue;
                }
                let d2 = self.dist2(p, self.pos[q]);
                if d2 < best {
                    best = d2;
                }
            }
        }
        best.sqrt()
    }

    /// Self-slot and same-molecule excluded template atoms — shared by
    /// `probe`, `nearest`, and `block_kind`.
    #[inline]
    fn skip_neighbour(&self, slot: usize, q: usize, excluded: &[u32]) -> bool {
        q == slot
            || (self.mol_of[q] == self.mol_of[slot]
                && excluded.binary_search(&self.atom_of[q]).is_ok())
    }

    /// Scaled hard-core contact distance (Å) for the (`slot`, `q`) pair.
    /// `hard_scale` is dimensionless (`1.0` = full declared contact
    /// `radius_i + radius_j`).
    #[inline]
    fn hard_contact(&self, slot: usize, q: usize, hard_scale: F) -> F {
        (self.radius[slot] + self.radius[q]) * hard_scale
    }

    /// Hard-core hit kind, or `None` if skipped or outside the scaled core.
    #[inline]
    fn hard_core_kind(
        &self,
        slot: usize,
        q: usize,
        p: [F; 3],
        excluded: &[u32],
        hard_scale: F,
    ) -> Option<BlockKind> {
        if self.skip_neighbour(slot, q, excluded) {
            return None;
        }
        let contact = self.hard_contact(slot, q, hard_scale);
        if self.dist2(p, self.pos[q]) < contact * contact {
            if self.mol_of[q] == self.mol_of[slot] {
                Some(BlockKind::SelfBlocked)
            } else {
                Some(BlockKind::InterChain)
            }
        } else {
            None
        }
    }

    #[inline]
    fn dist2(&self, a: [F; 3], b: [F; 3]) -> F {
        let d = self.mic.apply([a[0] - b[0], a[1] - b[1], a[2] - b[2]]);
        d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    }
}
