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

use molrs::types::F;

/// A placed-atom field over an orthorhombic, optionally periodic box.
///
/// Atoms live in a flat *slot* space so one field can carry several species at
/// once: slot `s` belongs to molecule `mol_of[s]` and is that molecule's
/// template atom `atom_of[s]`.
#[derive(Debug, Clone)]
pub struct OverlapField {
    origin: [F; 3],
    length: [F; 3],
    pbc: [bool; 3],
    nc: [usize; 3],
    inv_cell: [F; 3],
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
    n_placed: usize,
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
        let nc: [usize; 3] =
            std::array::from_fn(|k| ((length[k] / cutoff).floor() as usize).max(1));
        let inv_cell = std::array::from_fn(|k| nc[k] as F / length[k]);
        let ncells = nc[0] * nc[1] * nc[2];
        Self {
            origin,
            length,
            pbc,
            nc,
            inv_cell,
            head: vec![-1; ncells],
            next: vec![-1; slots],
            cell_of: vec![-1; slots],
            pos: vec![[0.0; 3]; slots],
            radius,
            mol_of,
            atom_of,
            n_placed: 0,
            cell_count: vec![0; ncells],
            empty_cells: (0..ncells as u32).collect(),
            empty_pos: (0..ncells as i32).collect(),
        }
    }

    pub fn n_placed(&self) -> usize {
        self.n_placed
    }

    /// Position of a placed slot, wrapped into the box.
    pub fn position(&self, slot: usize) -> [F; 3] {
        self.pos[slot]
    }

    pub fn is_placed(&self, slot: usize) -> bool {
        self.cell_of[slot] >= 0
    }

    /// Wrap a point into the box along every periodic axis.
    #[inline]
    pub fn wrap(&self, p: [F; 3]) -> [F; 3] {
        std::array::from_fn(|k| {
            if self.pbc[k] {
                self.origin[k] + (p[k] - self.origin[k]).rem_euclid(self.length[k])
            } else {
                p[k]
            }
        })
    }

    #[inline]
    fn cell_index(&self, p: [F; 3]) -> usize {
        let c: [usize; 3] = std::array::from_fn(|k| {
            let f = ((p[k] - self.origin[k]) * self.inv_cell[k]).floor() as isize;
            if self.pbc[k] {
                f.rem_euclid(self.nc[k] as isize) as usize
            } else {
                // A non-periodic axis must clamp, not alias: wrapping would
                // file an out-of-box atom in the cell on the opposite face,
                // where no probe stencil near its true position would ever
                // find it — silently voiding the constructive guarantee.
                f.clamp(0, self.nc[k] as isize - 1) as usize
            }
        });
        (c[2] * self.nc[1] + c[1]) * self.nc[0] + c[0]
    }

    /// Insert `slot` at `p` (wrapped on insertion).
    pub fn insert(&mut self, slot: usize, p: [F; 3]) {
        debug_assert!(self.cell_of[slot] < 0, "slot already placed");
        let p = self.wrap(p);
        let cell = self.cell_index(p);
        self.pos[slot] = p;
        self.next[slot] = self.head[cell];
        self.head[cell] = slot as i32;
        self.cell_of[slot] = cell as i32;
        self.n_placed += 1;
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
        self.n_placed -= 1;
        self.cell_count[cell] -= 1;
        if self.cell_count[cell] == 0 {
            self.empty_pos[cell] = self.empty_cells.len() as i32;
            self.empty_cells.push(cell as u32);
        }
    }

    /// Number of cells currently holding no atom.
    pub fn n_empty_cells(&self) -> usize {
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
        let c = self.unflat(self.empty_cells[pick] as usize);
        Some(std::array::from_fn(|k| {
            self.origin[k] + (c[k] as F + u[k + 1]) * self.length[k] / self.nc[k] as F
        }))
    }

    /// Score a candidate position for `slot`.
    ///
    /// `p` is a lab-frame position in Å. `excluded` lists same-molecule
    /// template atom indices held at template geometry (sorted). `hard_scale`
    /// is dimensionless: `1.0` is full declared contact (`radius_i + radius_j`);
    /// the growth driver walks it down by 0.97 per softening rung to
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
        let c0 = self.unflat(self.cell_index(p));
        let mut seen = [usize::MAX; 27];
        let mut n_seen = 0usize;
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let Some(cell) = self.shift(c0, [dx, dy, dz]) else {
                        continue;
                    };
                    if seen[..n_seen].contains(&cell) {
                        continue;
                    }
                    seen[n_seen] = cell;
                    n_seen += 1;
                    let mut cur = self.head[cell];
                    while cur >= 0 {
                        let q = cur as usize;
                        cur = self.next[q];
                        // Skip the probing slot's own (old) position, like
                        // `nearest` does — a relax pass may score a new
                        // candidate before removing the old placement.
                        if self.skip_neighbour(slot, q, excluded) {
                            continue;
                        }
                        let d2 = self.dist2(p, self.pos[q]);
                        let contact = self.hard_contact(slot, q, hard_scale);
                        if d2 < contact * contact {
                            return Probe::Blocked;
                        }
                        // The soft shell is an *inter*-molecular crowding
                        // bias. Same-molecule non-excluded pairs get the
                        // hard core only (no self-threading), never the soft
                        // charge: intrachain conformational statistics belong
                        // to the priors, and charging them here would
                        // double-count sterics and systematically stretch
                        // the chains.
                        if self.mol_of[q] == mol {
                            continue;
                        }
                        let shell = contact + soft_shell;
                        let shell2 = shell * shell;
                        if d2 < shell2 {
                            let gap = (shell2 - d2) / shell2;
                            penalty += gap * gap;
                            if penalty > penalty_cap {
                                // Clamp so equally-hopeless candidates compare
                                // equal instead of by list-traversal order.
                                return Probe::Room(penalty_cap);
                            }
                        }
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
        let c0 = self.unflat(self.cell_index(p));
        let mut seen = [usize::MAX; 27];
        let mut n_seen = 0usize;
        let mut inter = None;
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let Some(cell) = self.shift(c0, [dx, dy, dz]) else {
                        continue;
                    };
                    if seen[..n_seen].contains(&cell) {
                        continue;
                    }
                    seen[n_seen] = cell;
                    n_seen += 1;
                    let mut cur = self.head[cell];
                    while cur >= 0 {
                        let q = cur as usize;
                        cur = self.next[q];
                        match self.hard_core_kind(slot, q, p, excluded, hard_scale) {
                            Some(BlockKind::SelfBlocked) => {
                                return Some(BlockKind::SelfBlocked);
                            }
                            Some(BlockKind::InterChain) => inter = Some(BlockKind::InterChain),
                            None => {}
                        }
                    }
                }
            }
        }
        inter
    }

    /// Smallest non-excluded distance from `p` to any placed atom.
    pub fn nearest(&self, slot: usize, p: [F; 3], excluded: &[u32]) -> F {
        let p = self.wrap(p);
        let c0 = self.unflat(self.cell_index(p));
        let mut best = F::INFINITY;
        let mut seen = [usize::MAX; 27];
        let mut n_seen = 0usize;
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let Some(cell) = self.shift(c0, [dx, dy, dz]) else {
                        continue;
                    };
                    if seen[..n_seen].contains(&cell) {
                        continue;
                    }
                    seen[n_seen] = cell;
                    n_seen += 1;
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
        let mut s = 0.0;
        for k in 0..3 {
            let mut d = a[k] - b[k];
            if self.pbc[k] {
                d -= (d / self.length[k]).round() * self.length[k];
            }
            s += d * d;
        }
        s
    }

    #[inline]
    fn unflat(&self, cell: usize) -> [usize; 3] {
        let x = cell % self.nc[0];
        let y = (cell / self.nc[0]) % self.nc[1];
        let z = cell / (self.nc[0] * self.nc[1]);
        [x, y, z]
    }

    /// Neighbour cell at `c + d`, or `None` when it falls outside a
    /// non-periodic axis.
    #[inline]
    fn shift(&self, c: [usize; 3], d: [i64; 3]) -> Option<usize> {
        let mut out = [0usize; 3];
        for k in 0..3 {
            let n = self.nc[k] as i64;
            let v = c[k] as i64 + d[k];
            out[k] = if self.pbc[k] {
                v.rem_euclid(n) as usize
            } else if v < 0 || v >= n {
                return None;
            } else {
                v as usize
            };
        }
        Some((out[2] * self.nc[1] + out[1]) * self.nc[0] + out[0])
    }
}
