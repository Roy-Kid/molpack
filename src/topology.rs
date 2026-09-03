//! The template bond graph, read once and shared.
//!
//! A *template* is the single reference copy of a molecular species: the atom
//! positions and the bond list that every packed copy of that species is a
//! transform of. It reaches this crate as a `molrs` [`Frame`] — a container of
//! named blocks, each block a set of typed columns. Two blocks matter here:
//! `atoms`, whose float columns `x` / `y` / `z` hold one coordinate triple per
//! atom (in ångström, Å — the crate's length unit), and `bonds`, whose
//! unsigned columns `atomi` / `atomj` hold, per bond, the two 0-based atom
//! indices it joins.
//!
//! A template is therefore geometry plus connectivity. This leaf owns the
//! connectivity half: [`Topology`] reads a [`Frame`] into an atom count, the
//! bond list as written, and a compressed-sparse-row (CSR) adjacency — every
//! atom's neighbor list concatenated into one flat `neighbors` array, plus an
//! `offsets` array whose entry `a` records where atom `a`'s neighbors begin,
//! so `neighbors[offsets[a]..offsets[a + 1]]` is that atom's slice and the
//! whole graph costs two allocations instead of one per atom.
//! [`frame_positions`] reads the geometry half from the same `atoms` block, so
//! both halves refuse the same unreadable frames with the same
//! [`TopologyError`].
//!
//! Two representation choices here are load-bearing for growth and must not be
//! "tidied":
//!
//! - **Adjacency follows bond-file insertion order, never sorted.** The CSR
//!   arrays are the flattening of `adj[i].push(j); adj[j].push(i)` walked over
//!   the bonds block in file order. Growth's tree decomposition takes the
//!   *first* neighbor of an atom (`pick_ref` / `bfs_order` in
//!   `src/grow/internal.rs`), so sorting the neighbor slices would reshape the
//!   growth tree and move every grown coordinate.
//! - **[`Topology::exclusions`] includes the root atom.** Each list holds every
//!   atom reachable from its root in at most `depth` bonds — the ball filled
//!   out by a breadth-first search (BFS: visit the root, then all its
//!   neighbors, then their unvisited neighbors, and so on, so every atom is
//!   first reached along a shortest bond path) — root included, sorted
//!   ascending. Growth hands these lists to the overlap field as skip sets; a
//!   root-exclusive list would score an atom against itself.
//!
//! This module is a crate-root leaf: it depends on `molrs` and `std` only, so
//! every consumer (growth today, further stages later) shares one bond graph
//! rather than each parsing the frame again.

use std::collections::VecDeque;
use std::fmt;

use molrs::store::frame::Frame;
use molrs::types::F;

/// A template's bond graph: atom count, bonds as written, CSR adjacency.
///
/// Construct with [`Topology::from_frame`]. The value owns its arrays and is
/// cheap to clone; there is no cache and nothing to invalidate.
///
/// ```no_run
/// use molpack::Topology;
///
/// # fn demo(frame: &molrs::store::frame::Frame) -> Result<(), molpack::TopologyError> {
/// let topo = Topology::from_frame(frame)?;
/// topo.require_connected()?;
/// let skip = topo.exclusions(3); // 1-2 / 1-3 / 1-4 partners, root included
/// # let _ = skip;
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct Topology {
    /// Number of atoms in the template.
    natoms: usize,
    /// Bonds in bonds-block order, orientation as written, self-loops dropped.
    bonds: Vec<(u32, u32)>,
    /// CSR row starts, `natoms + 1` entries.
    offsets: Vec<u32>,
    /// CSR neighbor entries, per atom in bond-file insertion order.
    neighbors: Vec<u32>,
}

impl Topology {
    /// Read the bond graph of a template frame.
    ///
    /// The atom count comes from [`frame_positions`], so both halves of a
    /// template accept and refuse exactly the same frames. Connectivity is
    /// *not* checked here — see [`Topology::require_connected`].
    ///
    /// This leaf never reports an atom-count error, because it does not know
    /// its consumer's minimum (growth needs three atoms for its rigid seed;
    /// another stage may be happy with one). A bondless two-atom template is
    /// therefore [`TopologyError::NoBonds`] here, and the "too small" verdict
    /// is the consumer's to add.
    ///
    /// # Errors
    ///
    /// - [`TopologyError::NoAtomsBlock`] — the frame has no `atoms` block with
    ///   `x` / `y` / `z` float columns, exactly as in [`frame_positions`].
    /// - [`TopologyError::BondOutOfRange`] — a bond names an endpoint outside
    ///   the template's atom range.
    /// - [`TopologyError::NoBonds`] — the `bonds` block is missing, lacks the
    ///   `atomi` / `atomj` unsigned columns, is empty, or holds nothing but
    ///   self-loops (a self-loop `(a, a)` is dropped silently).
    pub fn from_frame(frame: &Frame) -> Result<Topology, TopologyError> {
        Self::from_frame_with_positions(frame).map(|(topo, _xyz)| topo)
    }

    /// The same read as [`Topology::from_frame`], also handing back the
    /// coordinates it had to read to learn the atom count.
    ///
    /// The coordinates are in Å (the crate's length unit, as carried by the
    /// frame), in atom order. A consumer that needs both halves of a template
    /// (growth's internal-coordinate decomposition, the lattice backbone
    /// analysis) takes this one and walks the `atoms` block once.
    ///
    /// # Errors
    ///
    /// The same as [`Topology::from_frame`], raised in the same order.
    pub fn from_frame_with_positions(
        frame: &Frame,
    ) -> Result<(Topology, Vec<[F; 3]>), TopologyError> {
        // The coordinates themselves belong to the caller's geometry pass; what
        // the bond graph needs from them is the atom count, under the same
        // readability rule, so that both halves reject the same frames.
        let xyz = frame_positions(frame)?;
        let natoms = xyz.len();
        let bonds = frame_bonds(frame, natoms)?;
        let (offsets, neighbors) = build_adjacency(natoms, &bonds);
        let topo = Topology {
            natoms,
            bonds,
            offsets,
            neighbors,
        };
        Ok((topo, xyz))
    }

    /// Number of atoms in the template.
    pub fn natoms(&self) -> usize {
        self.natoms
    }

    /// The bonds in bonds-block order, each keeping the orientation it was
    /// written with (a ring closure stays `(5, 0)`; it is not normalized).
    /// Self-loops are already dropped.
    ///
    /// Each entry is a pair of 0-based indices into the template's atom order
    /// — the same order [`frame_positions`] returns coordinates in.
    pub fn bonds(&self) -> &[(u32, u32)] {
        &self.bonds
    }

    /// The neighbors of `atom`, in the order the bonds block introduced them.
    ///
    /// `atom` is the 0-based index of the atom in the template's atom order,
    /// and the returned entries are the 0-based indices of the atoms bonded to
    /// it. Never sorted: the first entry is the neighbor growth's tree walk
    /// picks as its reference atom (see the module docs).
    ///
    /// # Panics
    ///
    /// If `atom >= self.natoms()`.
    pub fn neighbors(&self, atom: usize) -> &[u32] {
        let start = self.offsets[atom] as usize;
        let end = self.offsets[atom + 1] as usize;
        &self.neighbors[start..end]
    }

    /// Reject a template in which some atom carries no bond at all (graph
    /// degree zero — the bonds block never names it).
    ///
    /// This is the cheap screen only. A template that splits into two *bonded*
    /// components passes here and is caught by the consumer's own traversal
    /// (growth's breadth-first search in `src/grow/internal.rs`), which reports
    /// the same [`TopologyError::Disconnected`].
    ///
    /// # Errors
    ///
    /// [`TopologyError::Disconnected`] if any atom has degree zero.
    pub fn require_connected(&self) -> Result<(), TopologyError> {
        if (0..self.natoms).any(|atom| self.neighbors(atom).is_empty()) {
            return Err(TopologyError::Disconnected);
        }
        Ok(())
    }

    /// Per-atom same-molecule partners within `depth` bonds (sorted, includes
    /// the root atom).
    ///
    /// Entry `a` lists every atom within `depth` bonds of atom `a` — the
    /// breadth-first-search ball of that radius — ascending, with `a` itself as
    /// a member. The partners are named by their bond distance: a 1-2 partner
    /// is a bonded neighbor, a 1-3 partner shares an angle with the root, and a
    /// 1-4 partner is three bonds away and swings with the torsion between
    /// them — so `depth = 3`, the all-atom default, exempts exactly the pairs
    /// whose distance the template's bonds, angles and torsion prior already
    /// fix. Callers use these lists as overlap-field skip sets, where the root
    /// must be skipped too.
    pub fn exclusions(&self, depth: usize) -> Vec<Vec<u32>> {
        let n = self.natoms;
        let mut out = Vec::with_capacity(n);
        let mut dist = vec![usize::MAX; n];
        for root in 0..n {
            let mut touched = vec![root];
            dist[root] = 0;
            let mut queue = VecDeque::new();
            queue.push_back(root);
            while let Some(a) = queue.pop_front() {
                if dist[a] == depth {
                    continue;
                }
                for &b in self.neighbors(a) {
                    let b = b as usize;
                    if dist[b] == usize::MAX {
                        dist[b] = dist[a] + 1;
                        touched.push(b);
                        queue.push_back(b);
                    }
                }
            }
            let mut list: Vec<u32> = touched.iter().map(|&a| a as u32).collect();
            list.sort_unstable();
            out.push(list);
            // `dist` is reused across roots; only the touched entries are dirty.
            for a in touched {
                dist[a] = usize::MAX;
            }
        }
        out
    }
}

/// Read a template frame's coordinates, in atom order.
///
/// The returned triples are coordinates in Å (the crate's length unit, as
/// carried by the frame — nothing is converted), one per atom, in the `atoms`
/// block's row order. This is the geometry counterpart of
/// [`Topology::from_frame`] and shares its error type, so a consumer reading
/// both halves of a template reports one kind of failure.
///
/// # Errors
///
/// [`TopologyError::NoAtomsBlock`] if the frame carries no `atoms` block, or
/// that block is missing any of the `x`, `y`, `z` float columns.
pub fn frame_positions(frame: &Frame) -> Result<Vec<[F; 3]>, TopologyError> {
    let atoms = frame.get("atoms").ok_or(TopologyError::NoAtomsBlock)?;
    let x = atoms.get_float("x").ok_or(TopologyError::NoAtomsBlock)?;
    let y = atoms.get_float("y").ok_or(TopologyError::NoAtomsBlock)?;
    let z = atoms.get_float("z").ok_or(TopologyError::NoAtomsBlock)?;
    Ok(x.iter()
        .zip(y.iter())
        .zip(z.iter())
        .map(|((&a, &b), &c)| [a, b, c])
        .collect())
}

/// The `bonds` block as a bond list, validated against `n` atoms.
fn frame_bonds(frame: &Frame, n: usize) -> Result<Vec<(u32, u32)>, TopologyError> {
    let block = frame.get("bonds").ok_or(TopologyError::NoBonds)?;
    let i = block.get_uint("atomi").ok_or(TopologyError::NoBonds)?;
    let j = block.get_uint("atomj").ok_or(TopologyError::NoBonds)?;
    let mut out = Vec::with_capacity(i.len());
    for (&a, &b) in i.iter().zip(j.iter()) {
        let (a, b) = (a as usize, b as usize);
        if a >= n || b >= n {
            return Err(TopologyError::BondOutOfRange { a, b, n });
        }
        if a != b {
            out.push((a as u32, b as u32));
        }
    }
    if out.is_empty() {
        return Err(TopologyError::NoBonds);
    }
    Ok(out)
}

/// CSR adjacency of `bonds`, filled in bond-file order so each atom's slice
/// replays `adj[i].push(j); adj[j].push(i)` — see the module docs on why this
/// order is never sorted away.
fn build_adjacency(n: usize, bonds: &[(u32, u32)]) -> (Vec<u32>, Vec<u32>) {
    let mut offsets = vec![0u32; n + 1];
    for &(a, b) in bonds {
        offsets[a as usize + 1] += 1;
        offsets[b as usize + 1] += 1;
    }
    for atom in 0..n {
        offsets[atom + 1] += offsets[atom];
    }
    let mut cursor = offsets.clone();
    let mut neighbors = vec![0u32; 2 * bonds.len()];
    for &(a, b) in bonds {
        neighbors[cursor[a as usize] as usize] = b;
        cursor[a as usize] += 1;
        neighbors[cursor[b as usize] as usize] = a;
        cursor[b as usize] += 1;
    }
    (offsets, neighbors)
}

/// Why a template's bond graph cannot be read, or does not qualify.
///
/// Owned by this leaf so that consumers wrap it (`GrowError::Topology(..)`)
/// instead of restating its variants; the message text is what the user sees.
#[derive(Debug, Clone)]
pub enum TopologyError {
    /// The template frame has no readable `atoms` block.
    NoAtomsBlock,
    /// The template frame has no `bonds` block (or an empty one).
    NoBonds,
    /// A bond references an atom index outside the template.
    BondOutOfRange {
        /// First endpoint of the offending bond.
        a: usize,
        /// Second endpoint of the offending bond.
        b: usize,
        /// Atom count of the template.
        n: usize,
    },
    /// The template's bond graph does not connect all atoms.
    Disconnected,
}

impl fmt::Display for TopologyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TopologyError::NoAtomsBlock => {
                write!(f, "the template frame has no readable atoms block")
            }
            TopologyError::NoBonds => write!(
                f,
                "the template frame carries no bonds; growth needs the bond graph — pack \
                 this target with GenCanPack or supply connectivity"
            ),
            TopologyError::BondOutOfRange { a, b, n } => write!(
                f,
                "bond ({a}, {b}) references an atom outside the template (natoms = {n})"
            ),
            TopologyError::Disconnected => {
                write!(f, "the template's bond graph does not connect all atoms")
            }
        }
    }
}

impl std::error::Error for TopologyError {}
