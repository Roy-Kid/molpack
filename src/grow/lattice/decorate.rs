//! Decoration: an on-lattice walk becomes a continuous all-atom chain.
//!
//! The lattice decides only the TORSION SEQUENCE; every bond length, bond
//! angle and rigid torsion is rebuilt verbatim from the template through the
//! existing [`InternalTree`] recipes, so the decorated molecule carries the
//! template's bonded geometry exactly. Because template bonds differ from
//! the uniform lattice bond (PEO C–O 1.41 Å vs C–C 1.53 Å), the decorated
//! chain drifts off its lattice path with chain length — the lattice's
//! excluded-volume guarantee degrades accordingly, and the shared objective
//! (plus the seeded GENCAN push-off) is the honest backstop
//! (lattice-growth-phase spec, risk 3).

use molrs::store::frame::Frame;
use molrs::types::F;

use crate::grow::GrowError;
use crate::grow::internal::{InternalTree, dihedral, nerf};
use crate::grow::lattice::saw::DiamondLattice;

const PI: F = std::f64::consts::PI as F;
const TWO_PI: F = std::f64::consts::TAU as F;

fn wrap_pi(x: F) -> F {
    let mut v = x % TWO_PI;
    if v > PI {
        v -= TWO_PI;
    } else if v <= -PI {
        v += TWO_PI;
    }
    v
}

/// The template's heavy-atom backbone, oriented root-end first, plus the
/// per-backbone-atom hook into the tree's steps.
pub(crate) struct Backbone {
    /// Template atom indices of the backbone, in walk order.
    pub(crate) atoms: Vec<usize>,
    /// For backbone position `j ≥ 3`: `(step, local site index, var, offset)`
    /// of the site that places `atoms[j]`.
    hooks: Vec<Option<(usize, usize, usize, F)>>,
    /// Backbone bond lengths in walk order (Å).
    pub(crate) bonds: Vec<F>,
    /// Mean backbone bond length in the template (Å) — sets the lattice
    /// constant.
    pub(crate) mean_bond: F,
}

fn is_hydrogen(el: &str) -> bool {
    el.eq_ignore_ascii_case("h")
}

/// Extract and validate the backbone of `frame` against `tree`.
///
/// v1 scope (lattice-growth-phase spec): a linear sp³ heavy-atom chain.
/// Branched heavy atoms and non-rotatable interior backbone bonds are
/// named rejections, never silent degradation.
pub(crate) fn analyze_backbone(frame: &Frame, tree: &InternalTree) -> Result<Backbone, GrowError> {
    let (topo, xyz) = crate::grow::topology_for_growth(frame)?;
    let n = topo.n_atoms();
    let elements: Vec<String> = frame
        .get("atoms")
        .and_then(|b| b.get_string("element"))
        .map(|c| c.iter().cloned().collect())
        .unwrap_or_else(|| vec!["X".to_string(); n]);

    let heavy: Vec<bool> = elements.iter().map(|e| !is_hydrogen(e)).collect();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for [i, j] in topo.bonds() {
        if heavy[i] && heavy[j] {
            adj[i].push(j);
            adj[j].push(i);
        }
    }

    let heavies: Vec<usize> = (0..n).filter(|&i| heavy[i]).collect();
    if heavies.len() < 4 {
        return Err(GrowError::NonTetrahedralTemplate(format!(
            "{} heavy atoms; the lattice walk needs a backbone of at least 4",
            heavies.len()
        )));
    }
    for &i in &heavies {
        if adj[i].len() > 2 {
            return Err(GrowError::NonTetrahedralTemplate(format!(
                "heavy atom {i} has {} heavy neighbours — branched templates \
                 are staged for the lattice branch phase",
                adj[i].len()
            )));
        }
        if adj[i].is_empty() {
            return Err(GrowError::NonTetrahedralTemplate(format!(
                "heavy atom {i} is detached from the heavy-atom backbone"
            )));
        }
    }
    let ends: Vec<usize> = heavies
        .iter()
        .copied()
        .filter(|&i| adj[i].len() == 1)
        .collect();
    if ends.len() != 2 {
        return Err(GrowError::NonTetrahedralTemplate(format!(
            "heavy-atom graph has {} endpoints, expected 2 (linear backbone)",
            ends.len()
        )));
    }

    // Walk the path from one end.
    let mut path = vec![ends[0]];
    let mut prev = usize::MAX;
    let mut cur = ends[0];
    while path.len() < heavies.len() {
        let next = *adj[cur]
            .iter()
            .find(|&&q| q != prev)
            .expect("path walk: interior atom has a fresh neighbour");
        path.push(next);
        prev = cur;
        cur = next;
    }

    // Per-atom site lookup: (step k, local index).
    let mut site_of: Vec<Option<(usize, usize)>> = vec![None; n];
    for k in 0..tree.n_steps() {
        for (li, s) in tree.step_sites(k).iter().enumerate() {
            site_of[s.atom] = Some((k, li));
        }
    }

    // Orient the path so that BFS parents run along it: the site of the
    // second atom (whichever end has one) must name the first as refs[0].
    let oriented = |p: &[usize]| -> bool {
        for w in p.windows(2) {
            if let Some((k, li)) = site_of[w[1]] {
                return tree.step_sites(k)[li].refs[0] == w[0];
            }
        }
        true
    };
    let atoms: Vec<usize> = if oriented(&path) {
        path
    } else {
        let rev: Vec<usize> = path.into_iter().rev().collect();
        if !oriented(&rev) {
            return Err(GrowError::NonTetrahedralTemplate(
                "backbone orientation does not follow the tree's BFS parents".to_string(),
            ));
        }
        rev
    };

    // Hooks: for j ≥ 3 the placing site must own a torsion variable and its
    // reference frame must be the backbone itself.
    let mut hooks: Vec<Option<(usize, usize, usize, F)>> = vec![None; atoms.len()];
    for j in 3..atoms.len() {
        let Some((k, li)) = site_of[atoms[j]] else {
            return Err(GrowError::NonTetrahedralTemplate(format!(
                "backbone atom {} sits in the rigid seed too deep in the chain",
                atoms[j]
            )));
        };
        let site = &tree.step_sites(k)[li];
        let Some((var, offset)) = site.follows else {
            return Err(GrowError::NonTetrahedralTemplate(format!(
                "backbone bond into atom {} is not rotatable (non-sp³ backbone?)",
                atoms[j]
            )));
        };
        if site.refs[0] != atoms[j - 1] || site.refs[1] != atoms[j - 2] {
            return Err(GrowError::NonTetrahedralTemplate(format!(
                "backbone atom {}'s reference frame leaves the backbone",
                atoms[j]
            )));
        }
        hooks[j] = Some((k, li, var, offset));
    }

    let bond_lengths: Vec<F> = atoms
        .windows(2)
        .map(|w| {
            let d = [
                xyz[w[0]][0] - xyz[w[1]][0],
                xyz[w[0]][1] - xyz[w[1]][1],
                xyz[w[0]][2] - xyz[w[1]][2],
            ];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .collect();
    let mean_bond = bond_lengths.iter().sum::<F>() / bond_lengths.len() as F;

    Ok(Backbone {
        atoms,
        hooks,
        bonds: bond_lengths,
        mean_bond,
    })
}

/// Build one chain's coordinates from its lattice walk: torsion variables
/// solved so the backbone dihedrals equal the walk's, everything else from
/// the template, then the whole molecule rigidly aligned onto the walk's
/// first three sites.
pub(crate) fn decorate_chain(
    tree: &InternalTree,
    bb: &Backbone,
    lat: &DiamondLattice,
    walk: &[[i64; 3]],
    track_tweak: F,
    coords: &mut [[F; 3]],
) {
    debug_assert_eq!(walk.len(), bb.atoms.len());
    // Targets: the walk's DIRECTION sequence at the chain's own bond
    // metric. Commensuration stretches the lattice by a few percent, so
    // raw site positions are longitudinally unreachable for a chain with
    // template-exact bonds — torsion tracking can only correct
    // perpendicular to the path. Rebuilding the path step-by-step with the
    // template's j-th backbone bond length keeps every target reachable
    // while the dihedral sequence (directions!) stays exactly the lattice's
    // RIS decision.
    let w: Vec<[F; 3]> = {
        let raw: Vec<[F; 3]> = walk.iter().map(|&p| lat.to_continuum(p)).collect();
        let mut out = Vec::with_capacity(raw.len());
        out.push(raw[0]);
        for j in 1..raw.len() {
            let d = [
                raw[j][0] - raw[j - 1][0],
                raw[j][1] - raw[j - 1][1],
                raw[j][2] - raw[j - 1][2],
            ];
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let b = bb.bonds[j - 1] / len;
            let prev: [F; 3] = out[j - 1];
            out.push([prev[0] + b * d[0], prev[1] + b * d[1], prev[2] + b * d[2]]);
        }
        out
    };

    // Backbone position placed by each step (j ≥ 3), if any.
    let mut bb_of_step: Vec<Option<usize>> = vec![None; tree.n_steps()];
    for (j, hook) in bb.hooks.iter().enumerate() {
        if let Some((k, _, _, _)) = hook {
            bb_of_step[*k] = Some(j);
        }
    }

    // Defaults: every variable starts at its template value.
    let mut vars = vec![0.0 as F; tree.n_vars()];
    for k in 0..tree.n_steps() {
        if let Some(v) = tree.step_var(k) {
            vars[v] = tree.template_var(k);
        }
    }

    let sub = |a: [F; 3], b: [F; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [F; 3], b: [F; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let frame_of = |p0: [F; 3], p1: [F; 3], p2: [F; 3]| -> [[F; 3]; 3] {
        let norm = |v: [F; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let e1v = sub(p1, p0);
        let n1 = norm(e1v);
        let e1 = [e1v[0] / n1, e1v[1] / n1, e1v[2] / n1];
        let u = sub(p2, p0);
        let d = dot(u, e1);
        let uo = [u[0] - d * e1[0], u[1] - d * e1[1], u[2] - d * e1[2]];
        let n2 = norm(uo);
        let e2 = [uo[0] / n2, uo[1] / n2, uo[2] / n2];
        let e3 = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        [e1, e2, e3] // rows
    };
    let rot = |p: [F; 3], r: &[[F; 3]; 3]| {
        [
            r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2],
            r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2],
            r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2],
        ]
    };

    // Phase 1: build the seed and every step before the first hooked one at
    // template torsions. That fixes (bb0, bb1, bb2), and with them the rigid
    // build→walk transform — so the walk targets can be brought into the
    // BUILD frame before any tracking decision is made (tracking in the
    // wrong frame chases points ~a box away).
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    tree.place_seed([0.0; 3], &identity, coords);
    let k0 = bb_of_step
        .iter()
        .position(|h| h.is_some())
        .unwrap_or(tree.n_steps());
    for k in 0..k0 {
        tree.place_step(k, &vars, coords);
    }

    let (b0, b1, b2) = (bb.atoms[0], bb.atoms[1], bb.atoms[2]);
    let fc = frame_of(coords[b0], coords[b1], coords[b2]);
    let fw = frame_of(w[0], w[1], w[2]);
    // R = fwᵀ · fc maps build → walk; both frames are row-orthonormal.
    let mut r = [[0.0 as F; 3]; 3];
    for (i, row) in r.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = fw[0][i] * fc[0][j] + fw[1][i] * fc[1][j] + fw[2][i] * fc[2][j];
        }
    }
    let c0 = coords[b0];
    // Targets in the build frame: wb = Rᵀ·(w − w0) + c0.
    let mut rt = [[0.0 as F; 3]; 3];
    for (i, row) in rt.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = r[j][i];
        }
    }
    let wb: Vec<[F; 3]> = w
        .iter()
        .map(|&p| {
            let rd = rot(sub(p, w[0]), &rt);
            [c0[0] + rd[0], c0[1] + rd[1], c0[2] + rd[2]]
        })
        .collect();

    // Phase 2: hooked steps — backbone torsions from the walk, optionally
    // tweaked toward the (build-frame) target sites.
    for (k, hooked) in bb_of_step.iter().enumerate().skip(k0) {
        if let Some(j) = *hooked {
            let (_, li, var, offset) = bb.hooks[j].expect("hooked step");
            let site = &tree.step_sites(k)[li];
            // Trial with site torsion 0: the backbone dihedral is linear in
            // the site torsion with slope 1 (same rotation axis), so one
            // trial solves it exactly.
            let pos0 = nerf(
                coords[site.refs[2]],
                coords[site.refs[1]],
                coords[site.refs[0]],
                site.bond,
                site.angle,
                0.0,
            );
            let beta = dihedral(
                coords[bb.atoms[j - 3]],
                coords[bb.atoms[j - 2]],
                coords[bb.atoms[j - 1]],
                pos0,
            );
            let want = dihedral(wb[j - 3], wb[j - 2], wb[j - 1], wb[j]);
            let phi_lat = want - beta;
            let mut phi = phi_lat;
            if track_tweak > 0.0 {
                // Path tracking: the placed atom sweeps a circle about the
                // (g, p) axis; pick the site torsion closest to this atom's
                // target site, clamped to ±track_tweak around the exact RIS
                // state. Bounded per-step correction ⇒ bounded drift ⇒ the
                // lattice's inter-chain distance guarantee survives
                // decoration.
                let p = coords[site.refs[0]];
                let g = coords[site.refs[1]];
                let q1 = nerf(
                    coords[site.refs[2]],
                    g,
                    p,
                    site.bond,
                    site.angle,
                    std::f64::consts::FRAC_PI_2 as F,
                );
                let axis = sub(p, g);
                let alen = dot(axis, axis).sqrt();
                let u = [axis[0] / alen, axis[1] / alen, axis[2] / alen];
                let d0 = sub(pos0, p);
                let along = dot(d0, u);
                let o = [
                    p[0] + along * u[0],
                    p[1] + along * u[1],
                    p[2] + along * u[2],
                ];
                let a_vec = sub(pos0, o);
                let b_vec = sub(q1, o);
                let tgt = sub(wb[j], o);
                let phi_best = (dot(tgt, b_vec)).atan2(dot(tgt, a_vec));
                let delta = wrap_pi(phi_best - phi_lat).clamp(-track_tweak, track_tweak);
                phi = phi_lat + delta;
            }
            vars[var] = wrap_pi(phi - offset);
        }
        tree.place_step(k, &vars, coords);
    }

    // Forward transform: the whole chain into walk space.
    for p in coords.iter_mut() {
        let rd = rot(sub(*p, c0), &r);
        *p = [w[0][0] + rd[0], w[0][1] + rd[1], w[0][2] + rd[2]];
    }
}
