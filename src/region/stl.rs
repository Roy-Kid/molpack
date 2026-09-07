//! Closed triangle mesh as a [`Region`](super::Region), loaded from STL.
//!
//! Coordinates after construction are Å. `scale` on [`StlRegion::from_file`]
//! / [`StlRegion::from_bytes`] is Å per file unit (`1.0` means the file is
//! already Å). [`StlRegion::from_triangles`] takes Å vertices and has no
//! scale. Sign of `signed_distance` is molpack's: negative inside.
//!
//! Reading the file and welding its corners is molrs's job
//! ([`molrs::io::mesh::parse_stl`] into a [`TriMesh`], which also answers the
//! two questions this region gates on). What is left here is the part that is
//! actually molpack's: turning a closed surface into a signed distance.
//!
//! Both queries go through a [`Bvh`](super::bvh::Bvh) built once at
//! construction. It changes only which triangles are visited, so the answers
//! — and the goldens below — are the same ones the linear scan gave.

use std::path::Path;

use molrs::io::mesh::{parse_stl, read_stl};
use molrs::spatial::{TriMesh, mesh::DEGENERATE_AREA2};
use molrs::types::F;

use super::bvh::Bvh;
use super::{Aabb, CellDeclaration, Region};

const EPS: F = 1e-9;
const RAY: [F; 3] = [1.0, std::f64::consts::SQRT_2, std::f64::consts::PI];

/// Why an STL / triangle soup cannot become a [`StlRegion`].
#[derive(Debug, Clone)]
pub enum StlError {
    /// Filesystem failure while reading `path`.
    Io { path: String, message: String },
    /// Bytes are neither a length-matched binary STL nor ASCII `solid`.
    Parse { detail: String },
    /// No triangles after parse.
    Empty,
    /// `scale` was non-finite or `<= 0`.
    InvalidScale { scale: F },
    /// Triangle `index` has vanishing area after scaling.
    DegenerateTriangle { index: usize },
    /// Some directed half-edges have no opposite (not a closed 2-manifold).
    NotWatertight { unpaired: usize },
}

impl std::fmt::Display for StlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StlError::Io { path, message } => {
                write!(f, "could not read STL {path}: {message}")
            }
            StlError::Parse { detail } => write!(f, "STL parse failed: {detail}"),
            StlError::Empty => write!(f, "STL mesh has no triangles"),
            StlError::InvalidScale { scale } => {
                write!(f, "STL scale must be finite and > 0, got {scale}")
            }
            StlError::DegenerateTriangle { index } => {
                write!(f, "STL triangle {index} is degenerate (zero area)")
            }
            StlError::NotWatertight { unpaired } => write!(
                f,
                "STL mesh is not watertight: {unpaired} directed half-edge(s) lack an opposite"
            ),
        }
    }
}

impl std::error::Error for StlError {}

/// Watertight triangle mesh as a geometric [`Region`].
///
/// Interior is the even-odd filling of the triangle set. Distance magnitude
/// is the Euclidean closest point on any triangle (Eberly).
#[derive(Debug, Clone)]
pub struct StlRegion {
    triangles: Vec<[[F; 3]; 3]>,
    /// Spatial index over `triangles`, in the same indexing.
    bvh: Bvh,
    aabb: Aabb,
}

impl StlRegion {
    /// Å triangles. Same watertight / degeneracy gates as the STL loaders.
    pub fn from_triangles(triangles: &[[[F; 3]; 3]]) -> Result<Self, StlError> {
        Self::from_mesh(TriMesh::from_triangles(triangles))
    }

    /// Gate a welded mesh on what containment actually needs: something to be
    /// inside of, faces with a normal, and a closed surface. Without the last
    /// one the even-odd ray test calls points inside that are not, and the
    /// packing run confines molecules to a shape nobody drew.
    fn from_mesh(mesh: TriMesh) -> Result<Self, StlError> {
        if mesh.is_empty() {
            return Err(StlError::Empty);
        }
        if let Some(index) = mesh.first_degenerate_face(DEGENERATE_AREA2) {
            return Err(StlError::DegenerateTriangle { index });
        }
        let unpaired = mesh.unpaired_half_edges();
        if unpaired > 0 {
            return Err(StlError::NotWatertight { unpaired });
        }
        let (min, max) = mesh.aabb().ok_or(StlError::Empty)?;
        let triangles = mesh.to_triangles();
        let bvh = Bvh::build(&triangles);
        Ok(Self {
            triangles,
            bvh,
            aabb: Aabb { min, max },
        })
    }

    /// ASCII or binary STL. `scale` is Å per file unit; must be finite and `> 0`.
    pub fn from_bytes(bytes: &[u8], scale: F) -> Result<Self, StlError> {
        check_scale(scale)?;
        let mesh = parse_stl(bytes).map_err(|e| StlError::Parse {
            detail: e.to_string(),
        })?;
        Self::from_mesh(mesh.scaled(scale))
    }

    /// Read a path, then gate it like [`from_bytes`].
    ///
    /// Not gated on molpack's `io` feature: the reader is molrs's, and molpack
    /// takes molrs on its default feature union (which includes `io`), so this
    /// is available in every molpack build — the wheel included, which builds
    /// without molpack's own `io`.
    pub fn from_file(path: impl AsRef<Path>, scale: F) -> Result<Self, StlError> {
        let path = path.as_ref();
        check_scale(scale)?;
        // molrs reports both failures as `io::Error`; the kind is what says
        // whether the file was unreadable or its contents were not an STL.
        let mesh = read_stl(path).map_err(|e| match e.kind() {
            std::io::ErrorKind::InvalidData => StlError::Parse {
                detail: e.to_string(),
            },
            _ => StlError::Io {
                path: path.display().to_string(),
                message: e.to_string(),
            },
        })?;
        Self::from_mesh(mesh.scaled(scale))
    }

    fn unsigned_and_closest(&self, x: &[F; 3]) -> (F, [F; 3]) {
        let (d2, c) = self.bvh.nearest(x, |i| {
            let t = &self.triangles[i as usize];
            closest_point_triangle(*x, t[0], t[1], t[2])
        });
        (d2.sqrt(), c)
    }

    fn even_odd_inside(&self, x: &[F; 3]) -> bool {
        // `RAY` has no zero component by construction, so the slab test's
        // reciprocal is finite and the tree never has to special-case an axis.
        let rn = (RAY[0] * RAY[0] + RAY[1] * RAY[1] + RAY[2] * RAY[2]).sqrt();
        let dir = [RAY[0] / rn, RAY[1] / rn, RAY[2] / rn];
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        let hits = self.bvh.count_hits(x, &inv, |i| {
            let t = &self.triangles[i as usize];
            ray_hits_triangle(*x, dir, t[0], t[1], t[2])
        });
        hits % 2 == 1
    }
}

fn check_scale(scale: F) -> Result<(), StlError> {
    if !scale.is_finite() || scale <= 0.0 {
        Err(StlError::InvalidScale { scale })
    } else {
        Ok(())
    }
}

impl Region for StlRegion {
    fn contains(&self, x: &[F; 3]) -> bool {
        self.signed_distance(x) <= 0.0
    }

    fn signed_distance(&self, x: &[F; 3]) -> F {
        let (d, _) = self.unsigned_and_closest(x);
        if d < EPS {
            0.0
        } else if self.even_odd_inside(x) {
            -d
        } else {
            d
        }
    }

    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        let (d, c) = self.unsigned_and_closest(x);
        if d < EPS {
            return [0.0; 3];
        }
        let s = if self.even_odd_inside(x) { -1.0 } else { 1.0 };
        [
            s * (x[0] - c[0]) / d,
            s * (x[1] - c[1]) / d,
            s * (x[2] - c[2]) / d,
        ]
    }

    fn declared_cell(&self) -> Option<CellDeclaration> {
        None
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some(self.aabb)
    }

    fn is_closed_mesh(&self) -> bool {
        true
    }
}

fn sub(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn mul(a: [F; 3], s: F) -> [F; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: [F; 3], b: [F; 3]) -> F {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Ericson, Real-Time Collision Detection, §5.1.5.
fn closest_point_triangle(p: [F; 3], a: [F; 3], b: [F; 3], c: [F; 3]) -> [F; 3] {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add(a, mul(ab, v));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add(a, mul(ac, w));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, mul(sub(c, b), w));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    add(a, add(mul(ab, v), mul(ac, w)))
}

/// Möller–Trumbore, half-open barycentric domain.
fn ray_hits_triangle(orig: [F; 3], dir: [F; 3], v0: [F; 3], v1: [F; 3], v2: [F; 3]) -> bool {
    let e1 = sub(v1, v0);
    let e2 = sub(v2, v0);
    let pvec = cross(dir, e2);
    let det = dot(e1, pvec);
    if det.abs() < EPS {
        return false;
    }
    let inv = 1.0 / det;
    let tvec = sub(orig, v0);
    let u = dot(tvec, pvec) * inv;
    if !(0.0..1.0).contains(&u) {
        return false;
    }
    let qvec = cross(tvec, e1);
    let v = dot(dir, qvec) * inv;
    if v < 0.0 || u + v >= 1.0 {
        return false;
    }
    let t = dot(e2, qvec) * inv;
    t > EPS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::region::{RegionExt, RegionRestraint};
    use crate::restraint::AtomRestraint;

    fn cube_tris(lo: [F; 3], hi: [F; 3]) -> Vec<[[F; 3]; 3]> {
        let [x0, y0, z0] = lo;
        let [x1, y1, z1] = hi;
        let p = |x, y, z| [x, y, z];
        vec![
            // -x
            [p(x0, y0, z0), p(x0, y0, z1), p(x0, y1, z1)],
            [p(x0, y0, z0), p(x0, y1, z1), p(x0, y1, z0)],
            // +x
            [p(x1, y0, z0), p(x1, y1, z0), p(x1, y1, z1)],
            [p(x1, y0, z0), p(x1, y1, z1), p(x1, y0, z1)],
            // -y
            [p(x0, y0, z0), p(x1, y0, z0), p(x1, y0, z1)],
            [p(x0, y0, z0), p(x1, y0, z1), p(x0, y0, z1)],
            // +y
            [p(x0, y1, z0), p(x0, y1, z1), p(x1, y1, z1)],
            [p(x0, y1, z0), p(x1, y1, z1), p(x1, y1, z0)],
            // -z
            [p(x0, y0, z0), p(x0, y1, z0), p(x1, y1, z0)],
            [p(x0, y0, z0), p(x1, y1, z0), p(x1, y0, z0)],
            // +z
            [p(x0, y0, z1), p(x1, y0, z1), p(x1, y1, z1)],
            [p(x0, y0, z1), p(x1, y1, z1), p(x0, y1, z1)],
        ]
    }

    fn unit_cube() -> StlRegion {
        StlRegion::from_triangles(&cube_tris([0.0; 3], [1.0; 3])).expect("cube")
    }

    #[test]
    fn unit_cube_sdf_goldens() {
        let r = unit_cube();
        assert!(r.contains(&[0.5, 0.5, 0.5]));
        assert!((r.signed_distance(&[0.5, 0.5, 0.5]) + 0.5).abs() < 1e-12);
        let face = r.signed_distance(&[1.0, 0.5, 0.5]);
        assert!(r.contains(&[1.0, 0.5, 0.5]));
        assert!(face.abs() < 1e-9);
        assert!(!r.contains(&[2.0, 0.5, 0.5]));
        assert!((r.signed_distance(&[2.0, 0.5, 0.5]) - 1.0).abs() < 1e-12);
        assert!((r.signed_distance(&[0.5, 0.5, -1.0]) - 1.0).abs() < 1e-12);
        let corner = r.signed_distance(&[3.0, 3.0, 3.0]);
        assert!((corner - 2.0 * 3.0_f64.sqrt()).abs() < 1e-12);
        let g = r.signed_distance_grad(&[2.0, 0.5, 0.5]);
        assert!((g[0] - 1.0).abs() < 1e-9);
        assert!(g[1].abs() < 1e-9 && g[2].abs() < 1e-9);
        let bb = r.bounding_box().unwrap();
        assert_eq!(bb.min, [0.0; 3]);
        assert_eq!(bb.max, [1.0; 3]);
        assert!(r.declared_cell().is_none());
        assert!(r.is_closed_mesh());
        assert!(RegionRestraint(unit_cube()).is_closed_mesh());
        assert!(!crate::region::InsideBoxRegion::new([0.0; 3], [1.0; 3]).is_closed_mesh());
        assert!(
            unit_cube()
                .and(crate::region::InsideBoxRegion::new([0.0; 3], [1.0; 3]))
                .is_closed_mesh()
        );
    }

    #[test]
    fn nested_cavity_even_odd() {
        let mut tris = cube_tris([-2.0; 3], [2.0; 3]);
        tris.extend(cube_tris([-1.0; 3], [1.0; 3]));
        let r = StlRegion::from_triangles(&tris).expect("nested");
        assert!(!r.contains(&[0.0; 3]));
        assert!((r.signed_distance(&[0.0; 3]) - 1.0).abs() < 1e-9);
        assert!(r.contains(&[1.5, 0.0, 0.0]));
        assert!((r.signed_distance(&[1.5, 0.0, 0.0]) + 0.5).abs() < 1e-9);
        assert!(!r.contains(&[3.0, 0.0, 0.0]));
        assert!((r.signed_distance(&[3.0, 0.0, 0.0]) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn named_rejects() {
        assert!(matches!(
            StlRegion::from_triangles(&[]),
            Err(StlError::Empty)
        ));
        assert!(matches!(
            StlRegion::from_bytes(b"solid x\nendsolid x\n", 0.0),
            Err(StlError::InvalidScale { .. })
        ));
        let deg = [[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]];
        assert!(matches!(
            StlRegion::from_triangles(&deg),
            Err(StlError::DegenerateTriangle { index: 0 })
        ));
        let mut open = cube_tris([0.0; 3], [1.0; 3]);
        open.pop();
        assert!(matches!(
            StlRegion::from_triangles(&open),
            Err(StlError::NotWatertight { .. })
        ));
    }

    #[test]
    fn restraint_lift() {
        let r = unit_cube();
        let rest = RegionRestraint(r.clone());
        assert_eq!(rest.f(&[0.5, 0.5, 0.5], 1.0, 1.0), 0.0);
        assert!((rest.f(&[2.0, 0.5, 0.5], 1.0, 1.0) - 1.0).abs() < 1e-12);
        assert_eq!(r.into_restraint().f(&[0.5, 0.5, 0.5], 1.0, 1.0), 0.0);
    }

    #[test]
    fn scale_doubles_the_cube() {
        let unit = cube_tris([0.0; 3], [1.0; 3]);
        let ascii = ascii_cube(&unit);
        let r = StlRegion::from_bytes(ascii.as_bytes(), 2.0).expect("scale");
        assert!((r.signed_distance(&[1.0, 1.0, 1.0]) + 1.0).abs() < 1e-9);
    }

    #[test]
    fn ascii_and_binary_match_triangles() {
        let tris = cube_tris([0.0; 3], [1.0; 3]);
        let from_t = StlRegion::from_triangles(&tris).unwrap();
        let ascii = StlRegion::from_bytes(ascii_cube(&tris).as_bytes(), 1.0).unwrap();
        let bin = StlRegion::from_bytes(&binary_cube(&tris), 1.0).unwrap();
        let q = [0.5, 0.5, 0.5];
        assert!((from_t.signed_distance(&q) - ascii.signed_distance(&q)).abs() < 1e-9);
        assert!((from_t.signed_distance(&q) - bin.signed_distance(&q)).abs() < 1e-9);
        let solid_hdr = {
            let mut b = binary_cube(&tris);
            b[..5].copy_from_slice(b"solid");
            b
        };
        let via_solid = StlRegion::from_bytes(&solid_hdr, 1.0).unwrap();
        assert!((from_t.signed_distance(&q) - via_solid.signed_distance(&q)).abs() < 1e-9);
    }

    #[test]
    fn from_file_roundtrip() {
        let tris = cube_tris([0.0; 3], [1.0; 3]);
        let bytes = binary_cube(&tris);
        let dir = std::env::temp_dir();
        let path = dir.join("molpack-stl-region-cube.stl");
        std::fs::write(&path, &bytes).unwrap();
        let r = StlRegion::from_file(&path, 1.0).unwrap();
        assert!((r.signed_distance(&[0.5, 0.5, 0.5]) + 0.5).abs() < 1e-9);
        let _ = std::fs::remove_file(&path);
    }

    fn ascii_cube(tris: &[[[F; 3]; 3]]) -> String {
        let mut s = String::from("solid cube\n");
        for t in tris {
            s.push_str("  facet normal 0 0 0\n    outer loop\n");
            for v in t {
                s.push_str(&format!("      vertex {} {} {}\n", v[0], v[1], v[2]));
            }
            s.push_str("    endloop\n  endfacet\n");
        }
        s.push_str("endsolid cube\n");
        s
    }

    fn binary_cube(tris: &[[[F; 3]; 3]]) -> Vec<u8> {
        let mut b = vec![0u8; 80];
        b.extend_from_slice(&(tris.len() as u32).to_le_bytes());
        for t in tris {
            b.extend_from_slice(&[0u8; 12]);
            for v in t {
                for vk in v {
                    b.extend_from_slice(&(*vk as f32).to_le_bytes());
                }
            }
            b.extend_from_slice(&0u16.to_le_bytes());
        }
        b
    }
}
