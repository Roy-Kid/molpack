//! The one geometric restraint: stay inside a molrs region.
//!
//! Geometry is molrs's ([`molrs::core`]: shapes, meshes, unions
//! of spheres, and their `And` / `Or` / `Not` compositions). What is
//! molpack's is the penalty that turns a region's signed distance into a
//! term of the shared objective, and that is all this type is.

use std::sync::Arc;

use molrs::core::Region;
use molrs::op::F;

use super::AtomRestraint;

/// Lift any molrs [`Region`] to a quadratic exterior penalty:
///
/// ```text
/// penalty(x) = scale · max(0, distance(x))²
/// ```
///
/// Zero everywhere inside the region and on its boundary; outside it grows
/// with the square of the region's signed distance, so its gradient is
/// `2 · scale · distance · ∇distance`, the region's outward direction.
/// The penalty is quadratic in the length by which the boundary is missed,
/// so it consumes `scale` like the `.inp` box and plane kernels (the
/// two-scale contract on [`AtomRestraint`](crate::AtomRestraint)): `precision = 0.01` reads as "within
/// 0.1 Å" of any region, and a `Cuboid` reproduces the box kernel's value.
///
/// The region is shared: `Arc<dyn Region + Send + Sync>` is what molrs's
/// own combinators take and what a region handed across the Python
/// boundary resolves to, so there is one lift and no generic twin.
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
/// use molpack::RegionRestraint;
/// use molpack::AtomRestraint;
/// use molrs::core::Sphere;
/// use ndarray::array;
///
/// let ball = RegionRestraint(Arc::new(Sphere::new(array![0.0, 0.0, 0.0], 5.0)));
/// assert_eq!(ball.f(&[1.0, 0.0, 0.0], 1.0, 1.0), 0.0);
/// assert_eq!(ball.f(&[10.0, 0.0, 0.0], 1.0, 1.0), 25.0);
/// ```
#[derive(Clone, Debug)]
pub struct RegionRestraint(pub Arc<dyn Region + Send + Sync>);

impl RegionRestraint {
    /// Is the region's bounding box finite on every axis `shift` moves along?
    ///
    /// If it is, the atoms the region allows never leave one image along that
    /// shift, so the restraint reads the same there whatever the cell does.
    fn bounded_along(&self, shift: [F; 3]) -> bool {
        let b = self.0.bounds();
        (0..3).all(|k| shift[k] == 0.0 || (b[[k, 0]].is_finite() && b[[k, 1]].is_finite()))
    }
}

impl AtomRestraint for RegionRestraint {
    fn f(&self, x: &[F; 3], scale: F, _scale2: F) -> F {
        let v = self.0.distance(x).max(0.0);
        scale * v * v
    }

    fn fg(&self, x: &[F; 3], scale: F, _scale2: F, g: &mut [F; 3]) -> F {
        let d = self.0.distance(x);
        if d <= 0.0 {
            return 0.0;
        }
        let grad = self.0.distance_grad(x);
        let coeff = 2.0 * scale * d;
        for k in 0..3 {
            g[k] += coeff * grad[k];
        }
        scale * d * d
    }

    fn name(&self) -> &'static str {
        "RegionRestraint"
    }

    /// Both halves of the rule come from the region itself: it repeats along
    /// `shift` (molrs [`Region::repeats_along`]), or it is bounded along it.
    ///
    /// [`Region::repeats_along`]: molrs::core::Region::repeats_along
    fn holds_along(&self, shift: [F; 3]) -> bool {
        self.0.repeats_along(shift) || self.bounded_along(shift)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molrs::core::{AndRegion, Cuboid, NotRegion, Sphere};
    use ndarray::array;

    fn ball() -> RegionRestraint {
        RegionRestraint(Arc::new(Sphere::new(array![0.0, 0.0, 0.0], 5.0)))
    }

    #[test]
    fn penalty_zero_inside_and_on_the_boundary() {
        let r = ball();
        assert_eq!(r.f(&[1.0, 2.0, 0.0], 1.0, 1.0), 0.0);
        assert_eq!(r.f(&[5.0, 0.0, 0.0], 1.0, 1.0), 0.0);
        let mut g = [0.0; 3];
        assert_eq!(r.fg(&[1.0, 2.0, 0.0], 1.0, 1.0, &mut g), 0.0);
        assert_eq!(g, [0.0; 3]);
    }

    /// The Euclidean distance of the sphere makes the penalty exactly
    /// `scale · (d − r)²`; `scale2` is not consumed.
    #[test]
    fn penalty_positive_outside_is_quadratic_in_the_distance() {
        let r = ball();
        assert_eq!(r.f(&[10.0, 0.0, 0.0], 1.0, 1.0), 25.0);
        assert_eq!(r.f(&[10.0, 0.0, 0.0], 0.5, 1.0), 12.5);
        assert_eq!(r.f(&[10.0, 0.0, 0.0], 1.0, 0.01), 25.0);
    }

    #[test]
    fn gradient_matches_finite_difference_on_a_composition() {
        let region: Arc<dyn Region + Send + Sync> = Arc::new(AndRegion::new(
            Arc::new(Cuboid::new(array![0.0, 0.0, 0.0], array![10.0, 10.0, 10.0])),
            Arc::new(NotRegion::new(Arc::new(Sphere::new(
                array![5.0, 5.0, 5.0],
                2.0,
            )))),
        ));
        let r = RegionRestraint(region);
        for x in [[15.0, 5.0, 5.0], [5.0, 5.5, 5.0], [-2.0, 3.0, 4.0]] {
            let mut g = [0.0; 3];
            let v = r.fg(&x, 0.7, 1.0, &mut g);
            assert!((v - r.f(&x, 0.7, 1.0)).abs() < 1e-12);
            for k in 0..3 {
                let h = 1e-6;
                let mut xp = x;
                xp[k] += h;
                let mut xm = x;
                xm[k] -= h;
                let fd = (r.f(&xp, 0.7, 1.0) - r.f(&xm, 0.7, 1.0)) / (2.0 * h);
                assert!(
                    (g[k] - fd).abs() < 1e-6,
                    "axis {k}: {} vs {fd} at {x:?}",
                    g[k]
                );
            }
        }
    }

    #[test]
    fn gradient_accumulates_not_overwrite() {
        let r = ball();
        let mut g = [1.0, 2.0, 3.0];
        r.fg(&[10.0, 0.0, 0.0], 1.0, 1.0, &mut g);
        assert!((g[0] - 11.0).abs() < 1e-12);
        assert_eq!(g[1], 2.0);
        assert_eq!(g[2], 3.0);
    }

    #[test]
    fn contains_point_iff_penalty_is_zero() {
        let r = ball();
        for x in [
            [0.0; 3],
            [4.9, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [5.1, 0.0, 0.0],
            [9.0, 9.0, 9.0],
        ] {
            assert_eq!(r.0.contains_point(&x), r.f(&x, 1.0, 1.0) == 0.0, "at {x:?}");
        }
    }

    #[test]
    fn declares_nothing_and_is_parallel_safe() {
        let r = ball();
        assert!(r.is_parallel_safe());
        assert!(r.declared_cell().is_none());
        assert_eq!(r.name(), "RegionRestraint");
    }
}
