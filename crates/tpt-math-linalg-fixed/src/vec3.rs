//! 3D Euclidean vector with `f64` components.
//!
//! Used for lattice vectors, slip directions, slip-plane normals, and
//! grain centroids throughout `tpt-materials`.

use core::ops::{Add, AddAssign, Div, Index, Mul, Neg, Sub, SubAssign};

use serde::{Deserialize, Serialize};

use crate::SymMat3;

/// A column-style 3-vector.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Vec3 {
    /// Column-major storage: `[x, y, z]`.
    pub data: [f64; 3],
}

impl Vec3 {
    /// Construct from components.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { data: [x, y, z] }
    }

    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    /// The unit `x` vector.
    pub const X: Self = Self::new(1.0, 0.0, 0.0);

    /// The unit `y` vector.
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);

    /// The unit `z` vector.
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// Dot product `a · b`.
    #[inline]
    pub fn dot(self, rhs: Self) -> f64 {
        self.data[0] * rhs.data[0] + self.data[1] * rhs.data[1] + self.data[2] * rhs.data[2]
    }

    /// Cross product `a × b`.
    #[inline]
    pub fn cross(self, rhs: Self) -> Self {
        let [ax, ay, az] = self.data;
        let [bx, by, bz] = rhs.data;
        Self::new(ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx)
    }

    /// Squared L2 norm `a · a`.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    /// L2 norm `‖a‖`.
    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Returns `self / ‖self‖`.  Returns `ZERO` if `self` has zero norm.
    pub fn normalized(self) -> Self {
        let n = self.norm();
        if n == 0.0 {
            Self::ZERO
        } else {
            self * (1.0 / n)
        }
    }

    /// Symmetric outer product `sym(a ⊗ b) = ½(a ⊗ b + b ⊗ a)`.
    pub fn sym_outer(self, rhs: Self) -> SymMat3 {
        let [ax, ay, az] = self.data;
        let [bx, by, bz] = rhs.data;
        SymMat3::new(
            ax * bx,
            ay * by,
            az * bz,
            0.5 * (ax * by + ay * bx),
            0.5 * (ay * bz + az * by),
            0.5 * (ax * bz + az * bx),
        )
    }
}

impl From<[f64; 3]> for Vec3 {
    #[inline]
    fn from(data: [f64; 3]) -> Self {
        Self { data }
    }
}

impl From<Vec3> for [f64; 3] {
    #[inline]
    fn from(v: Vec3) -> Self {
        v.data
    }
}

impl Index<usize> for Vec3 {
    type Output = f64;

    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.data[i]
    }
}

impl Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(
            self.data[0] + rhs.data[0],
            self.data[1] + rhs.data[1],
            self.data[2] + rhs.data[2],
        )
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.data[0] += rhs.data[0];
        self.data[1] += rhs.data[1];
        self.data[2] += rhs.data[2];
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(
            self.data[0] - rhs.data[0],
            self.data[1] - rhs.data[1],
            self.data[2] - rhs.data[2],
        )
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.data[0] -= rhs.data[0];
        self.data[1] -= rhs.data[1];
        self.data[2] -= rhs.data[2];
    }
}

impl Neg for Vec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.data[0], -self.data[1], -self.data[2])
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.data[0] * rhs, self.data[1] * rhs, self.data[2] * rhs)
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f64) -> Self {
        Self::new(self.data[0] / rhs, self.data[1] / rhs, self.data[2] / rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn dot_and_cross() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert_relative_eq!(a.dot(b), 32.0);
        let c = a.cross(b);
        assert_relative_eq!(c.data[0], -3.0, max_relative = 1e-12);
        assert_relative_eq!(c.data[1], 6.0, max_relative = 1e-12);
        assert_relative_eq!(c.data[2], -3.0, max_relative = 1e-12);
        assert_relative_eq!(a.dot(c), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn normalise_zero_is_zero() {
        assert_eq!(Vec3::ZERO.normalized(), Vec3::ZERO);
        let n = Vec3::new(3.0, 0.0, 4.0).normalized();
        assert_relative_eq!(n.norm(), 1.0, epsilon = 1e-12);
    }
}
