//! 6-component stress / strain vector in Voigt notation.
//!
//! Storage: `[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]`.  The Cauchy stress and
//! the small-strain tensor share this representation; the conversion to /
//! from [`SymMat3`](crate::SymMat3) is lossless.

use core::ops::{Add, Index, Mul, Sub};

use serde::{Deserialize, Serialize};

use crate::SymMat3;

/// 6-component Voigt vector (symmetric stress or strain).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Vec6 {
    /// `[xx, yy, zz, xy, yz, xz]`.
    pub data: [f64; 6],
}

impl Vec6 {
    /// Construct from Voigt components.
    #[inline]
    pub const fn new(xx: f64, yy: f64, zz: f64, xy: f64, yz: f64, xz: f64) -> Self {
        Self {
            data: [xx, yy, zz, xy, yz, xz],
        }
    }

    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    /// Convert to a [`SymMat3`].
    pub fn to_sym_mat3(self) -> SymMat3 {
        let d = self.data;
        SymMat3::new(d[0], d[1], d[2], d[3], d[4], d[5])
    }

    /// Construct from a [`SymMat3`].
    pub fn from_sym_mat3(m: SymMat3) -> Self {
        Self { data: m.data }
    }

    /// Scalar multiply.
    pub fn scale(self, s: f64) -> Self {
        let d = self.data;
        Self::new(d[0] * s, d[1] * s, d[2] * s, d[3] * s, d[4] * s, d[5] * s)
    }

    /// Dot product `a · b` with the engineering-shear convention
    /// (shear components weighted by 2 so the result matches `A : B`).
    pub fn double_dot(self, rhs: Self) -> f64 {
        let a = self.data;
        let b = rhs.data;
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])
    }
}

impl Index<usize> for Vec6 {
    type Output = f64;
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.data[i]
    }
}

impl Add for Vec6 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let a = self.data;
        let b = rhs.data;
        Self::new(
            a[0] + b[0],
            a[1] + b[1],
            a[2] + b[2],
            a[3] + b[3],
            a[4] + b[4],
            a[5] + b[5],
        )
    }
}

impl Sub for Vec6 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let a = self.data;
        let b = rhs.data;
        Self::new(
            a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2],
            a[3] - b[3],
            a[4] - b[4],
            a[5] - b[5],
        )
    }
}

impl Mul<f64> for Vec6 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        self.scale(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn round_trip_sym_mat3() {
        let v = Vec6::new(1.0, 2.0, 3.0, 0.5, -0.2, 0.7);
        let m = v.to_sym_mat3();
        let v2 = Vec6::from_sym_mat3(m);
        assert_eq!(v, v2);
    }

    #[test]
    fn double_dot_matches_sym_mat3() {
        let a = Vec6::new(1.0, 2.0, 3.0, 0.5, -0.2, 0.7);
        let b = Vec6::new(4.0, 5.0, 6.0, 0.1, 0.2, 0.3);
        assert_relative_eq!(
            a.double_dot(b),
            a.to_sym_mat3().double_dot(b.to_sym_mat3()),
            epsilon = 1e-12
        );
    }
}
