//! Symmetric 3x3 matrix in Voigt storage.
//!
//! Storage is `[xx, yy, zz, xy, yz, xz]`, giving direct access to the six
//! independent components of a symmetric tensor.  This is the storage
//! used by Cauchy stress, strain, and the symmetric part of `s ⊗ n`
//! (Schmid tensor).

use core::ops::{Add, Index, Mul, Sub};

use serde::{Deserialize, Serialize};

/// Symmetric 3x3 matrix stored in 6-component Voigt order.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct SymMat3 {
    /// `[xx, yy, zz, xy, yz, xz]` in Voigt order.
    pub data: [f64; 6],
}

impl SymMat3 {
    /// Construct from Voigt components `[xx, yy, zz, xy, yz, xz]`.
    #[inline]
    pub const fn new(xx: f64, yy: f64, zz: f64, xy: f64, yz: f64, xz: f64) -> Self {
        Self {
            data: [xx, yy, zz, xy, yz, xz],
        }
    }

    /// The zero matrix.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    /// Identity matrix.
    pub const IDENTITY: Self = Self::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

    /// Component access in `(i, j)` matrix notation, with `i, j ∈ {0, 1, 2}`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let (a, b) = (i.min(j), i.max(j));
        self.data[voigt(a, b)]
    }

    /// Map an index `k ∈ 0..6` back to `(i, j)` matrix indices.
    pub fn voigt_to_ij(k: usize) -> (usize, usize) {
        match k {
            0 => (0, 0),
            1 => (1, 1),
            2 => (2, 2),
            3 => (0, 1),
            4 => (1, 2),
            5 => (0, 2),
            _ => panic!("Voigt index out of range: {k}"),
        }
    }

    /// Scalar multiply (returns a new matrix).
    #[inline]
    pub fn scale(self, s: f64) -> Self {
        let d = self.data;
        Self::new(d[0] * s, d[1] * s, d[2] * s, d[3] * s, d[4] * s, d[5] * s)
    }

    /// Multiply every Voigt component by `s` in place (helper for builders).
    #[inline]
    pub fn mul_val(self, s: f64) -> Self {
        self.scale(s)
    }

    /// Frobenius norm squared: `‖A‖_F^2 = sum_ij A_ij^2`.
    #[inline]
    pub fn frobenius_sq(self) -> f64 {
        let d = self.data;
        d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + 2.0 * (d[3] * d[3] + d[4] * d[4] + d[5] * d[5])
    }

    /// Frobenius norm.
    #[inline]
    pub fn frobenius(self) -> f64 {
        self.frobenius_sq().sqrt()
    }

    /// Trace.
    #[inline]
    pub fn trace(self) -> f64 {
        self.data[0] + self.data[1] + self.data[2]
    }

    /// Double contraction `A : B = sum_ij A_ij B_ij`.
    pub fn double_dot(self, rhs: Self) -> f64 {
        let a = self.data;
        let b = rhs.data;
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])
    }
}

#[inline]
fn voigt(i: usize, j: usize) -> usize {
    match (i, j) {
        (0, 0) => 0,
        (1, 1) => 1,
        (2, 2) => 2,
        (0, 1) | (1, 0) => 3,
        (1, 2) | (2, 1) => 4,
        (0, 2) | (2, 0) => 5,
        _ => panic!("Voigt index out of range: ({i}, {j})"),
    }
}

impl Index<usize> for SymMat3 {
    type Output = f64;
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.data[i]
    }
}

impl Add for SymMat3 {
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

impl Sub for SymMat3 {
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

impl Mul<f64> for SymMat3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        self.scale(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn double_dot_equals_trace_of_product() {
        let a = SymMat3::new(1.0, 2.0, 3.0, 0.5, -0.2, 0.7);
        let b = SymMat3::IDENTITY;
        let dd = a.double_dot(b);
        assert_relative_eq!(dd, a.trace(), epsilon = 1e-12);
    }

    #[test]
    fn get_returns_symmetric_components() {
        let m = SymMat3::new(1.0, 2.0, 3.0, 0.5, -0.2, 0.7);
        assert_relative_eq!(m.get(0, 1), 0.5);
        assert_relative_eq!(m.get(1, 0), 0.5);
        assert_relative_eq!(m.get(2, 0), 0.7);
    }

    #[test]
    fn frobenius() {
        // For identity, ‖I‖_F = √3.
        assert_relative_eq!(SymMat3::IDENTITY.frobenius(), 3_f64.sqrt(), epsilon = 1e-12);
    }
}
