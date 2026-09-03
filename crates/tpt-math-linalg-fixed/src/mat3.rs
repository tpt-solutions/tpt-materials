//! 3x3 matrix used for finite-deformation kinematics.
//!
//! Column-major storage so that `Mat3 * v` matches the right-multiplied
//! convention `v' = M v` and `Mat3 * Mat3` is matrix multiplication.

use core::ops::{Add, AddAssign, Index, Mul, Sub, SubAssign};

use serde::{Deserialize, Serialize};

use crate::Vec3;

/// 3x3 matrix in column-major storage.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Mat3 {
    /// 9 components in column-major order: `[m11, m21, m31, m12, m22, m32, m13, m23, m33]`.
    pub data: [f64; 9],
}

impl Mat3 {
    /// Construct from column-major components.
    #[inline]
    pub const fn new(
        m11: f64,
        m12: f64,
        m13: f64,
        m21: f64,
        m22: f64,
        m23: f64,
        m31: f64,
        m32: f64,
        m33: f64,
    ) -> Self {
        Self {
            data: [m11, m21, m31, m12, m22, m32, m13, m23, m33],
        }
    }

    /// Construct from row vectors.
    pub fn from_rows(r0: [f64; 3], r1: [f64; 3], r2: [f64; 3]) -> Self {
        Self::new(
            r0[0], r0[1], r0[2], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2],
        )
    }

    /// Identity matrix.
    pub const IDENTITY: Self = Self::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    /// Zero matrix.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    /// Component access `(i, j)`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[3 * j + i]
    }

    /// Set component `(i, j)`.
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[3 * j + i] = v;
    }

    /// Transpose.
    pub fn transpose(self) -> Self {
        let d = self.data;
        Self::new(d[0], d[3], d[6], d[1], d[4], d[7], d[2], d[5], d[8])
    }

    /// Trace.
    #[inline]
    pub fn trace(self) -> f64 {
        self.data[0] + self.data[4] + self.data[8]
    }

    /// Determinant.
    pub fn det(self) -> f64 {
        let d = self.data;
        d[0] * (d[4] * d[8] - d[5] * d[7]) - d[3] * (d[1] * d[8] - d[7] * d[2])
            + d[6] * (d[1] * d[5] - d[4] * d[2])
    }

    /// Inverse.  Panics on singular matrices.
    pub fn inverse(self) -> Self {
        let d = self.data;
        let det = self.det();
        if det.abs() < 1e-15 {
            panic!("Mat3 inverse: singular matrix");
        }
        let inv_det = 1.0 / det;
        // Column-major storage: data[3j+i] = m_ij.  Each slot of
        // `Mat3::new` corresponds to data[0]=m11, data[3]=m12,
        // data[6]=m13, data[1]=m21, data[4]=m22, data[7]=m23,
        // data[2]=m31, data[5]=m32, data[8]=m33.
        Self::new(
            (d[4] * d[8] - d[7] * d[5]) * inv_det,
            (d[6] * d[5] - d[3] * d[8]) * inv_det,
            (d[3] * d[7] - d[6] * d[4]) * inv_det,
            (d[7] * d[2] - d[1] * d[8]) * inv_det,
            (d[0] * d[8] - d[6] * d[2]) * inv_det,
            (d[6] * d[1] - d[0] * d[7]) * inv_det,
            (d[1] * d[5] - d[4] * d[2]) * inv_det,
            (d[3] * d[2] - d[0] * d[5]) * inv_det,
            (d[0] * d[4] - d[3] * d[1]) * inv_det,
        )
    }

    /// Symmetric part `½(M + Mᵀ)`.
    pub fn sym(self) -> Self {
        let t = self.transpose();
        // Column-major: data[3j+i] = m_ij.  t is computed via
        // transpose so t.data = [m11, m12, m13, m21, m22, m23, m31, m32, m33].
        let d = self.data;
        let m11 = (d[0] + t.data[0]) * 0.5;
        let m12 = (d[3] + t.data[1]) * 0.5;
        let m13 = (d[6] + t.data[2]) * 0.5;
        let m22 = (d[4] + t.data[4]) * 0.5;
        let m23 = (d[7] + t.data[5]) * 0.5;
        let m33 = (d[8] + t.data[8]) * 0.5;
        Self::new(m11, m12, m13, m12, m22, m23, m13, m23, m33)
    }

    /// Skew (anti-symmetric) part `½(M − Mᵀ)`.
    pub fn skew(self) -> Self {
        let t = self.transpose();
        let d = self.data;
        let m12 = (d[3] - t.data[1]) * 0.5;
        let m13 = (d[6] - t.data[2]) * 0.5;
        let m23 = (d[7] - t.data[5]) * 0.5;
        Self::new(0.0, m12, m13, -m12, 0.0, m23, -m13, -m23, 0.0)
    }

    /// Frobenius norm squared.
    #[inline]
    pub fn frobenius_sq(self) -> f64 {
        self.data.iter().map(|x| x * x).sum()
    }

    /// Frobenius norm.
    #[inline]
    pub fn frobenius(self) -> f64 {
        self.frobenius_sq().sqrt()
    }

    /// Scalar multiply.
    pub fn scale(self, s: f64) -> Self {
        let d = self.data;
        Self::new(
            d[0] * s,
            d[1] * s,
            d[2] * s,
            d[3] * s,
            d[4] * s,
            d[5] * s,
            d[6] * s,
            d[7] * s,
            d[8] * s,
        )
    }
}

impl Index<(usize, usize)> for Mat3 {
    type Output = f64;
    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        &self.data[3 * j + i]
    }
}

impl Add for Mat3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut out = [0.0_f64; 9];
        for i in 0..9 {
            out[i] = self.data[i] + rhs.data[i];
        }
        Self { data: out }
    }
}

impl AddAssign for Mat3 {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..9 {
            self.data[i] += rhs.data[i];
        }
    }
}

impl Sub for Mat3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut out = [0.0_f64; 9];
        for i in 0..9 {
            out[i] = self.data[i] - rhs.data[i];
        }
        Self { data: out }
    }
}

impl SubAssign for Mat3 {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..9 {
            self.data[i] -= rhs.data[i];
        }
    }
}

impl Mul<Mat3> for Mat3 {
    type Output = Self;
    fn mul(self, rhs: Mat3) -> Self {
        let a = self.data;
        let b = rhs.data;
        let mut out = [0.0_f64; 9];
        for j in 0..3 {
            for i in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += a[3 * k + i] * b[3 * j + k];
                }
                out[3 * j + i] = s;
            }
        }
        Self { data: out }
    }
}

impl Mul<Vec3> for Mat3 {
    type Output = Vec3;
    fn mul(self, v: Vec3) -> Vec3 {
        let d = self.data;
        Vec3::new(
            d[0] * v[0] + d[3] * v[1] + d[6] * v[2],
            d[1] * v[0] + d[4] * v[1] + d[7] * v[2],
            d[2] * v[0] + d[5] * v[1] + d[8] * v[2],
        )
    }
}

impl Mul<f64> for Mat3 {
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
    fn transpose_inverse() {
        let m = Mat3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0);
        let inv = m.inverse();
        let prod = m * inv;
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_relative_eq!(prod.get(i, j), expected, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn sym_skew_decompose() {
        let m = Mat3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        let s = m.sym();
        let k = m.skew();
        assert_relative_eq!((s + k).frobenius(), m.frobenius(), epsilon = 1e-10);
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(s.get(i, j), s.get(j, i), epsilon = 1e-12);
                assert_relative_eq!(k.get(i, j), -k.get(j, i), epsilon = 1e-12);
            }
        }
    }
}
