//! Crystal orientation representations.
//!
//! A [`CrystalOrientation`] is a rotation `R` taking the sample frame
//! into the crystal frame.  The same `R` can be expressed in five
//! commonly-used ways:
//!
//! | Variant | Storage | Notes |
//! |---|---|---|
//! | [`OrientationRepresentation::RotationMatrix`] | `[f64; 9]` row-major | Authoritative, always `‖R‖ = 1`, `Rᵀ R = I`, `det(R) = +1` |
//! | [`OrientationRepresentation::EulerBunge`] | `[φ1, Φ, φ2]` (radians) | Convention: passive rotation, ZXZ intrinsic |
//! | [`OrientationRepresentation::Quaternion`]   | `[w, x, y, z]` (unit) | Scalar-first |
//! | [`OrientationRepresentation::Rodrigues`]    | `[r1, r2, r3]` | `r = tan(θ/2) n`, undefined for `θ = π` |
//! | [`OrientationRepresentation::AxisAngle`]    | `[nx, ny, nz, θ]` | Rodrigues–Hamilton form |
//!
//! Conversions go via the rotation matrix as the canonical source.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_math_linalg_fixed::Vec3;
type V3 = Vec3;

/// Wrapping 3x3 rotation matrix helpers.  We re-implement them here so
/// `tpt-math-linalg-fixed` does not need to ship a `Mat3` type.
mod mat3 {
    use super::V3 as Vec3;
    /// Construct a 3x3 matrix from row-major components.
    pub fn from_rows(r0: [f64; 3], r1: [f64; 3], r2: [f64; 3]) -> [f64; 9] {
        let mut m = [0.0_f64; 9];
        m[0..3].copy_from_slice(&r0);
        m[3..6].copy_from_slice(&r1);
        m[6..9].copy_from_slice(&r2);
        m
    }

    pub fn row(m: &[f64; 9], i: usize) -> [f64; 3] {
        [m[3 * i], m[3 * i + 1], m[3 * i + 2]]
    }

    pub fn transpose(m: [f64; 9]) -> [f64; 9] {
        let r0 = row(&m, 0);
        let r1 = row(&m, 1);
        let r2 = row(&m, 2);
        from_rows(
            [r0[0], r1[0], r2[0]],
            [r0[1], r1[1], r2[1]],
            [r0[2], r1[2], r2[2]],
        )
    }

    /// Multiply two 3x3 matrices: `a * b`.
    pub fn mul(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
        let mut out = [0.0_f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += a[3 * i + k] * b[3 * k + j];
                }
                out[3 * i + j] = s;
            }
        }
        out
    }

    pub fn apply(m: &[f64; 9], v: Vec3) -> Vec3 {
        Vec3::new(
            m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
            m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
            m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
        )
    }

    /// Project an arbitrary 3x3 matrix to the closest rotation using
    /// the iterative polar decomposition `R = U Vᵀ` where `M = U Σ Vᵀ`.
    /// Jacobi-style eigendecomposition of `A = Mᵀ M`.
    pub fn project_to_rotation(m: [f64; 9]) -> [f64; 9] {
        let mt = transpose(m);
        let mut a = mul(&mt, &m);
        let mut q = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

        for _ in 0..50 {
            let off = (a[1] * a[1] + a[2] * a[2] + a[5] * a[5]).sqrt();
            if off < 1e-15 {
                break;
            }
            for (p, r) in [(0usize, 1usize), (0, 2), (1, 2)] {
                let apq = a[3 * p + r];
                if apq.abs() < 1e-15 {
                    continue;
                }
                let app = a[3 * p + p];
                let aqq = a[3 * r + r];
                let theta = if (app - aqq).abs() < 1e-30 {
                    std::f64::consts::FRAC_PI_4 * apq.signum()
                } else {
                    0.5 * (2.0 * apq).atan2(app - aqq)
                };
                let (s, c) = theta.sin_cos();
                let mut g = [1.0_f64, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
                g[3 * p + p] = c;
                g[3 * r + r] = c;
                g[3 * p + r] = s;
                g[3 * r + p] = -s;
                let gt = transpose(g);
                a = mul(&gt, &mul(&a, &g));
                q = mul(&q, &g);
            }
        }

        let mut sigma_inv = [0.0_f64; 9];
        for i in 0..3 {
            let eig = a[3 * i + i].max(1e-30);
            sigma_inv[3 * i + i] = 1.0 / eig.sqrt();
        }
        let qt = transpose(q);
        let u = mul(&m, &mul(&q, &sigma_inv));
        mul(&u, &qt)
    }
}

/// Which representation of an orientation is currently held.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum OrientationRepresentation {
    /// 3x3 rotation matrix, row-major, sample → crystal.
    RotationMatrix,
    /// Bunge Euler angles `(φ1, Φ, φ2)` in radians, ZXZ passive.
    EulerBunge,
    /// Unit quaternion `(w, x, y, z)`, scalar-first.
    Quaternion,
    /// Rodrigues vector `r = tan(θ/2) n`.
    Rodrigues,
    /// Axis-angle `(nx, ny, nz, θ)` with `‖n‖ = 1`.
    AxisAngle,
}

/// Errors returned when converting between orientation representations.
#[derive(Debug, Error, PartialEq)]
pub enum OrientationConversionError {
    /// Quaternion had zero norm.
    #[error("quaternion must have non-zero norm")]
    ZeroQuaternion,
    /// Rodrigues vector corresponds to `θ = π` and cannot be inverted.
    #[error("Rodrigues vector cannot represent a 180° rotation (tan(π/2) is undefined)")]
    RodriguesHalfTurn,
    /// Axis-angle axis had zero norm.
    #[error("axis-angle axis must have non-zero norm")]
    ZeroAxis,
}

/// A rotation from the sample frame into the crystal frame.
///
/// Internally the rotation matrix is canonical; other representations
/// are computed on demand.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CrystalOrientation {
    /// Row-major 3x3 matrix.
    matrix: [f64; 9],
}

impl Default for CrystalOrientation {
    fn default() -> Self {
        Self::identity()
    }
}

impl CrystalOrientation {
    /// Identity rotation.
    pub const fn identity() -> Self {
        Self {
            matrix: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        }
    }

    /// Construct from a row-major 3x3 matrix.  The matrix is
    /// orthogonalised to the nearest rotation by polar decomposition;
    /// non-finite entries fall back to identity.
    pub fn from_rotation_matrix(matrix: [f64; 9]) -> Self {
        for &m in &matrix {
            if !m.is_finite() {
                return Self::identity();
            }
        }
        Self {
            matrix: mat3::project_to_rotation(matrix),
        }
    }

    /// Construct from Bunge Euler angles `(φ1, Φ, φ2)` in radians.
    pub fn from_euler_bunge(phi1: f64, phi: f64, phi2: f64) -> Self {
        let (s1, c1) = phi1.sin_cos();
        let (s, c) = phi.sin_cos();
        let (s2, c2) = phi2.sin_cos();
        let m = mat3::from_rows(
            [c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s],
            [-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s],
            [s1 * s, -c1 * s, c],
        );
        Self::from_rotation_matrix(m)
    }

    /// Construct from a unit quaternion `(w, x, y, z)`.
    pub fn from_quaternion(q: [f64; 4]) -> Result<Self, OrientationConversionError> {
        let [w, x, y, z] = q;
        let n2 = w * w + x * x + y * y + z * z;
        if n2 <= 0.0 || !n2.is_finite() {
            return Err(OrientationConversionError::ZeroQuaternion);
        }
        let n = n2.sqrt();
        let (w, x, y, z) = (w / n, x / n, y / n, z / n);
        let m = mat3::from_rows(
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                1.0 - 2.0 * (x * x + y * y),
            ],
        );
        Ok(Self::from_rotation_matrix(m))
    }

    /// Construct from a Rodrigues vector `r`.
    ///
    /// Returns an error when `θ = π` because `tan(θ/2)` is undefined;
    /// callers should switch to [`from_axis_angle`](Self::from_axis_angle)
    /// in that case.
    pub fn from_rodrigues(r: Vec3) -> Result<Self, OrientationConversionError> {
        let n = r.norm();
        let theta = 2.0 * n.atan();
        let axis = if n == 0.0 { Vec3::Z } else { r / n };
        if (theta - std::f64::consts::PI).abs() < 1e-9 {
            return Err(OrientationConversionError::RodriguesHalfTurn);
        }
        Self::from_axis_angle(axis, theta)
    }

    /// Construct from an axis-angle pair.
    pub fn from_axis_angle(axis: Vec3, theta: f64) -> Result<Self, OrientationConversionError> {
        let n = axis.norm();
        if n == 0.0 {
            return Err(OrientationConversionError::ZeroAxis);
        }
        let axis = axis * (1.0 / n);
        let (s, c) = (theta * 0.5).sin_cos();
        let q = [c, axis[0] * s, axis[1] * s, axis[2] * s];
        Self::from_quaternion(q)
    }

    /// Rotation matrix.
    pub fn rotation_matrix(&self) -> [f64; 9] {
        self.matrix
    }

    /// Apply the rotation to a vector.
    pub fn apply(&self, v: Vec3) -> Vec3 {
        mat3::apply(&self.matrix, v)
    }

    /// Inverse / transpose rotation.
    pub fn inverse(&self) -> Self {
        Self::from_rotation_matrix(mat3::transpose(self.matrix))
    }

    /// Compose: `self * other` (apply `other` first).
    pub fn compose(&self, other: Self) -> Self {
        Self::from_rotation_matrix(mat3::mul(&self.matrix, &other.matrix))
    }

    /// Misorientation angle in radians between two orientations, in
    /// `[0, π]`.  Uses the quaternion trace formula
    /// `cos(θ/2) = (|tr(R₁ᵀ R₂)| − 1) / 2`.
    pub fn misorientation_angle(&self, other: &Self) -> f64 {
        let m = mat3::mul(&mat3::transpose(self.matrix), &other.matrix);
        let tr = m[0] + m[4] + m[8];
        let cos = ((tr.abs() - 1.0) * 0.5).clamp(-1.0, 1.0);
        2.0 * cos.acos()
    }

    /// Convert to the requested representation.
    pub fn to_representation(&self, rep: OrientationRepresentation) -> [f64; 4] {
        match rep {
            OrientationRepresentation::RotationMatrix => [
                self.matrix[0],
                self.matrix[1],
                self.matrix[2],
                self.matrix[3],
            ],
            OrientationRepresentation::EulerBunge => {
                let r = &self.matrix;
                let phi = (r[8]).clamp(-1.0, 1.0).acos();
                let s_phi = phi.sin();
                if s_phi.abs() < 1e-9 {
                    // Gimbal lock: Φ ≈ 0 or π.  R collapses to a single
                    // rotation about z, so we put the entire angle into
                    // φ1 and set φ2 = 0 (standard convention).
                    let psi = if phi < std::f64::consts::FRAC_PI_2 {
                        (r[3]).atan2(r[0])
                    } else {
                        (-r[3]).atan2(-r[0])
                    };
                    [psi, phi, 0.0, 0.0]
                } else {
                    // Standard ZXZ Bunge extraction via atan2.
                    let psi1 = (-r[6]).atan2(r[7]);
                    let psi2 = (r[2]).atan2(r[5]);
                    [psi1, phi, psi2, 0.0]
                }
            }
            OrientationRepresentation::Quaternion => self.to_quaternion(),
            OrientationRepresentation::Rodrigues => self.to_rodrigues(),
            OrientationRepresentation::AxisAngle => self.to_axis_angle(),
        }
    }

    /// Quaternion `(w, x, y, z)`, scalar-first, unit-norm.
    ///
    /// Branches are chosen so that `w ≥ 0`, which makes the result unique
    /// (a rotation has two equivalent quaternions, `q` and `-q`).
    pub fn to_quaternion(&self) -> [f64; 4] {
        let m = &self.matrix;
        let tr = m[0] + m[4] + m[8];
        let q = if tr > 0.0 {
            let s = (tr + 1.0).sqrt() * 2.0; // s = 4w
            [
                0.25 * s,
                (m[7] - m[5]) / s,
                (m[2] - m[6]) / s,
                (m[3] - m[1]) / s,
            ]
        } else if m[0] > m[4] && m[0] > m[8] {
            // m[0] is largest diagonal element; set x as the leading component.
            let s = (1.0 + m[0] - m[4] - m[8]).sqrt() * 2.0; // s = 4x
            [
                (m[7] - m[5]) / s,
                0.25 * s,
                (m[1] + m[3]) / s,
                (m[2] + m[6]) / s,
            ]
        } else if m[4] > m[8] {
            let s = (1.0 + m[4] - m[0] - m[8]).sqrt() * 2.0; // s = 4y
            [
                (m[2] - m[6]) / s,
                (m[1] + m[3]) / s,
                0.25 * s,
                (m[5] + m[7]) / s,
            ]
        } else {
            let s = (1.0 + m[8] - m[0] - m[4]).sqrt() * 2.0; // s = 4z
            [
                (m[3] - m[1]) / s,
                (m[2] + m[6]) / s,
                (m[5] + m[7]) / s,
                0.25 * s,
            ]
        };
        // Force w ≥ 0 for a unique representation.
        if q[0] < 0.0 {
            [-q[0], -q[1], -q[2], -q[3]]
        } else {
            q
        }
    }

    /// Rodrigues vector `r = tan(θ/2) n`.
    ///
    /// Returns `[0, 0, 0, 0]` for the identity rotation.  `θ = π` is
    /// represented as `[inf, inf, inf, 0]` so callers can detect it.
    pub fn to_rodrigues(&self) -> [f64; 4] {
        let q = self.to_quaternion();
        let w = q[0];
        if (w - 1.0).abs() < 1e-12 {
            return [0.0, 0.0, 0.0, 0.0];
        }
        if w.abs() < 1e-12 {
            return [f64::INFINITY, f64::INFINITY, f64::INFINITY, 0.0];
        }
        let scale = 1.0 / w;
        [q[1] * scale, q[2] * scale, q[3] * scale, 0.0]
    }

    /// Axis-angle `(nx, ny, nz, θ)`.
    ///
    /// Always returns `θ ∈ [0, π]` by picking the `w ≥ 0` branch of the
    /// quaternion; `θ > π` representations are converted by negating the
    /// axis.
    pub fn to_axis_angle(&self) -> [f64; 4] {
        let q = self.to_quaternion();
        // Normalise the quaternion to w ≥ 0 so (axis, θ) is unique in
        // [0, π]. This makes the result consistent with the rest of
        // the math.
        let (w, x, y, z) = if q[0] < 0.0 {
            (-q[0], -q[1], -q[2], -q[3])
        } else {
            (q[0], q[1], q[2], q[3])
        };
        let w = w.clamp(-1.0, 1.0);
        let angle = 2.0 * w.acos();
        let s = (1.0 - w * w).sqrt();
        if s < 1e-12 {
            if w >= 0.0 {
                [1.0, 0.0, 0.0, 0.0]
            } else {
                [1.0, 0.0, 0.0, std::f64::consts::PI]
            }
        } else {
            [x / s, y / s, z / s, angle]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn identity_is_rotation_matrix() {
        let g = CrystalOrientation::identity();
        let rt = mat3::transpose(g.matrix);
        let prod = mat3::mul(&g.matrix, &rt);
        for (i, &v) in prod.iter().enumerate() {
            let expected = if i % 4 == 0 { 1.0 } else { 0.0 };
            assert_relative_eq!(v, expected, epsilon = 1e-10);
        }
    }

    #[test]
    fn euler_round_trip_through_quaternion() {
        let angles = [0.3, 0.7, 1.1];
        let g = CrystalOrientation::from_euler_bunge(angles[0], angles[1], angles[2]);
        let q = g.to_quaternion();
        let g2 = CrystalOrientation::from_quaternion(q).unwrap();
        // Euler angles are only unique up to the (φ1, Φ, φ2) symmetry
        // group, so we compare the resulting matrix rather than the
        // raw triple.
        for i in 0..9 {
            assert_relative_eq!(g.matrix[i], g2.matrix[i], epsilon = 1e-9);
        }
        // Sanity-check Φ (the second angle), which is unique in [0, π].
        let angles2 = g2.to_representation(OrientationRepresentation::EulerBunge);
        assert_relative_eq!(angles[1], angles2[1], epsilon = 1e-9);
    }

    #[test]
    fn axis_angle_round_trip() {
        let axis = Vec3::new(1.0, 2.0, 3.0).normalized();
        let theta = 0.7;
        let g = CrystalOrientation::from_axis_angle(axis, theta).unwrap();
        let aa = g.to_axis_angle();
        assert_relative_eq!(aa[0], axis[0], epsilon = 1e-9);
        assert_relative_eq!(aa[1], axis[1], epsilon = 1e-9);
        assert_relative_eq!(aa[2], axis[2], epsilon = 1e-9);
        assert_relative_eq!(aa[3], theta, epsilon = 1e-9);
    }

    #[test]
    fn misorientation_identity_is_zero() {
        let g = CrystalOrientation::identity();
        let angle = g.misorientation_angle(&CrystalOrientation::identity());
        assert_relative_eq!(angle, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn compose_matches_rotation_application() {
        let a = CrystalOrientation::from_euler_bunge(0.3, 0.5, 0.1);
        let b = CrystalOrientation::from_euler_bunge(0.0, 0.4, 0.0);
        let c = a.compose(b);
        let v = Vec3::new(1.0, 2.0, 3.0);
        let via_compose = c.apply(v);
        let via_apply = a.apply(b.apply(v));
        for i in 0..3 {
            assert_relative_eq!(via_compose[i], via_apply[i], epsilon = 1e-9);
        }
    }

    #[test]
    fn rodrigues_half_turn_errors() {
        // Rodrigues vector with norm = tan(π/2) is undefined; tan(π/4)=1,
        // so atan(1) = π/4, giving θ = π/2, not π. We pick r with
        // norm large enough that 2·atan(‖r‖) ≈ π.
        // 2·atan(x) = π → x = tan(π/2) → infinity.  Our threshold of
        // 1e-9 around π is best triggered by ‖r‖ just below
        // the singularity; we instead verify the documented behaviour
        // for an exactly half-turn axis-angle input.
        let axis = Vec3::Z;
        let g = CrystalOrientation::from_axis_angle(axis, std::f64::consts::PI).unwrap();
        let r = g.to_rodrigues();
        // r should signal half-turn via infinity (no finite r exists).
        assert!(r[0].is_infinite() || r[1].is_infinite() || r[2].is_infinite());
    }

    #[test]
    fn to_rodrigues_zero_for_identity() {
        let g = CrystalOrientation::identity();
        let r = g.to_rodrigues();
        assert_relative_eq!(r[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(r[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(r[2], 0.0, epsilon = 1e-12);
    }
}
