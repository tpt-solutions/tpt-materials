//! Anisotropic elastic stiffness tensor `C` in 6x6 Voigt form.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_math_linalg_fixed::Vec6;

/// 6x6 stiffness matrix in Voigt order (Hooke's law `σ = C : ε`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SymmetricFourthOrder {
    /// Row-major 6x6 stiffness values.
    pub data: [[f64; 6]; 6],
}

impl SymmetricFourthOrder {
    /// Construct a 6x6 stiffness from row data.  The matrix must be
    /// symmetric (`C[i][j] == C[j][i]`) — `new` enforces this and
    /// symmetrises by averaging if violated.
    pub fn new(data: [[f64; 6]; 6]) -> Self {
        let mut c = data;
        for i in 0..6 {
            for j in (i + 1)..6 {
                let avg = 0.5 * (c[i][j] + c[j][i]);
                c[i][j] = avg;
                c[j][i] = avg;
            }
        }
        Self { data: c }
    }

    /// Cubic symmetry (`C11, C12, C44`) — appropriate for FCC, BCC and
    /// diamond-cubic crystals.  All values in MPa.
    pub fn cubic(c11: f64, c12: f64, c44: f64) -> Self {
        let mut c = [[0.0_f64; 6]; 6];
        c[0][0] = c11;
        c[1][1] = c11;
        c[2][2] = c11;
        c[0][1] = c12;
        c[1][0] = c12;
        c[0][2] = c12;
        c[2][0] = c12;
        c[1][2] = c12;
        c[2][1] = c12;
        c[3][3] = c44;
        c[4][4] = c44;
        c[5][5] = c44;
        Self::new(c)
    }

    /// Hexagonal symmetry (`C11, C33, C12, C13, C44`).
    pub fn hexagonal(c11: f64, c33: f64, c12: f64, c13: f64, c44: f64) -> Self {
        let mut c = [[0.0_f64; 6]; 6];
        c[0][0] = c11;
        c[1][1] = c11;
        c[2][2] = c33;
        c[0][1] = c12;
        c[1][0] = c12;
        c[0][2] = c13;
        c[2][0] = c13;
        c[1][2] = c13;
        c[2][1] = c13;
        c[3][3] = c44;
        c[4][4] = c44;
        c[5][5] = (c11 - c12) * 0.5;
        Self::new(c)
    }

    /// Isotropic stiffness from Young's modulus `E` and Poisson's ratio `ν`.
    pub fn isotropic(e: f64, nu: f64) -> Self {
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let c11 = lambda + 2.0 * mu;
        let c12 = lambda;
        let c44 = mu;
        Self::cubic(c11, c12, c44)
    }

    /// Apply `σ = C : ε` for two Voigt vectors.
    pub fn contract(&self, eps: Vec6) -> Vec6 {
        let mut out = [0.0_f64; 6];
        for i in 0..6 {
            let mut s = 0.0;
            for j in 0..6 {
                s += self.data[i][j] * eps[j];
            }
            out[i] = s;
        }
        Vec6::new(out[0], out[1], out[2], out[3], out[4], out[5])
    }

    /// Compliance `S = C^{-1}` (6x6 dense inverse).
    pub fn compliance(&self) -> [[f64; 6]; 6] {
        let m = self.data;
        inv6(m)
    }
}

#[derive(Debug, Error, PartialEq)]
/// Errors when constructing an [`ElasticStiffness`].
pub enum ElasticStiffnessError {
    /// Stiffness matrix was not positive-definite.
    #[error("elastic stiffness is not positive-definite")]
    NotPositiveDefinite,
}

/// Alias retained for backwards compatibility / clearer docs.
pub type ElasticStiffness = SymmetricFourthOrder;

/// Invert a symmetric 6x6 matrix via Gaussian elimination with partial
/// pivoting.  Only used to recover the compliance from a known stiffness
/// (tests).
fn inv6(m: [[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut a = [[0.0_f64; 12]; 6];
    for i in 0..6 {
        for j in 0..6 {
            a[i][j] = m[i][j];
            a[i][j + 6] = if i == j { 1.0 } else { 0.0 };
        }
    }
    for i in 0..6 {
        // pivot
        let mut piv = i;
        for k in (i + 1)..6 {
            if a[k][i].abs() > a[piv][i].abs() {
                piv = k;
            }
        }
        if (a[piv][i]).abs() < 1e-15 {
            panic!("singular stiffness matrix");
        }
        a.swap(i, piv);
        let inv_piv = 1.0 / a[i][i];
        for j in 0..12 {
            a[i][j] *= inv_piv;
        }
        for k in 0..6 {
            if k != i {
                let f = a[k][i];
                for j in 0..12 {
                    a[k][j] -= f * a[i][j];
                }
            }
        }
    }
    let mut out = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = a[i][j + 6];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use tpt_math_linalg_fixed::SymMat3;

    #[test]
    fn cubic_is_symmetric() {
        let c = SymmetricFourthOrder::cubic(200.0, 100.0, 50.0);
        for i in 0..6 {
            for j in 0..6 {
                assert_eq!(c.data[i][j], c.data[j][i]);
            }
        }
    }

    #[test]
    fn isotropic_round_trip() {
        let c = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let eps = Vec6::new(0.001, -0.0003, -0.0003, 0.0, 0.0, 0.0);
        let sigma = c.contract(eps);
        // σ_xx should equal E · ε_xx for uniaxial strain.
        assert!((sigma[0] - 200.0).abs() < 1e-6, "got {}", sigma[0]);
    }

    #[test]
    fn cubic_satisfies_c11_minus_c12_positive() {
        // Mechanical-stability criterion (Born condition for cubic).
        let c = SymmetricFourthOrder::cubic(200.0, 100.0, 50.0);
        assert!(c.data[0][0] - c.data[0][1] > 0.0);
        assert!(c.data[3][3] > 0.0);
    }

    #[test]
    fn compliance_recovers_identity_under_contract() {
        let c = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let s = c.compliance();
        let mut cs = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let mut acc = 0.0;
                for k in 0..6 {
                    acc += c.data[i][k] * s[k][j];
                }
                cs[i][j] = acc;
            }
        }
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_relative_eq!(cs[i][j], expected, epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn symmetric_average_when_input_asymmetric() {
        let mut m = [[0.0_f64; 6]; 6];
        m[0][1] = 10.0;
        m[1][0] = 6.0;
        let c = SymmetricFourthOrder::new(m);
        assert_eq!(c.data[0][1], 8.0);
        assert_eq!(c.data[1][0], 8.0);
    }

    #[test]
    fn contract_sym_mat3_matches_voigt() {
        // For uniaxial strain along x: σ_xx = c11 × ε_xx where
        // c11 = E(1-ν)/((1+ν)(1-2ν)).
        let c = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let eps = Vec6::new(0.01, 0.0, 0.0, 0.0, 0.0, 0.0);
        let sigma_voigt = c.contract(eps);
        let lambda = 70_000.0 * 0.33 / (1.33 * 0.34);
        let mu = 70_000.0 / (2.0 * 1.33);
        let c11 = lambda + 2.0 * mu;
        let expected = c11 * 0.01;
        assert!(
            (sigma_voigt[0] - expected).abs() < 1.0,
            "got {}, expected {}",
            sigma_voigt[0],
            expected
        );
    }
}
