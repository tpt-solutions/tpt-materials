//! Dilute (Maxwell) estimate for a non-interacting inclusion phase.
//!
//! For small `f_inclusion`, the effective stiffness is approximated by
//!
//! `C_eff = C_matrix + f_inclusion · (C_inclusion - C_matrix) : A`
//!
//! where `A` is the dilute strain-concentration tensor
//! (see [`tpt_mat_homogenization::dilute_strain_concentration`]).
//!
//! This estimate is accurate only for `f ≪ 1`; for moderate fractions,
//! use [`crate::mori_tanaka`].

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_math_linalg_fixed::Vec6;

use tpt_mat_homogenization::{dilute_strain_concentration, eshelby_spherical, EshelbySpherical};

/// Dilute effective stiffness for a two-phase composite.
///
/// `c_matrix` and `c_inclusion` are the matrix and inclusion
/// stiffnesses; `f_inclusion` is the inclusion volume fraction.
/// `eshelby` is the Eshelby tensor for an inclusion of the chosen
/// shape (typically spherical).
pub fn dilute_estimate(
    c_matrix: &SymmetricFourthOrder,
    c_inclusion: &SymmetricFourthOrder,
    f_inclusion: f64,
    eshelby: EshelbySpherical,
) -> SymmetricFourthOrder {
    let a = dilute_strain_concentration(c_matrix, c_inclusion, eshelby);
    // C_eff = C_m + f · (C_i - C_m) : A
    let mut dc = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            dc[i][j] = c_inclusion.data[i][j] - c_matrix.data[i][j];
        }
    }
    // D = (C_i - C_m) : A  (6x6 matrix product)
    let mut d = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut acc = 0.0;
            for k in 0..6 {
                acc += dc[i][k] * a.data[k][j];
            }
            d[i][j] = acc;
        }
    }
    // Effective stiffness.
    let mut c_eff = c_matrix.data;
    for i in 0..6 {
        for j in 0..6 {
            c_eff[i][j] += f_inclusion * d[i][j];
        }
    }
    SymmetricFourthOrder::new(c_eff)
}

/// Convenience constructor for spherical inclusions in an isotropic
/// matrix: derive `ν` from `E, ν` and build the Eshelby tensor.
pub fn eshelby_for_isotropic_sphere(matrix_e: f64, matrix_nu: f64) -> EshelbySpherical {
    eshelby_spherical(matrix_e, matrix_nu)
}

/// Apply the dilute estimate to a `Vec6` strain to get the average
/// stress response: `⟨σ⟩ = C_eff : ε`.  Useful as a sanity check.
pub fn dilute_response_stress(
    c_matrix: &SymmetricFourthOrder,
    c_inclusion: &SymmetricFourthOrder,
    f_inclusion: f64,
    eshelby: EshelbySpherical,
    eps: Vec6,
) -> Vec6 {
    let c_eff = dilute_estimate(c_matrix, c_inclusion, f_inclusion, eshelby);
    c_eff.contract(eps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn dilute_recovers_matrix_for_zero_inclusion() {
        let c0 = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c1 = SymmetricFourthOrder::isotropic(400_000.0, 0.25);
        let s = eshelby_for_isotropic_sphere(70_000.0, 0.33);
        let c = dilute_estimate(&c0, &c1, 0.0, s);
        for i in 0..6 {
            for j in 0..6 {
                assert_relative_eq!(c.data[i][j], c0.data[i][j], epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn dilute_response_stress_increases_with_stiff_inclusion() {
        let c0 = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c1 = SymmetricFourthOrder::isotropic(400_000.0, 0.25);
        let s = eshelby_for_isotropic_sphere(70_000.0, 0.33);
        let eps = Vec6::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let sigma = dilute_response_stress(&c0, &c1, 0.2, s, eps);
        let sigma_pure = c0.contract(eps);
        // Stiff inclusion ⇒ σ_xx must exceed the pure-matrix value.
        assert!(sigma[0] > sigma_pure[0]);
    }
}
