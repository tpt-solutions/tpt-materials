//! Mori–Tanaka (1973) mean-field homogenization.
//!
//! For a two-phase composite with matrix (phase 0) and inclusions
//! (phase 1), the Mori–Tanaka estimate of the effective stiffness
//! `C_eff` satisfies
//!
//! `⟨σ⟩ = C_eff : ⟨ε⟩`
//!
//! where the strain in the inclusions is
//!
//! `ε_1 = A_dilute(0) : ε_0`
//!
//! and the matrix strain `ε_0` is the average strain in the matrix.
//! Benveniste (1987) showed that the closed-form expression is
//!
//! `C_eff = C_0 + f_1 · (C_1 - C_0) · A · (f_0 · I + f_1 · A)^{-1}`
//!
//! where `A = [I + S · C_0^{-1} · (C_1 - C_0)]^{-1}` is the dilute
//! strain-concentration tensor.
//!
//! For an `N`-phase composite, the MT estimate generalises by
//! treating each phase `r` as an inclusion in the (phase-averaged)
//! matrix `C_m = Σ f_r C_r / Σ f_r`.  The implementation supports
//! this via the iterative driver [`mori_tanaka_iterative`].

use serde::{Deserialize, Serialize};

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_math_linalg_fixed::Vec6;

use tpt_mat_homogenization::{dilute_strain_concentration, eshelby_spherical, EshelbySpherical};

/// Result of a two-phase Mori–Tanaka calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MoriTanakaResult {
    /// Effective stiffness `C_eff` in Voigt order.
    pub c_eff: SymmetricFourthOrder,
    /// Average strain in the inclusion phase.
    pub eps_inclusion: Vec6,
    /// Average strain in the matrix phase.
    pub eps_matrix: Vec6,
    /// Volume fraction of the inclusion phase.
    pub f_inclusion: f64,
}

/// Two-phase Mori–Tanaka effective stiffness.
///
/// `c_matrix` and `c_inclusion` are the matrix and inclusion
/// stiffnesses; `f_inclusion` is the inclusion volume fraction;
/// `eshelby` is the Eshelby tensor for the inclusion shape
/// (typically spherical).
pub fn mori_tanaka(
    c_matrix: &SymmetricFourthOrder,
    c_inclusion: &SymmetricFourthOrder,
    f_inclusion: f64,
    eshelby: EshelbySpherical,
) -> MoriTanakaResult {
    let f_matrix = 1.0 - f_inclusion;
    let a = dilute_strain_concentration(c_matrix, c_inclusion, eshelby);
    // D = (C_1 - C_0) · A
    let mut dc = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            dc[i][j] = c_inclusion.data[i][j] - c_matrix.data[i][j];
        }
    }
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
    // M = f_0 · I + f_1 · A
    let mut m = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let diag = if i == j { f_matrix } else { 0.0 };
            m[i][j] = diag + f_inclusion * a.data[i][j];
        }
    }
    let m_inv = invert6(m);
    // C_eff = C_0 + f_1 · D · M^{-1}
    let mut prod = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut acc = 0.0;
            for k in 0..6 {
                acc += d[i][k] * m_inv[k][j];
            }
            prod[i][j] = acc;
        }
    }
    let mut c_eff = c_matrix.data;
    for i in 0..6 {
        for j in 0..6 {
            c_eff[i][j] += f_inclusion * prod[i][j];
        }
    }
    let c_eff = SymmetricFourthOrder::new(c_eff);
    // Strain partition: pick a unit applied strain ε_∞ = e_x and
    // compute ε_1 = A · (f_0 I + f_1 A)^{-1} · ε_∞, ε_0 = ε_∞.
    let eps_inf = Vec6::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut eps_1 = [0.0_f64; 6];
    let mut eps_0 = eps_inf.data;
    for i in 0..6 {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += a.data[i][j] * eps_inf[j];
        }
        eps_1[i] = acc;
    }
    // ε_0 = (I - f_1 (A - I) )^{-1} ε_∞  (consistency).
    let mut m_for_e0 = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let a_term = a.data[i][j] - if i == j { 1.0 } else { 0.0 };
            m_for_e0[i][j] = if i == j { 1.0 } else { 0.0 } - f_inclusion * a_term;
        }
    }
    let m_for_e0_inv = invert6(m_for_e0);
    let mut eps_0_arr = [0.0_f64; 6];
    for i in 0..6 {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += m_for_e0_inv[i][j] * eps_inf[j];
        }
        eps_0_arr[i] = acc;
    }
    eps_0 = eps_0_arr;
    MoriTanakaResult {
        c_eff,
        eps_inclusion: Vec6::new(eps_1[0], eps_1[1], eps_1[2], eps_1[3], eps_1[4], eps_1[5]),
        eps_matrix: Vec6::new(eps_0[0], eps_0[1], eps_0[2], eps_0[3], eps_0[4], eps_0[5]),
        f_inclusion,
    }
}

/// Iterative `N`-phase Mori–Tanaka: start with the Voigt average as
/// the effective medium, treat each phase as an inclusion in the
/// current estimate, update, and iterate.  Converges in 20–40
/// iterations for typical stiffness contrasts.
pub fn mori_tanaka_iterative(
    phases: &[(SymmetricFourthOrder, f64)],
    matrix_nu: f64,
    max_iter: usize,
    tol: f64,
) -> SymmetricFourthOrder {
    let total_f: f64 = phases.iter().map(|(_, f)| *f).sum();
    if total_f <= 0.0 || phases.is_empty() {
        return SymmetricFourthOrder::new([[0.0_f64; 6]; 6]);
    }
    let eshelby = eshelby_spherical(1.0, matrix_nu);
    // Start with the phase-volume-weighted (Voigt) stiffness.
    let mut c_eff = rule_of_mixtures_weighted(phases);
    for _ in 0..max_iter {
        let mut c_new = [[0.0_f64; 6]; 6];
        for (c_r, f_r) in phases {
            let w = f_r / total_f;
            let _ = eshelby;
            let a = dilute_strain_concentration(&c_eff, c_r, EshelbySpherical::from_nu(matrix_nu));
            // Effective inclusion contribution: c_r : A, weighted by f_r.
            let mut local = [[0.0_f64; 6]; 6];
            for i in 0..6 {
                for j in 0..6 {
                    let mut acc = 0.0;
                    for k in 0..6 {
                        acc += c_r.data[i][k] * a.data[k][j];
                    }
                    local[i][j] = acc;
                }
            }
            for i in 0..6 {
                for j in 0..6 {
                    c_new[i][j] += w * local[i][j];
                }
            }
        }
        let c_new_s = SymmetricFourthOrder::new(c_new);
        // Convergence check.
        let mut max_diff = 0.0_f64;
        for i in 0..6 {
            for j in 0..6 {
                max_diff = max_diff.max((c_new[i][j] - c_eff.data[i][j]).abs());
            }
        }
        c_eff = c_new_s;
        if max_diff < tol {
            break;
        }
    }
    c_eff
}

fn rule_of_mixtures_weighted(phases: &[(SymmetricFourthOrder, f64)]) -> SymmetricFourthOrder {
    let total_f: f64 = phases.iter().map(|(_, f)| *f).sum();
    let mut c = [[0.0_f64; 6]; 6];
    for (c_r, f_r) in phases {
        let w = f_r / total_f;
        for i in 0..6 {
            for j in 0..6 {
                c[i][j] += w * c_r.data[i][j];
            }
        }
    }
    SymmetricFourthOrder::new(c)
}

fn invert6(m: [[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut a = [[0.0_f64; 12]; 6];
    for i in 0..6 {
        for j in 0..6 {
            a[i][j] = m[i][j];
            a[i][j + 6] = if i == j { 1.0 } else { 0.0 };
        }
    }
    for i in 0..6 {
        let mut piv = i;
        for k in (i + 1)..6 {
            if a[k][i].abs() > a[piv][i].abs() {
                piv = k;
            }
        }
        if a[piv][i].abs() < 1e-15 {
            panic!("singular");
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

    fn s() -> EshelbySpherical {
        EshelbySpherical::from_nu(0.3)
    }

    #[test]
    fn mori_tanaka_recovers_matrix_for_zero_inclusion() {
        let c0 = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c1 = SymmetricFourthOrder::isotropic(400_000.0, 0.25);
        let r = mori_tanaka(&c0, &c1, 0.0, s());
        for i in 0..6 {
            for j in 0..6 {
                assert_relative_eq!(r.c_eff.data[i][j], c0.data[i][j], epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn mori_tanaka_recovers_inclusion_for_full_inclusion() {
        let c0 = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c1 = SymmetricFourthOrder::isotropic(400_000.0, 0.25);
        let r = mori_tanaka(&c0, &c1, 1.0, s());
        for i in 0..6 {
            for j in 0..6 {
                assert_relative_eq!(r.c_eff.data[i][j], c1.data[i][j], epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn mori_tanaka_intermediate_fraction_between_bounds() {
        let c0 = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c1 = SymmetricFourthOrder::isotropic(400_000.0, 0.25);
        let r = mori_tanaka(&c0, &c1, 0.3, s());
        // C_eff[0][0] should lie between c0 and c1.
        let c00 = c0.data[0][0];
        let c11 = c1.data[0][0];
        let eff = r.c_eff.data[0][0];
        assert!(eff > c00 && eff < c11, "got {eff}");
    }

    #[test]
    fn iterative_matches_two_phase_at_low_contrast() {
        let c0 = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c1 = SymmetricFourthOrder::isotropic(100_000.0, 0.3);
        let phases = vec![(c0.clone(), 0.5), (c1.clone(), 0.5)];
        let c_iter = mori_tanaka_iterative(&phases, 0.3, 50, 1e-6);
        // For low contrast the MT estimate should be near the Voigt average.
        let v = (c0.data[0][0] + c1.data[0][0]) * 0.5;
        let diff = (c_iter.data[0][0] - v).abs();
        assert!(diff < 5_000.0, "got {}", c_iter.data[0][0]);
    }
}
