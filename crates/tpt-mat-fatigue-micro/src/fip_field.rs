//! Compute the fatigue-indicator parameter (FIP) field from a
//! converged [`CpFemResult`].
//!
//! Given the per-grain (or per-integration-point) state — stress
//! tensor, accumulated shear per slip system — this driver returns a
//! scalar FIP per grain under the chosen
//! [`FatigueCriterion`](super::FatigueCriterion).

use tpt_mat_crystal_plasticity::CpFemResult;
use tpt_math_linalg_fixed::{Mat3, Vec6};

use super::{fatemi_socie_fip, findley_fip, swt_fip, FatigueCriterion};

/// Compute the FIP field.
///
/// Returns a `Vec<f64>` of length `result.accumulated_shear.len()`
/// (one value per grain / integration point).
pub fn fatigue_indicator_parameter(result: &CpFemResult, criterion: FatigueCriterion) -> Vec<f64> {
    let n = result.accumulated_shear.len();
    let mut fip = Vec::with_capacity(n);
    for i in 0..n {
        let v = match criterion {
            FatigueCriterion::Findley { k } => {
                let s = result.stresses[i];
                let sym = s.to_sym_mat3_data();
                let s_max = principal_max(&sym);
                let s_min = principal_min(&sym);
                findley_fip(s_max, s_min, k)
            }
            FatigueCriterion::FatemiSocie { sigma_y, k } => {
                let s = result.stresses[i];
                let sym = s.to_sym_mat3_data();
                let s_max = principal_max(&sym);
                let strain = result.strains[i];
                let eps_half = 0.5 * strain_norm(&strain.to_sym_mat3_data());
                fatemi_socie_fip(eps_half, s_max, sigma_y, k)
            }
            FatigueCriterion::SmithWatsonTopper => {
                let s = result.stresses[i];
                let sym = s.to_sym_mat3_data();
                let s_max = principal_max(&sym);
                let s_min = principal_min(&sym);
                let strain = result.strains[i];
                let eps_total = 2.0 * strain_norm(&strain.to_sym_mat3_data());
                swt_fip(s_max, s_min, eps_total)
            }
            FatigueCriterion::CrystallographicSlip {
                critical_accumulated_shear,
            } => {
                // FIP = max_slip(γ_acc / γ_c), saturated at 1.0.
                let acc = &result.accumulated_shear[i];
                let max_acc = acc.iter().copied().fold(0.0_f64, f64::max);
                (max_acc / critical_accumulated_shear)
                    .min(f64::INFINITY)
                    .max(0.0)
            }
        };
        fip.push(v);
    }
    fip
}

/// Was initiation reached on this grain at this FIP?
pub fn initiation_reached(criterion: FatigueCriterion, fip: f64) -> bool {
    match criterion {
        FatigueCriterion::CrystallographicSlip { .. } => fip >= 1.0,
        // Other criteria require Coffin–Manson calibration; we use
        // an FIP ≥ 1 as a proxy for "unit-damage per cycle".
        _ => fip >= 1.0,
    }
}

fn principal_max(s: &[[f64; 3]; 3]) -> f64 {
    let eig = eigenvalues_symmetric(s);
    eig.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn principal_min(s: &[[f64; 3]; 3]) -> f64 {
    let eig = eigenvalues_symmetric(s);
    eig.iter().copied().fold(f64::INFINITY, f64::min)
}

fn strain_norm(s: &[[f64; 3]; 3]) -> f64 {
    let mut acc = 0.0_f64;
    for i in 0..3 {
        acc += s[i][i] * s[i][i];
    }
    for i in 0..3 {
        for j in (i + 1)..3 {
            acc += 2.0 * s[i][j] * s[i][j];
        }
    }
    acc.sqrt()
}

/// Approximate eigenvalues of a symmetric 3×3 matrix by the iterative
/// Jacobi method (3×3 only, so cheap).
fn eigenvalues_symmetric(a: &[[f64; 3]; 3]) -> [f64; 3] {
    let mut m = *a;
    for _ in 0..40 {
        let mut off = 0.0_f64;
        for i in 0..3 {
            for j in (i + 1)..3 {
                off += m[i][j] * m[i][j];
            }
        }
        if off < 1.0e-24 {
            break;
        }
        for p in 0..3 {
            for q in (p + 1)..3 {
                if m[p][q].abs() < 1.0e-18 {
                    continue;
                }
                let theta = if (m[p][p] - m[q][q]).abs() < 1.0e-24 {
                    core::f64::consts::FRAC_PI_4
                } else {
                    0.5 * (((2.0 * m[p][q]) / (m[p][p] - m[q][q])).atan())
                };
                let c = theta.cos();
                let s = theta.sin();
                let app = c * c * m[p][p] + 2.0 * c * s * m[p][q] + s * s * m[q][q];
                let aqq = s * s * m[p][p] - 2.0 * c * s * m[p][q] + c * c * m[q][q];
                let apq = (c * c - s * s) * m[p][q] + c * s * (m[q][q] - m[p][p]);
                m[p][p] = app;
                m[q][q] = aqq;
                m[p][q] = apq;
                m[q][p] = apq;
                for i in 0..3 {
                    if i != p && i != q {
                        let aip = c * m[i][p] + s * m[i][q];
                        let aiq = -s * m[i][p] + c * m[i][q];
                        m[i][p] = aip;
                        m[i][q] = aiq;
                        m[p][i] = aip;
                        m[q][i] = aiq;
                    }
                }
            }
        }
    }
    [m[0][0], m[1][1], m[2][2]]
}

// `Vec6` helpers local to this module.
trait Vec6Ext {
    fn to_sym_mat3_data(&self) -> [[f64; 3]; 3];
}
impl Vec6Ext for Vec6 {
    fn to_sym_mat3_data(&self) -> [[f64; 3]; 3] {
        let d = self.data;
        [[d[0], d[3], d[5]], [d[3], d[1], d[4]], [d[5], d[4], d[2]]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_result(n_grains: usize, gamma: f64) -> CpFemResult {
        let slip = vec![0.0_f64; 12];
        CpFemResult {
            stresses: vec![Vec6::new(100.0, 0.0, 0.0, 0.0, 0.0, 0.0); n_grains],
            strains: vec![Vec6::new(0.001, 0.0, 0.0, 0.0, 0.0, 0.0); n_grains],
            slip_rates: vec![slip; n_grains],
            lattice_rotations: vec![Mat3::IDENTITY; n_grains],
            accumulated_shear: vec![vec![gamma; 12]; n_grains],
            reaction_force: Vec::new(),
        }
    }

    #[test]
    fn findley_field_constant_loading() {
        let crit = FatigueCriterion::Findley { k: 0.3 };
        let res = dummy_result(5, 0.0);
        let fip = fatigue_indicator_parameter(&res, crit);
        assert_eq!(fip.len(), 5);
        for v in &fip {
            assert!(*v > 0.0);
        }
    }

    #[test]
    fn crystallographic_slip_critical_field() {
        let crit = FatigueCriterion::CrystallographicSlip {
            critical_accumulated_shear: 0.15,
        };
        let res = dummy_result(3, 0.20);
        let fip = fatigue_indicator_parameter(&res, crit);
        for v in &fip {
            assert!(*v > 1.0);
        }
    }

    #[test]
    fn initiation_reached_threshold() {
        assert!(initiation_reached(
            FatigueCriterion::Findley { k: 0.3 },
            1.5
        ));
        assert!(!initiation_reached(
            FatigueCriterion::Findley { k: 0.3 },
            0.5
        ));
        assert!(initiation_reached(
            FatigueCriterion::CrystallographicSlip {
                critical_accumulated_shear: 0.10
            },
            1.0
        ));
    }
}
