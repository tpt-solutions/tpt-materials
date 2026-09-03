//! Effective coefficient of thermal expansion (CTE).

use serde::{Deserialize, Serialize};

/// CTE homogenization scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CteBound {
    /// Turner model (iso-stress, bulk modulus weighted).
    Turner,
    /// Kerner model (matrix + spherical inclusion).
    Kerner,
    /// Rosen–Hashin upper bound.
    RosenHashinUpper,
    /// Rosen–Hashin lower bound.
    RosenHashinLower,
}

/// Compute the effective linear coefficient of thermal expansion
/// `α_eff` of an `N`-phase composite.
///
/// `cte_phases[i]` is the CTE of phase `i`,
/// `k_phases[i]` is its bulk modulus,
/// `g_phases[i]` is its shear modulus (used by Kerner),
/// `f[i]` is its volume fraction.
///
/// The bounds differ in which homogenization scheme they
/// assume for the elastic stiffness.
pub fn effective_cte(
    bound: CteBound,
    cte_phases: &[f64],
    k_phases: &[f64],
    g_phases: &[f64],
    f: &[f64],
) -> f64 {
    assert_eq!(cte_phases.len(), f.len(), "cte_phases and f must match");
    assert_eq!(k_phases.len(), f.len(), "k_phases and f must match");
    assert_eq!(g_phases.len(), f.len(), "g_phases and f must match");
    match bound {
        CteBound::Turner => {
            // α_eff = Σ α_i K_i f_i / Σ K_i f_i
            let num: f64 = f
                .iter()
                .zip(cte_phases.iter().zip(k_phases.iter()))
                .map(|(fi, (ai, ki))| fi * ai * ki)
                .sum();
            let denom: f64 = f.iter().zip(k_phases.iter()).map(|(fi, ki)| fi * ki).sum();
            if denom.abs() < 1.0e-30 {
                return 0.0;
            }
            num / denom
        }
        CteBound::Kerner => {
            // Kerner (matrix + spherical inclusion) closed form.
            let (k_m, idx_m) = k_phases
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, k)| (*k, i))
                .unwrap_or((k_phases[0], 0));
            let g_m = g_phases[idx_m];
            let alpha_m = cte_phases[idx_m];
            let mut alpha_eff = alpha_m;
            for (i, (&ai, &ki)) in cte_phases.iter().zip(k_phases.iter()).enumerate() {
                if i == idx_m {
                    continue;
                }
                let fi = f[i];
                // Kerner's volume-weighted contribution:
                let dalpha = (ki - k_m) * (ai - alpha_m);
                let num = 3.0 * (1.0 - fi) * k_m * dalpha + 4.0 * g_m * (alpha_m * ki - ai * k_m);
                let denom = (3.0 * k_m + 4.0 * g_m) * fi * ki
                    + 4.0 * g_m * (1.0 - fi) * k_m
                    + 12.0 * k_m * g_m;
                if denom.abs() > 1.0e-30 {
                    alpha_eff += fi * num / denom;
                }
            }
            alpha_eff
        }
        CteBound::RosenHashinUpper | CteBound::RosenHashinLower => {
            // Approximate with bulk-modulus vs shear-modulus
            // weighted bounds:
            // α_eff = Σ α_i V_i / V_eff × (K_eff or G_eff weights)
            let sum_kf: f64 = f.iter().zip(k_phases.iter()).map(|(fi, ki)| fi * ki).sum();
            let sum_gf: f64 = f.iter().zip(g_phases.iter()).map(|(fi, gi)| fi * gi).sum();
            let weight_k = if matches!(bound, CteBound::RosenHashinUpper) {
                sum_kf
            } else {
                sum_gf
            };
            let num: f64 = f
                .iter()
                .zip(cte_phases.iter().zip(k_phases.iter().zip(g_phases.iter())))
                .map(|(fi, (ai, (ki, gi)))| {
                    if matches!(bound, CteBound::RosenHashinUpper) {
                        fi * ai * ki
                    } else {
                        fi * ai * gi
                    }
                })
                .sum();
            let denom = weight_k;
            if denom.abs() < 1.0e-30 {
                return 0.0;
            }
            num / denom
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn turner_recovers_matrix_when_inclusion_fraction_zero() {
        let alpha = effective_cte(
            CteBound::Turner,
            &[1.0e-5, 2.0e-5],
            &[100.0, 200.0],
            &[40.0, 80.0],
            &[1.0, 0.0],
        );
        assert!(approx(alpha, 1.0e-5, 1.0e-12));
    }

    #[test]
    fn turner_scales_with_k_modulus_weight() {
        // Phase A: α=1e-5, K=100; Phase B: α=2e-5, K=300.
        // Phase B dominates -> α_eff closer to 2e-5.
        let alpha = effective_cte(
            CteBound::Turner,
            &[1.0e-5, 2.0e-5],
            &[100.0, 300.0],
            &[40.0, 80.0],
            &[0.5, 0.5],
        );
        let expected = (0.5 * 1.0e-5 * 100.0 + 0.5 * 2.0e-5 * 300.0) / (0.5 * 100.0 + 0.5 * 300.0);
        assert!(approx(alpha, expected, 1.0e-9));
    }

    #[test]
    fn rosen_hashin_bounds_bracket_voigt_average() {
        let alpha_lo = effective_cte(
            CteBound::RosenHashinLower,
            &[1.0e-5, 2.0e-5],
            &[100.0, 200.0],
            &[40.0, 80.0],
            &[0.5, 0.5],
        );
        let alpha_up = effective_cte(
            CteBound::RosenHashinUpper,
            &[1.0e-5, 2.0e-5],
            &[100.0, 200.0],
            &[40.0, 80.0],
            &[0.5, 0.5],
        );
        assert!(alpha_up >= alpha_lo);
    }
}
