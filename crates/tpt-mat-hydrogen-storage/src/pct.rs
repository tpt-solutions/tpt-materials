//! Pressure–composition–temperature (PCT) isotherms.

use serde::{Deserialize, Serialize};

const R_GAS: f64 = 8.314_462;

/// PCT isotherm parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PctParams {
    /// Plateau pressure at the reference temperature
    /// `P_0(T_ref)` (Pa).
    pub p0_at_tref: f64,
    /// Reference temperature (K).
    pub t_ref_k: f64,
    /// Enthalpy of desorption `ΔH` (J/mol H₂).
    pub delta_h_j_per_mol: f64,
    /// Entropy of desorption `ΔS` (J/(mol H₂·K)).
    pub delta_s_j_per_mol_k: f64,
    /// Maximum hydrogen-to-metal ratio `H/M_max`.
    pub h_to_m_max: f64,
    /// Slope of the lower plateau branch `m_low` (Pa per H/M).
    pub slope_low: f64,
    /// Slope of the upper plateau branch `m_up` (Pa per H/M).
    pub slope_up: f64,
}

/// One point on a PCT isotherm.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PctPoint {
    /// Hydrogen-to-metal ratio `H/M` (`0 ≤ H/M ≤ H/M_max`).
    pub h_to_m: f64,
    /// Equilibrium pressure (Pa).
    pub pressure_pa: f64,
}

/// Van 't Hoff equilibrium pressure
/// `ln(P / P_ref) = -ΔH / R (1/T − 1/T_ref)`
/// at temperature `t_k`, where `P_ref = P_0(T_ref)` is the
/// plateau pressure at the reference temperature.
pub fn van_t_hoff_pressure(p: &PctParams, t_k: f64) -> f64 {
    let exponent = -p.delta_h_j_per_mol / R_GAS * (1.0 / t_k - 1.0 / p.t_ref_k);
    (p.p0_at_tref * exponent.exp()).max(1.0)
}

/// Compute the PCT isotherm as a `Vec<PctPoint>` at temperature
/// `t_k` with `n_points` samples uniformly spaced in
/// `H/M ∈ [0, H/M_max]`.
pub fn pct_isotherm(p: &PctParams, t_k: f64, n_points: usize) -> Vec<PctPoint> {
    if n_points < 2 {
        return Vec::new();
    }
    let p0 = van_t_hoff_pressure(p, t_k);
    let mut pts = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let h_to_m = p.h_to_m_max * (i as f64) / ((n_points - 1) as f64);
        // Two-branch plateau:
        // lower branch (H/M < 0.5 H/M_max): P = P0 + m_low * (H/M - 0.5 H/M_max)
        // upper branch (H/M > 0.5 H/M_max): P = P0 + m_up * (H/M - 0.5 H/M_max)
        let mid = 0.5 * p.h_to_m_max;
        let slope = if h_to_m < mid { p.slope_low } else { p.slope_up };
        let pressure = (p0 + slope * (h_to_m - mid)).max(1.0);
        pts.push(PctPoint { h_to_m, pressure_pa: pressure });
    }
    pts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn lani5_params() -> PctParams {
        // LaNi5: ΔH_des ≈ +30 kJ/mol (endothermic desorption,
        // so increasing T raises the plateau pressure).
        PctParams {
            p0_at_tref: 2.0e5, // ~2 bar at 298 K
            t_ref_k: 298.0,
            delta_h_j_per_mol: 30_000.0,
            delta_s_j_per_mol_k: 100.0,
            h_to_m_max: 6.0,
            slope_low: 1.0e5,
            slope_up: 1.0e5,
        }
    }

    #[test]
    fn pct_returns_correct_length() {
        let p = lani5_params();
        let iso = pct_isotherm(&p, 298.0, 21);
        assert_eq!(iso.len(), 21);
    }

    #[test]
    fn van_t_hoff_recovers_reference_pressure() {
        let p = lani5_params();
        let p_eq = van_t_hoff_pressure(&p, p.t_ref_k);
        assert!(approx(p_eq, p.p0_at_tref, 1.0e-6));
    }

    #[test]
    fn pressure_at_plateau_midpoint_equals_van_t_hoff() {
        let p = lani5_params();
        let p0 = van_t_hoff_pressure(&p, 350.0);
        let iso = pct_isotherm(&p, 350.0, 13);
        let mid = iso[6]; // H/M = 3.0 = 0.5 * 6
        assert!(approx(mid.pressure_pa, p0, 1.0e-3));
    }

    #[test]
    fn pressure_decreases_with_temperature_at_fixed_h_to_m() {
        // For exothermic desorption, increasing T raises P_eq.
        let p = lani5_params();
        let iso_lo = pct_isotherm(&p, 300.0, 11);
        let iso_hi = pct_isotherm(&p, 350.0, 11);
        for (a, b) in iso_lo.iter().zip(iso_hi.iter()) {
            assert!(b.pressure_pa > a.pressure_pa);
        }
    }
}