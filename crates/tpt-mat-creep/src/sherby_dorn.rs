//! Sherby–Dorn (1958) temperature-compensated time parameter.
//!
//! `θ_D = t_r · exp(-Q / (R T))`
//!
//! where `Q` is the activation energy and `T` the absolute
//! temperature.  The Dorn parameter `θ_D` is a *time* that
//! collapses data from different temperatures onto a single
//! master curve `σ(θ_D)`.

use serde::{Deserialize, Serialize};

/// Sherby–Dorn parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SherbyDornParams {
    /// Activation energy `Q` (kJ/mol).
    pub activation_energy_kj_per_mol: f64,
    /// Universal gas constant `R kJ/(mol·K)`.
    pub r_kj_per_mol_k: f64,
}

impl Default for SherbyDornParams {
    fn default() -> Self {
        Self {
            activation_energy_kj_per_mol: 200.0,
            r_kj_per_mol_k: 8.314_462_618e-3,
        }
    }
}

/// Compute the Dorn parameter `θ_D = t_r · exp(-Q / (R T))`.
/// Returns `0` for `t_r ≤ 0` or `T ≤ 0`.
pub fn sherby_dorn_parameter(
    params: &SherbyDornParams,
    t_r: f64,
    t_kelvin: f64,
) -> f64 {
    if t_r <= 0.0 || t_kelvin <= 0.0 {
        return 0.0;
    }
    let exponent = -params.activation_energy_kj_per_mol / (params.r_kj_per_mol_k * t_kelvin);
    t_r * exponent.exp()
}

/// Inverse: rupture time `t_r = θ_D · exp(Q / (R T))`.
pub fn rupture_time_from_dorn(
    params: &SherbyDornParams,
    theta_d: f64,
    t_kelvin: f64,
) -> f64 {
    if theta_d <= 0.0 || t_kelvin <= 0.0 {
        return f64::INFINITY;
    }
    let exponent = params.activation_energy_kj_per_mol / (params.r_kj_per_mol_k * t_kelvin);
    theta_d * exponent.exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn dorn_round_trip() {
        let p = SherbyDornParams::default();
        let t = 1.0e4;
        let temp = 1000.0;
        let theta = sherby_dorn_parameter(&p, t, temp);
        let t_back = rupture_time_from_dorn(&p, theta, temp);
        assert_relative_eq!(t, t_back, max_relative = 1e-9);
    }

    #[test]
    fn dorn_increases_with_temperature() {
        let p = SherbyDornParams::default();
        let d_low = sherby_dorn_parameter(&p, 1.0e4, 800.0);
        let d_high = sherby_dorn_parameter(&p, 1.0e4, 1200.0);
        // Higher T ⇒ smaller Q/(RT) ⇒ larger exp(-Q/RT) ⇒ larger θ_D.
        assert!(d_high > d_low);
    }

    #[test]
    fn dorn_is_zero_for_nonpositive_input() {
        let p = SherbyDornParams::default();
        assert_eq!(sherby_dorn_parameter(&p, 0.0, 1000.0), 0.0);
        assert_eq!(sherby_dorn_parameter(&p, 1.0, 0.0), 0.0);
    }
}