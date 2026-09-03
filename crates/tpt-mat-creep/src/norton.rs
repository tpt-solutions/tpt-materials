//! Norton–Bailey power-law creep.
//!
//! Steady-state creep rate (often called *secondary* or *minimum*
//! creep rate):
//!
//! `ε̇_ss = A · σ^n · exp(-Q / (R T))`
//!
//! where `A` is the creep constant, `σ` the uniaxial stress, `n`
//! the stress exponent (`3 ≤ n ≤ 8` for dislocation creep,
//! `n = 1` for Nabarro–Herring / diffusional creep, `n > 8` for
//! power-law breakdown), `Q` the activation energy, and `T` the
//! absolute temperature.

use serde::{Deserialize, Serialize};

/// Norton–Bailey parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NortonBaileyParams {
    /// Creep constant `A` (units `(1/s) / (MPa^n · exp(-Q/RT))` —
    /// the user is responsible for consistency).
    pub a: f64,
    /// Stress exponent `n` (typically 3–8 for metals).
    pub n: f64,
    /// Activation energy `Q` (kJ/mol).
    pub activation_energy_kj_per_mol: f64,
    /// Universal gas constant (kJ/(mol·K)).
    pub r_kj_per_mol_k: f64,
}

impl Default for NortonBaileyParams {
    fn default() -> Self {
        Self {
            a: 1.0e-5,
            n: 5.0,
            activation_energy_kj_per_mol: 200.0,
            r_kj_per_mol_k: 8.314_462_618e-3,
        }
    }
}

/// Steady-state creep rate `ε̇_ss` at stress `σ` (MPa) and
/// temperature `T` (K).  Returns `0` for `σ ≤ 0` or `T ≤ 0`.
pub fn norton_creep_rate(params: &NortonBaileyParams, sigma: f64, t_kelvin: f64) -> f64 {
    if sigma <= 0.0 || t_kelvin <= 0.0 {
        return 0.0;
    }
    let exponent = -params.activation_energy_kj_per_mol / (params.r_kj_per_mol_k * t_kelvin);
    params.a * sigma.powf(params.n) * exponent.exp()
}

/// Rupture time `t_r` estimated from the steady-state creep rate
/// via `t_r = ε_r / ε̇_ss` (engineering estimate).  Returns `∞`
/// for `ε̇_ss ≤ 0`.  `ε_r` is the rupture strain (typically
/// `0.05–0.30`).
pub fn rupture_time_norton(
    params: &NortonBaileyParams,
    sigma: f64,
    t_kelvin: f64,
    rupture_strain: f64,
) -> f64 {
    let rate = norton_creep_rate(params, sigma, t_kelvin);
    if rate <= 0.0 || rupture_strain <= 0.0 {
        return f64::INFINITY;
    }
    rupture_strain / rate
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn creep_rate_is_zero_below_zero_stress() {
        let p = NortonBaileyParams::default();
        assert_eq!(norton_creep_rate(&p, -1.0, 1000.0), 0.0);
        assert_eq!(norton_creep_rate(&p, 100.0, -1.0), 0.0);
    }

    #[test]
    fn creep_rate_scales_as_sigma_pow_n() {
        let p = NortonBaileyParams::default();
        let r1 = norton_creep_rate(&p, 100.0, 1000.0);
        let r2 = norton_creep_rate(&p, 200.0, 1000.0);
        let ratio = r2 / r1;
        let expected = 2.0_f64.powf(p.n);
        assert_relative_eq!(ratio, expected, max_relative = 1e-9);
    }

    #[test]
    fn creep_rate_increases_with_temperature() {
        let p = NortonBaileyParams::default();
        let r_low = norton_creep_rate(&p, 100.0, 800.0);
        let r_high = norton_creep_rate(&p, 100.0, 1200.0);
        assert!(r_high > r_low);
    }

    #[test]
    fn rupture_time_inverse_of_steady_state_rate() {
        let p = NortonBaileyParams::default();
        let rate = norton_creep_rate(&p, 100.0, 1000.0);
        let t_r = rupture_time_norton(&p, 100.0, 1000.0, 0.1);
        assert_relative_eq!(t_r, 0.1 / rate, epsilon = 1e-9);
    }
}
