//! Larson–Miller (1952) time-temperature parameter.
//!
//! `LMP = T (log t_r + C_LM)`
//!
//! where `t_r` is the rupture time (hours), `T` the absolute
//! temperature (Kelvin), and `C_LM ≈ 20` the Larson–Miller
//! constant.  The LMP is used to extrapolate short-term
//! high-temperature creep data to long rupture times at lower
//! temperatures by collecting all data on a single master curve
//! `LMP(σ)` of stress vs. LMP.

use serde::{Deserialize, Serialize};

/// Larson–Miller parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LarsonMillerParams {
    /// Larson–Miller constant `C_LM` (typically 20).
    pub c_lm: f64,
}

impl Default for LarsonMillerParams {
    fn default() -> Self {
        Self { c_lm: 20.0 }
    }
}

/// Compute the Larson–Miller parameter `LMP = T (log t_r + C_LM)`
/// from a rupture time `t_r` in *hours* (the historical unit) and
/// temperature `T` (K).  Returns `0` for `t_r ≤ 0` or `T ≤ 0`.
pub fn larson_miller_parameter(
    params: &LarsonMillerParams,
    t_hours: f64,
    t_kelvin: f64,
) -> f64 {
    if t_hours <= 0.0 || t_kelvin <= 0.0 {
        return 0.0;
    }
    t_kelvin * (t_hours.log10() + params.c_lm)
}

/// Inverse: rupture time `t_r` from LMP and temperature.
pub fn rupture_time_from_lmp(
    params: &LarsonMillerParams,
    lmp: f64,
    t_kelvin: f64,
) -> f64 {
    if t_kelvin <= 0.0 {
        return f64::INFINITY;
    }
    let log_t = lmp / t_kelvin - params.c_lm;
    if log_t <= 0.0 {
        return f64::INFINITY;
    }
    10_f64.powf(log_t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn lmp_round_trip() {
        let p = LarsonMillerParams::default();
        let t = 1000.0_f64; // hours
        let temp = 873.0_f64; // 600 °C in K
        let lmp = larson_miller_parameter(&p, t, temp);
        let t_back = rupture_time_from_lmp(&p, lmp, temp);
        assert_relative_eq!(t, t_back, max_relative = 1e-9);
    }

    #[test]
    fn lmp_increases_with_temperature() {
        let p = LarsonMillerParams::default();
        let t = 1000.0;
        let l1 = larson_miller_parameter(&p, t, 800.0);
        let l2 = larson_miller_parameter(&p, t, 1000.0);
        assert!(l2 > l1);
    }

    #[test]
    fn lmp_increases_with_rupture_time() {
        let p = LarsonMillerParams::default();
        let l1 = larson_miller_parameter(&p, 100.0, 1000.0);
        let l2 = larson_miller_parameter(&p, 10_000.0, 1000.0);
        assert!(l2 > l1);
    }

    #[test]
    fn lmp_is_zero_for_nonpositive_input() {
        let p = LarsonMillerParams::default();
        assert_eq!(larson_miller_parameter(&p, 0.0, 1000.0), 0.0);
        assert_eq!(larson_miller_parameter(&p, 1.0, 0.0), 0.0);
    }
}