//! Basquin (1910) high-cycle fatigue S-N law.
//!
//! `σ_a = σ_f' (2 N_f)^b`
//!
//! where `σ_a` is the stress amplitude, `σ_f'` the fatigue strength
//! coefficient (≈ fracture strength for many metals), `b` the
//! Basquin exponent (typically negative), and `N_f` the cycles to
//! failure.  Inverting gives
//!
//! `N_f = ½ (σ_a / σ_f')^(1/b)`.
//!
//! The **endurance limit** `σ_e` (Basquin strength at `N_e = 1e6` or
//! `1e7` cycles) is computed from the same parameters.

use serde::{Deserialize, Serialize};

/// Basquin S-N parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BasquinParams {
    /// Fatigue strength coefficient `σ_f'` (MPa).
    pub sigma_f_prime: f64,
    /// Fatigue strength exponent `b` (typically `-0.12 ≤ b ≤ -0.05`).
    pub b: f64,
    /// Reference cycles for the endurance limit (default `1e6`).
    pub endurance_cycles: f64,
}

impl Default for BasquinParams {
    fn default() -> Self {
        Self {
            sigma_f_prime: 1000.0,
            b: -0.1,
            endurance_cycles: 1.0e6,
        }
    }
}

/// Cycles to failure `N_f` for a fully-reversed stress amplitude
/// `σ_a`.  Returns `∞` for `σ_a ≤ 0`.  For `σ_a > σ_f'` the Basquin
/// law gives `N_f < 0.5` (i.e. failure in the first half-cycle),
/// which we clamp to `0.5`.
pub fn cycles_to_failure_basquin(params: &BasquinParams, sigma_a: f64) -> f64 {
    if sigma_a <= 0.0 {
        return f64::INFINITY;
    }
    if sigma_a >= params.sigma_f_prime {
        return 0.5;
    }
    let ratio = sigma_a / params.sigma_f_prime;
    0.5 * ratio.powf(1.0 / params.b)
}

/// Stress amplitude `σ_a` at `N_f` cycles.
pub fn fatigue_strength_at(params: &BasquinParams, n_f: f64) -> f64 {
    if n_f <= 0.0 {
        return params.sigma_f_prime;
    }
    let two_n = 2.0 * n_f;
    params.sigma_f_prime * two_n.powf(params.b)
}

/// Endurance limit `σ_e` at `params.endurance_cycles`.
pub fn endurance_limit(params: &BasquinParams) -> f64 {
    fatigue_strength_at(params, params.endurance_cycles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn basquin_round_trip() {
        let p = BasquinParams::default();
        let n_f = 1.0e5;
        let sa = fatigue_strength_at(&p, n_f);
        let n_f_2 = cycles_to_failure_basquin(&p, sa);
        assert_relative_eq!(n_f, n_f_2, max_relative = 1e-9);
    }

    #[test]
    fn basquin_strength_at_one_cycle_is_sigma_f_prime() {
        let p = BasquinParams::default();
        // At N_f = 0.5, 2 N_f = 1, so σ_a = σ_f'.
        let sa = fatigue_strength_at(&p, 0.5);
        assert_relative_eq!(sa, p.sigma_f_prime, epsilon = 1e-9);
    }

    #[test]
    fn endurance_limit_decreases_with_negative_b() {
        let mut p = BasquinParams::default();
        p.endurance_cycles = 1.0e6;
        let e0 = endurance_limit(&p);
        p.b = -0.05;
        let e1 = endurance_limit(&p);
        // Less negative b → shallower slope → higher endurance limit.
        assert!(e1 > e0);
    }

    #[test]
    fn cycles_to_failure_is_infinite_for_nonpositive_amplitude() {
        let p = BasquinParams::default();
        assert!(cycles_to_failure_basquin(&p, 0.0).is_infinite());
        assert!(cycles_to_failure_basquin(&p, -100.0).is_infinite());
    }

    #[test]
    fn cycles_to_failure_is_clamped_below_05_for_overload() {
        let p = BasquinParams::default();
        let n = cycles_to_failure_basquin(&p, p.sigma_f_prime * 2.0);
        assert!((n - 0.5).abs() < 1e-9);
    }
}