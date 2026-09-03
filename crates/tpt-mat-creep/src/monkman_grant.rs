//! Monkman–Grant (1956) empirical relation between steady-state
//! creep rate and rupture time.
//!
//! `t_r · ε̇_ss^m = C_mg`
//!
//! where `m ∈ [0.7, 1.0]` (typically `m ≈ 0.9`) and `C_mg` is the
//! Monkman–Grant constant.  Useful for extrapolating short-term
//! creep-rate data to long rupture times.

use serde::{Deserialize, Serialize};

/// Monkman–Grant parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MonkmanGrantParams {
    /// Monkman–Grant constant `C_mg`.
    pub c_mg: f64,
    /// Monkman–Grant exponent `m` (typically `≈ 0.9`).
    pub m: f64,
}

impl Default for MonkmanGrantParams {
    fn default() -> Self {
        Self {
            c_mg: 0.05,
            m: 0.9,
        }
    }
}

/// Compute the Monkman–Grant product `t_r · ε̇_ss^m` for the
/// supplied rupture time and steady-state creep rate.  For
/// consistency with `t_r · ε̇_ss^m = C_mg`, the function returns
/// the *left-hand side*.
pub fn monkman_grant_product(
    _params: &MonkmanGrantParams,
    steady_state_creep_rate: f64,
    rupture_time: f64,
) -> f64 {
    if steady_state_creep_rate <= 0.0 || rupture_time <= 0.0 {
        return 0.0;
    }
    steady_state_creep_rate.powf(_params.m) * rupture_time
}

/// Rupture time `t_r = (C_mg / ε̇_ss^m)^(1 / (1 - m))` from a
/// known steady-state creep rate.  Returns `∞` if `ε̇_ss ≤ 0` or
/// `m ≥ 1`.
pub fn rupture_time_monkman_grant(
    params: &MonkmanGrantParams,
    steady_state_creep_rate: f64,
) -> f64 {
    if steady_state_creep_rate <= 0.0 {
        return f64::INFINITY;
    }
    if params.m >= 1.0 {
        return f64::INFINITY;
    }
    let ratio = params.c_mg / steady_state_creep_rate.powf(params.m);
    ratio.powf(1.0 / (1.0 - params.m))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn rupture_time_recovers_c_mg_at_unit_rate() {
        let p = MonkmanGrantParams::default();
        let t_r = rupture_time_monkman_grant(&p, 1.0);
        // At ε̇_ss = 1, t_r = (C_mg)^(1/(1-m)).
        let expected = p.c_mg.powf(1.0 / (1.0 - p.m));
        assert_relative_eq!(t_r, expected, max_relative = 1e-9);
    }

    #[test]
    fn rupture_time_increases_with_decreasing_creep_rate() {
        let p = MonkmanGrantParams::default();
        let t_r_high = rupture_time_monkman_grant(&p, 1.0e-3);
        let t_r_low = rupture_time_monkman_grant(&p, 1.0e-5);
        assert!(t_r_low > t_r_high);
    }

    #[test]
    fn rupture_time_is_infinite_for_nonpositive_rate() {
        let p = MonkmanGrantParams::default();
        assert!(rupture_time_monkman_grant(&p, 0.0).is_infinite());
        assert!(rupture_time_monkman_grant(&p, -1.0e-6).is_infinite());
    }

    #[test]
    fn rupture_time_is_infinite_for_m_ge_1() {
        let p = MonkmanGrantParams { c_mg: 1.0, m: 1.0 };
        assert!(rupture_time_monkman_grant(&p, 1.0e-3).is_infinite());
    }
}