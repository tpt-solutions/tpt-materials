//! θ-projection (Wilshire–Burt 1982) full creep-curve description.
//!
//! `ε(t) = θ_1 (1 - exp(-θ_2 t)) + θ_3 (exp(θ_4 t) - 1)`
//!
//! - `θ_1` (positive): the *primary* (decelerating) creep extent.
//! - `θ_2` (positive): the *primary* decay rate.
//! - `θ_3` (positive): the *tertiary* (accelerating) creep extent.
//! - `θ_4` (positive): the *tertiary* growth rate.
//!
//! The strain rate is then
//!
//! `ε̇(t) = θ_1 θ_2 exp(-θ_2 t) + θ_3 θ_4 exp(θ_4 t)`
//!
//! which is the sum of a decaying primary component and a growing
//! tertiary component.  Rupture occurs (approximately) at
//! `t_r ≈ ln(1 + 1/θ_3) / θ_4` for `θ_1 → 0`.

use serde::{Deserialize, Serialize};

/// θ-projection parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThetaProjectionParams {
    /// Primary creep extent `θ_1` (dimensionless).
    pub theta_1: f64,
    /// Primary decay rate `θ_2` (1/s).
    pub theta_2: f64,
    /// Tertiary creep extent `θ_3` (dimensionless).
    pub theta_3: f64,
    /// Tertiary growth rate `θ_4` (1/s).
    pub theta_4: f64,
}

impl Default for ThetaProjectionParams {
    fn default() -> Self {
        Self {
            theta_1: 0.02,
            theta_2: 0.1,
            theta_3: 0.05,
            theta_4: 0.001,
        }
    }
}

/// Creep strain `ε(t)` at time `t` (seconds).
pub fn theta_projection_creep_strain(params: &ThetaProjectionParams, t: f64) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    let primary = params.theta_1 * (1.0 - (-params.theta_2 * t).exp());
    let tertiary = params.theta_3 * ((params.theta_4 * t).exp() - 1.0);
    primary + tertiary
}

/// Creep strain rate `ε̇(t)` at time `t` (seconds).
pub fn theta_projection_creep_rate(params: &ThetaProjectionParams, t: f64) -> f64 {
    if t <= 0.0 {
        return params.theta_1 * params.theta_2;
    }
    let primary = params.theta_1 * params.theta_2 * (-params.theta_2 * t).exp();
    let tertiary = params.theta_3 * params.theta_4 * (params.theta_4 * t).exp();
    primary + tertiary
}

/// Approximate rupture time: the largest `t` for which
/// `ε̇(t) / ε̇(0) ≤ 100` (tertiary dominates by two orders of
/// magnitude).  A pragmatic choice; the precise definition
/// depends on the application.
pub fn rupture_time_theta(params: &ThetaProjectionParams) -> f64 {
    if params.theta_4 <= 0.0 {
        return f64::INFINITY;
    }
    let log_threshold = (100.0_f64).ln() / params.theta_4;
    log_threshold.max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn creep_strain_at_zero_time_is_zero() {
        let p = ThetaProjectionParams::default();
        assert_eq!(theta_projection_creep_strain(&p, 0.0), 0.0);
    }

    #[test]
    fn creep_strain_is_monotonic() {
        let p = ThetaProjectionParams::default();
        let e0 = theta_projection_creep_strain(&p, 100.0);
        let e1 = theta_projection_creep_strain(&p, 200.0);
        let e2 = theta_projection_creep_strain(&p, 1000.0);
        assert!(e0 <= e1);
        assert!(e1 <= e2);
    }

    #[test]
    fn creep_rate_drops_then_grows() {
        let p = ThetaProjectionParams::default();
        let r0 = theta_projection_creep_rate(&p, 0.0);
        let r_mid = theta_projection_creep_rate(&p, 500.0);
        // The primary decay reduces the rate initially.
        assert!(r_mid < r0);
        let r_late = theta_projection_creep_rate(&p, 5000.0);
        // Eventually tertiary dominates.
        assert!(r_late > r_mid);
    }

    #[test]
    fn rupture_time_is_finite() {
        let p = ThetaProjectionParams::default();
        let t_r = rupture_time_theta(&p);
        assert!(t_r.is_finite());
        assert!(t_r > 0.0);
    }
}