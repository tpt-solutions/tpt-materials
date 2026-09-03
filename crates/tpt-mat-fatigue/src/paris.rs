//! Paris (1963) fatigue crack-growth law.
//!
//! `da/dN = C (ΔK)^m`
//!
//! where `a` is the crack length, `ΔK = K_max - K_min` the
//! stress-intensity range, and `C, m` the Paris parameters
//! (typically `m ∈ [2.5, 6]` for metals).
//!
//! `ΔK = Y σ √(π a)` for a centre-cracked plate of length `2a`
//! with stress range `σ` and geometry factor `Y ≈ 1` for
//! `a/W ≤ 0.5` (`W` = plate width).
//!
//! We expose three variants:
//! - `crack_growth_rate`: the bare Paris law.
//! - `paris_lifetime`: numerical integration from `a_0` to `a_c`.
//! - `critical_crack_length`: the crack length at which
//!   `ΔK = K_c` (fracture toughness).

use serde::{Deserialize, Serialize};

/// Paris-law parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ParisParams {
    /// Pre-exponential factor `C` (units `(m / cycle) / (MPa √m)^m`).
    pub c: f64,
    /// Paris exponent `m` (typically 2.5–6).
    pub m: f64,
    /// Geometry factor `Y` for the `ΔK` calculation.  Default `1.0`.
    pub y: f64,
    /// Fracture toughness `K_c` (MPa √m).  Used by
    /// `critical_crack_length`.
    pub k_c: f64,
}

impl Default for ParisParams {
    fn default() -> Self {
        Self {
            c: 1.0e-11,
            m: 3.0,
            y: 1.0,
            k_c: 50.0,
        }
    }
}

/// Paris crack-growth rate `da/dN` for a stress-intensity range
/// `ΔK` (MPa √m).  Returns `0` for `ΔK ≤ 0`.
pub fn crack_growth_rate(params: &ParisParams, delta_k: f64) -> f64 {
    if delta_k <= 0.0 {
        return 0.0;
    }
    params.c * delta_k.powf(params.m)
}

/// Critical crack length `a_c` at which `ΔK = K_c` for the supplied
/// stress range `Δσ`.  Returns `0` if `Δσ ≤ 0` or `Y = 0`.
pub fn critical_crack_length(params: &ParisParams, delta_sigma: f64) -> f64 {
    if delta_sigma <= 0.0 || params.y <= 0.0 {
        return 0.0;
    }
    let k = params.k_c;
    let k_over = k / (params.y * delta_sigma);
    k_over * k_over / std::f64::consts::PI
}

/// Number of cycles to grow a crack from `a_0` to `a_f` under a
/// constant stress range `Δσ`.  Uses numerical integration
/// (`N = ∫ da / (da/dN)`) with 1000 midpoint-rule samples.
pub fn paris_lifetime(
    params: &ParisParams,
    delta_sigma: f64,
    a_0: f64,
    a_f: f64,
    n_steps: usize,
) -> f64 {
    if a_f <= a_0 || delta_sigma <= 0.0 || params.y <= 0.0 {
        return f64::INFINITY;
    }
    let n_steps = n_steps.max(10);
    let da = (a_f - a_0) / n_steps as f64;
    let mut n = 0.0_f64;
    for i in 0..n_steps {
        let a_mid = a_0 + (i as f64 + 0.5) * da;
        let delta_k = params.y * delta_sigma * (std::f64::consts::PI * a_mid).sqrt();
        let dadn = crack_growth_rate(params, delta_k);
        if dadn <= 0.0 {
            continue;
        }
        n += da / dadn;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn crack_growth_rate_is_zero_below_threshold() {
        let p = ParisParams::default();
        assert_eq!(crack_growth_rate(&p, 0.0), 0.0);
        assert_eq!(crack_growth_rate(&p, -10.0), 0.0);
    }

    #[test]
    fn crack_growth_rate_matches_definition() {
        let p = ParisParams::default();
        let dk: f64 = 20.0;
        let expected = p.c * dk.powf(p.m);
        assert_relative_eq!(crack_growth_rate(&p, dk), expected, epsilon = 1e-12);
    }

    #[test]
    fn critical_crack_length_scales_with_inverse_sigma_squared() {
        let p = ParisParams::default();
        let a1 = critical_crack_length(&p, 100.0);
        let a2 = critical_crack_length(&p, 200.0);
        // a ∝ 1 / σ².
        assert!((a1 / a2 - 4.0).abs() < 1e-9);
    }

    #[test]
    fn paris_lifetime_is_finite_for_growth() {
        let p = ParisParams::default();
        let n = paris_lifetime(&p, 100.0, 1.0e-3, 5.0e-3, 1000);
        assert!(n.is_finite());
        assert!(n > 0.0);
    }

    #[test]
    fn paris_lifetime_is_infinite_for_decreasing_crack() {
        let p = ParisParams::default();
        let n = paris_lifetime(&p, 100.0, 5.0e-3, 1.0e-3, 1000);
        assert!(n.is_infinite());
    }
}
