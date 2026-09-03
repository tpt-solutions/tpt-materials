//! Kachanov (1958) creep-damage law.
//!
//! `Ḋ = A σ^n / (1 - D)^k`
//!
//! with rupture at `D = 1`.  This is the classical creep-damage
//! model for creep rupture (tertiary creep dominated by void
//! growth).

use serde::{Deserialize, Serialize};

/// Kachanov creep-damage parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct KachanovParams {
    /// Damage-rate constant `A`.
    pub a: f64,
    /// Stress exponent `n` (typically 3–8).
    pub n: f64,
    /// Damage-coupling exponent `k` (typically 1–3).
    pub k: f64,
}

impl Default for KachanovParams {
    fn default() -> Self {
        Self {
            a: 1.0e-12,
            n: 5.0,
            k: 2.0,
        }
    }
}

/// Damage rate `Ḋ` at stress `σ` and current damage `D`.
pub fn kachanov_damage_rate(params: &KachanovParams, sigma: f64, damage: f64) -> f64 {
    if damage >= 1.0 {
        return f64::INFINITY;
    }
    let one_minus_d = 1.0 - damage;
    if one_minus_d <= 0.0 {
        return f64::INFINITY;
    }
    params.a * sigma.powf(params.n) * one_minus_d.powf(-params.k)
}

/// Closed-form rupture time `t_r = (1 + k) / (A σ^n)` derived by
/// integrating `Ḋ` from `D = 0` to `D = 1`.  Rupture requires
/// `A σ^n > 0`.
pub fn kachanov_rupture_time(params: &KachanovParams, sigma: f64) -> f64 {
    if sigma <= 0.0 {
        return f64::INFINITY;
    }
    let stress_term = params.a * sigma.powf(params.n);
    if stress_term <= 0.0 {
        return f64::INFINITY;
    }
    (1.0 + params.k) / stress_term
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn damage_rate_is_infinite_at_rupture() {
        let p = KachanovParams::default();
        assert!(kachanov_damage_rate(&p, 100.0, 1.0).is_infinite());
        assert!(kachanov_damage_rate(&p, 100.0, 1.5).is_infinite());
    }

    #[test]
    fn damage_rate_increases_with_damage() {
        let p = KachanovParams::default();
        let r0 = kachanov_damage_rate(&p, 100.0, 0.0);
        let r5 = kachanov_damage_rate(&p, 100.0, 0.5);
        let r9 = kachanov_damage_rate(&p, 100.0, 0.9);
        assert!(r5 > r0);
        assert!(r9 > r5);
    }

    #[test]
    fn rupture_time_matches_closed_form() {
        let p = KachanovParams { a: 1.0e-6, n: 3.0, k: 2.0 };
        let t_r = kachanov_rupture_time(&p, 100.0);
        let expected = (1.0 + 2.0) / (1.0e-6 * 100.0_f64.powf(3.0));
        assert_relative_eq!(t_r, expected, max_relative = 1e-9);
    }

    #[test]
    fn rupture_time_is_infinite_for_zero_stress() {
        let p = KachanovParams::default();
        assert!(kachanov_rupture_time(&p, 0.0).is_infinite());
    }
}