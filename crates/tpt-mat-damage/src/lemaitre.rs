//! Lemaitre (1971) ductile continuum-damage law.
//!
//! `Ḋ = (σ̃_eq / σ_u)^s · (1 - D)^(-α) / B`
//!
//! where `σ̃_eq = σ_eq / (1 - D)` is the *effective* (undamaged)
//! equivalent stress, `σ_u` a reference stress, `s, α, B`
//! material parameters.  Rupture occurs at `D = 1`.

use serde::{Deserialize, Serialize};

/// Lemaitre damage-law parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LemaitreParams {
    /// Reference stress `σ_u` (MPa).
    pub sigma_u: f64,
    /// Stress exponent `s` (typically > 1).
    pub s: f64,
    /// Damage-coupling exponent `α` (typically ≥ 0).
    pub alpha: f64,
    /// Time constant `B`.
    pub b: f64,
}

impl Default for LemaitreParams {
    fn default() -> Self {
        Self {
            sigma_u: 500.0,
            s: 2.0,
            alpha: 1.0,
            b: 1.0e4,
        }
    }
}

/// Damage rate `Ḋ` at the current damage state.
pub fn lemaitre_damage_rate(
    params: &LemaitreParams,
    sigma_eq: f64,
    damage: f64,
) -> f64 {
    if damage >= 1.0 {
        return f64::INFINITY;
    }
    let one_minus_d = 1.0 - damage;
    if one_minus_d <= 0.0 {
        return f64::INFINITY;
    }
    let sigma_eff = sigma_eq / one_minus_d;
    let stress_term = (sigma_eff / params.sigma_u).powf(params.s);
    let damage_term = one_minus_d.powf(-params.alpha);
    stress_term * damage_term / params.b
}

/// Forward-Euler advance of the damage state by a time step `Δt`.
/// Clamps `D ∈ [0, 1)` and returns the new damage.
pub fn lemaitre_damage_step(
    params: &LemaitreParams,
    sigma_eq: f64,
    damage: f64,
    dt: f64,
) -> f64 {
    if damage >= 1.0 {
        return 1.0;
    }
    let rate = lemaitre_damage_rate(params, sigma_eq, damage);
    let new_d = damage + rate * dt;
    new_d.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn damage_rate_is_infinite_at_rupture() {
        let p = LemaitreParams::default();
        assert!(lemaitre_damage_rate(&p, 100.0, 1.0).is_infinite());
        assert!(lemaitre_damage_rate(&p, 100.0, 1.5).is_infinite());
    }

    #[test]
    fn damage_step_clamps_to_unit_at_rupture() {
        let p = LemaitreParams::default();
        let new_d = lemaitre_damage_step(&p, 200.0, 0.99, 1.0e3);
        assert!(new_d <= 1.0);
        assert!(new_d >= 0.99);
    }

    #[test]
    fn damage_grows_monotonically_under_constant_stress() {
        let p = LemaitreParams::default();
        let mut d = 0.0;
        for _ in 0..10 {
            d = lemaitre_damage_step(&p, 200.0, d, 1.0);
            if d >= 1.0 {
                break;
            }
        }
        assert!(d > 0.0);
    }

    #[test]
    fn damage_rate_is_zero_for_zero_stress() {
        let p = LemaitreParams::default();
        assert_eq!(lemaitre_damage_rate(&p, 0.0, 0.0), 0.0);
    }
}