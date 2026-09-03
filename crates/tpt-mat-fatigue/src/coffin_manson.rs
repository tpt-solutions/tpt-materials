//! Coffin–Manson (1954/1953) low-cycle fatigue law.
//!
//! `Δε_p / 2 = ε_f' (2 N_f)^c`
//!
//! where `Δε_p / 2` is the plastic strain amplitude, `ε_f'` the
//! fatigue ductility coefficient (≈ the monotonic fracture strain),
//! `c` the Coffin–Manson exponent (typically `-0.7 ≤ c ≤ -0.5`),
//! and `N_f` the cycles to failure.

use serde::{Deserialize, Serialize};

/// Coffin–Manson parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CoffinMansonParams {
    /// Fatigue ductility coefficient `ε_f'` (dimensionless).
    pub epsilon_f_prime: f64,
    /// Fatigue ductility exponent `c` (typically negative).
    pub c: f64,
}

impl Default for CoffinMansonParams {
    fn default() -> Self {
        Self {
            epsilon_f_prime: 0.3,
            c: -0.6,
        }
    }
}

/// Cycles to failure from a plastic strain amplitude `Δε_p / 2`.
pub fn cycles_to_failure_coffin_manson(
    params: &CoffinMansonParams,
    plastic_strain_amplitude: f64,
) -> f64 {
    if plastic_strain_amplitude <= 0.0 {
        return f64::INFINITY;
    }
    if params.epsilon_f_prime <= 0.0 {
        return f64::INFINITY;
    }
    let ratio = plastic_strain_amplitude / params.epsilon_f_prime;
    0.5 * ratio.powf(1.0 / params.c)
}

/// Plastic strain amplitude `Δε_p / 2` at `N_f` cycles.
pub fn plastic_strain_at(params: &CoffinMansonParams, n_f: f64) -> f64 {
    if n_f <= 0.0 {
        return params.epsilon_f_prime;
    }
    let two_n = 2.0 * n_f;
    params.epsilon_f_prime * two_n.powf(params.c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn coffin_manson_round_trip() {
        let p = CoffinMansonParams::default();
        let n_f = 500.0;
        let epa = plastic_strain_at(&p, n_f);
        let n_f_2 = cycles_to_failure_coffin_manson(&p, epa);
        assert_relative_eq!(n_f, n_f_2, max_relative = 1e-9);
    }

    #[test]
    fn plastic_strain_at_half_cycle_is_epsilon_f_prime() {
        let p = CoffinMansonParams::default();
        let epa = plastic_strain_at(&p, 0.5);
        assert_relative_eq!(epa, p.epsilon_f_prime, epsilon = 1e-9);
    }

    #[test]
    fn cycles_to_failure_is_infinite_for_nonpositive_input() {
        let p = CoffinMansonParams::default();
        assert!(cycles_to_failure_coffin_manson(&p, 0.0).is_infinite());
        assert!(cycles_to_failure_coffin_manson(&p, -0.01).is_infinite());
    }

    #[test]
    fn plastic_strain_decreases_with_increasing_n() {
        let p = CoffinMansonParams::default();
        let e1 = plastic_strain_at(&p, 10.0);
        let e2 = plastic_strain_at(&p, 1000.0);
        // c < 0 ⇒ more cycles ⇒ smaller allowable plastic strain.
        assert!(e1 > e2);
    }
}
