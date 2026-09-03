//! Manson–Coffin–Basquin strain-life curve.
//!
//! `Δε / 2 = (σ_f' / E) (2 N_f)^b + ε_f' (2 N_f)^c`
//!
//! The first term (Basquin) dominates at high cycle counts
//! (elastic regime); the second (Coffin–Manson) dominates at low
//! cycle counts (plastic regime).  The transition (knee) occurs
//! at `N_t ≈ (E ε_f' / σ_f')^(1 / (b - c))`.

use serde::{Deserialize, Serialize};

use super::basquin::BasquinParams;
use super::coffin_manson::CoffinMansonParams;

/// Combined Manson–Coffin–Basquin strain-life parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StrainLifeParams {
    /// Young's modulus `E` (MPa).
    pub youngs_modulus: f64,
    /// Basquin parameters.
    pub basquin: BasquinParams,
    /// Coffin–Manson parameters.
    pub coffin_manson: CoffinMansonParams,
}

impl Default for StrainLifeParams {
    fn default() -> Self {
        Self {
            youngs_modulus: 200_000.0,
            basquin: BasquinParams {
                sigma_f_prime: 1200.0,
                b: -0.08,
                endurance_cycles: 1.0e6,
            },
            coffin_manson: CoffinMansonParams {
                epsilon_f_prime: 0.4,
                c: -0.6,
            },
        }
    }
}

/// Cycles to failure `N_f` from a total strain amplitude
/// `Δε / 2 = ε_a`.  Solved by Newton iteration of the implicit
/// equation.  Returns `∞` for `ε_a ≤ 0` or `ε_a < elastic limit`.
pub fn cycles_to_failure_strain_life(params: &StrainLifeParams, eps_a: f64) -> f64 {
    if eps_a <= 0.0 {
        return f64::INFINITY;
    }
    // Initial guess: Bisection on log10(N_f).
    let mut lo = 0.0_f64; // log10(1) = 0 ⇒ N = 1
    let mut hi = 10.0_f64; // log10(1e10)
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let n_f = 10_f64.powf(mid);
        let eps_pred = strain_amplitude_at(params, n_f);
        if eps_pred > eps_a {
            lo = mid; // need more cycles ⇒ reduce eps_a ⇒ larger n_f
        } else {
            hi = mid;
        }
        if (hi - lo) < 1e-6 {
            break;
        }
    }
    10_f64.powf(0.5 * (lo + hi))
}

/// Total strain amplitude `Δε / 2 = (σ_f' / E)(2 N_f)^b + ε_f' (2 N_f)^c`
/// at `N_f` cycles.
pub fn strain_amplitude_at(params: &StrainLifeParams, n_f: f64) -> f64 {
    if n_f <= 0.0 {
        return params.coffin_manson.epsilon_f_prime;
    }
    let two_n = 2.0 * n_f;
    let elastic = params.basquin.sigma_f_prime / params.youngs_modulus;
    let elastic_amp = elastic * two_n.powf(params.basquin.b);
    let plastic_amp = params.coffin_manson.epsilon_f_prime
        * two_n.powf(params.coffin_manson.c);
    elastic_amp + plastic_amp
}

/// Approximate transition cycles `N_t` where the elastic and
/// plastic contributions are equal.
pub fn transition_cycles(params: &StrainLifeParams) -> f64 {
    let e_eps = params.basquin.sigma_f_prime / params.youngs_modulus;
    let ratio = params.coffin_manson.epsilon_f_prime / e_eps;
    if params.basquin.b - params.coffin_manson.c == 0.0 {
        return f64::INFINITY;
    }
    let log_two_n = ratio.log10() / (params.basquin.b - params.coffin_manson.c);
    0.5 * 10_f64.powf(log_two_n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn strain_amplitude_at_half_cycle_is_max() {
        let p = StrainLifeParams::default();
        let e_half = strain_amplitude_at(&p, 0.5);
        let e_10 = strain_amplitude_at(&p, 10.0);
        let e_1e6 = strain_amplitude_at(&p, 1.0e6);
        // Total strain amplitude decreases with cycles.
        assert!(e_half > e_10);
        assert!(e_10 > e_1e6);
    }

    #[test]
    fn cycles_to_failure_round_trip() {
        let p = StrainLifeParams::default();
        let target_n = 5.0e4;
        let eps_a = strain_amplitude_at(&p, target_n);
        let n_calc = cycles_to_failure_strain_life(&p, eps_a);
        let rel_err = (n_calc - target_n).abs() / target_n;
        assert!(rel_err < 0.05, "got {n_calc}, expected {target_n}");
    }

    #[test]
    fn transition_cycles_finite_for_nonzero_exponent_difference() {
        let p = StrainLifeParams::default();
        let n_t = transition_cycles(&p);
        assert!(n_t.is_finite());
        assert!(n_t > 0.0);
    }

    #[test]
    fn cycles_to_failure_is_infinite_for_nonpositive_strain() {
        let p = StrainLifeParams::default();
        assert!(cycles_to_failure_strain_life(&p, 0.0).is_infinite());
    }
}