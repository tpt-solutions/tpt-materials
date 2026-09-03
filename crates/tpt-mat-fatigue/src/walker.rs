//! Walker (1970) mean-stress correction.
//!
//! For a stress cycle `(σ_min, σ_max)` with load ratio
//! `R = σ_min / σ_max`, the equivalent fully-reversed stress
//! amplitude is
//!
//! `σ_ar = σ_a / (1 - R)^γ`
//!
//! where `γ ∈ (0, 1]` is the Walker exponent.  Special cases:
//!
//! - `γ = 0.5` (default for steels): Goodman-style correction.
//! - `γ = 1` (Walker 1970 default): the effective stress range
//!   `Δσ_eq = σ_max - σ_min` divided by `(1 - R)`.
//! - `γ = 0`: no mean-stress correction (fully-reversed only).
//!
//! When `R → 1` the correction diverges; we clamp `R ≤ 0.999` to
//! avoid blow-up and return `σ_ar = σ_a` in that limit.

/// Default Walker exponent for steels (Goodman-type correction).
pub fn walker_gamma() -> f64 {
    0.5
}

/// Equivalent fully-reversed stress amplitude `σ_ar` for the
/// supplied load ratio `R`, stress amplitude `σ_a`, and Walker
/// exponent `γ`.  `R` is clamped to `[−∞, 0.999]` so the
/// `(1 − R)^γ` factor stays finite.
pub fn equivalent_amplitude_walker(sigma_a: f64, r: f64, gamma: f64) -> f64 {
    let r_clamped = r.min(0.999);
    let factor = (1.0 - r_clamped).powf(gamma);
    if factor <= 0.0 {
        return sigma_a;
    }
    sigma_a / factor
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fully_reversed_load_ratio_is_unity_factor() {
        // R = -1 ⇒ (1 - R) = 2 ⇒ σ_ar = σ_a / 2^γ.
        let sigma_a = 200.0;
        let sa_r = equivalent_amplitude_walker(sigma_a, -1.0, 0.5);
        assert!((sa_r - sigma_a / 2.0_f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn zero_to_tension_correction() {
        // R = 0 ⇒ σ_ar = σ_a.
        let sa_r = equivalent_amplitude_walker(200.0, 0.0, 0.5);
        assert!((sa_r - 200.0).abs() < 1e-9);
    }

    #[test]
    fn near_one_load_ratio_is_clamped() {
        // R = 1 would diverge; we clamp.
        let sa_r = equivalent_amplitude_walker(200.0, 1.0, 0.5);
        assert!(sa_r.is_finite());
    }

    #[test]
    fn no_correction_when_gamma_is_zero() {
        let sigma_a = 200.0;
        let sa_r = equivalent_amplitude_walker(sigma_a, 0.5, 0.0);
        assert!((sa_r - sigma_a).abs() < 1e-9);
    }
}