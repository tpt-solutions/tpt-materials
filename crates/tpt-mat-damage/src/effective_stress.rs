//! Effective-stress concept and stress triaxiality.

/// Effective Cauchy stress `σ̃ = σ / (1 - D)` for a scalar damage
/// variable `D ∈ [0, 1]`.  Returns `∞` for `D ≥ 1` (rupture).
pub fn effective_stress(sigma: f64, damage: f64) -> f64 {
    let denom = 1.0 - damage;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    sigma / denom
}

/// Effective Young's modulus `Ẽ = E (1 - D)` (Lemaitre 1971).
pub fn effective_youngs_modulus(youngs: f64, damage: f64) -> f64 {
    if damage < 0.0 || damage > 1.0 {
        return 0.0;
    }
    youngs * (1.0 - damage)
}

/// Stress triaxiality `η = σ_m / σ_eq` where `σ_m = tr(σ)/3` is
/// the hydrostatic stress and `σ_eq` the von-Mises equivalent
/// (for principal stresses `σ_1, σ_2, σ_3`).
pub fn stress_triaxiality(sigma_1: f64, sigma_2: f64, sigma_3: f64) -> f64 {
    let sigma_m = (sigma_1 + sigma_2 + sigma_3) / 3.0;
    let s1 = sigma_1 - sigma_m;
    let s2 = sigma_2 - sigma_m;
    let s3 = sigma_3 - sigma_m;
    let sigma_eq = (1.5 * (s1 * s1 + s2 * s2 + s3 * s3)).sqrt();
    if sigma_eq == 0.0 {
        return 0.0;
    }
    sigma_m / sigma_eq
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn effective_stress_is_infinite_at_rupture() {
        assert!(effective_stress(100.0, 1.0).is_infinite());
        assert!(effective_stress(100.0, 1.5).is_infinite());
    }

    #[test]
    fn effective_stress_amplifies_stress_with_damage() {
        let e0 = effective_stress(100.0, 0.0);
        let e5 = effective_stress(100.0, 0.5);
        let e9 = effective_stress(100.0, 0.9);
        assert!((e0 - 100.0).abs() < 1e-9);
        assert!((e5 - 200.0).abs() < 1e-9);
        assert!((e9 - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn effective_youngs_modulus_decreases_with_damage() {
        let e0 = effective_youngs_modulus(200_000.0, 0.0);
        let e5 = effective_youngs_modulus(200_000.0, 0.5);
        assert!((e0 - 200_000.0).abs() < 1e-9);
        assert!((e5 - 100_000.0).abs() < 1e-9);
    }

    #[test]
    fn stress_triaxiality_is_two_thirds_for_equibiaxial_tension() {
        let eta = stress_triaxiality(100.0, 100.0, 0.0);
        // σ_m = 200/3, σ_eq = 100 ⇒ η = 2/3.
        assert!((eta - 2.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn stress_triaxiality_is_zero_for_uniaxial_tension() {
        let eta = stress_triaxiality(100.0, 0.0, 0.0);
        assert!((eta - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn stress_triaxiality_is_zero_for_pure_shear() {
        let eta = stress_triaxiality(100.0, -100.0, 0.0);
        assert!(eta.abs() < 1e-9);
    }
}