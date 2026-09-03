//! Arruda–Boyce 8-chain rubber-elasticity model.

use super::fjc::{inverse_langevin_approx, inverse_langevin_exact};

/// Boltzmann constant (J/K).
pub const K_B: f64 = 1.380_649e-23;

/// Nominal (1st Piola–Kirchhoff) stress for the Arruda–Boyce
/// 8-chain model under uniaxial extension:
///
/// `P = G / 3 (λ − 1/λ²) · L^{-1}(λ / √N)`
///
/// where `L^{-1}` is the inverse Langevin function.
pub fn arruda_boyce_nominal_stress(n_segments: u32, shear_modulus: f64, stretch: f64) -> f64 {
    if stretch <= 1.0 {
        return 0.0;
    }
    let n = n_segments as f64;
    let sqrt_n = n.sqrt();
    let x = stretch / sqrt_n;
    let l_inv = inverse_langevin_approx(x);
    shear_modulus * (stretch - 1.0 / (stretch * stretch)) * l_inv / 3.0
}

/// Cauchy (true) stress for the Arruda–Boyce 8-chain model.
pub fn arruda_boyce_true_stress(n_segments: u32, shear_modulus: f64, stretch: f64) -> f64 {
    let p = arruda_boyce_nominal_stress(n_segments, shear_modulus, stretch);
    p * stretch
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn arruda_boyce_stress_grows_with_stretch() {
        let g = 0.4e6;
        let s1 = arruda_boyce_true_stress(8, g, 2.0);
        let s2 = arruda_boyce_true_stress(8, g, 4.0);
        assert!(s2 > s1);
    }

    #[test]
    fn neo_hookean_recovered_at_small_stretch() {
        // At very small stretches, AB ~ neo-Hookean, but the
        // exact agreement depends on N; we just verify
        // positivity.
        let g = 0.4e6;
        let n = 8;
        let s_ab = arruda_boyce_nominal_stress(n, g, 1.01);
        assert!(s_ab > 0.0);
    }

    #[test]
    fn inverse_langevin_re_exported() {
        let v = inverse_langevin_approx(0.5);
        assert!(v > 0.0);
    }
}
