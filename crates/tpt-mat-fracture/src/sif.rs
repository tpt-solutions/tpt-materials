//! Stress intensity factors for classical crack geometries.

use serde::{Deserialize, Serialize};

/// Mode mixity tag for a stress-intensity factor.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Mode {
    /// Opening mode (`K_I`).
    Opening,
    /// In-plane shear (`K_II`).
    ShearInPlane,
    /// Out-of-plane (anti-plane) shear (`K_III`).
    ShearOutOfPlane,
}

/// A stress-intensity factor value associated with a crack mode.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StressIntensityFactor {
    /// Crack mode.
    pub mode: Mode,
    /// Stress-intensity factor magnitude (Pa·m^0.5).
    pub value: f64,
}

/// Stress intensity factor for a centre crack of length `2a` in
/// an infinite plate under remote tension `σ`:
///
/// `K_I = σ √(π a)`.
///
/// Returns 0 for `a < 0`.
pub fn k_i_centre_crack(stress: f64, half_crack_length: f64) -> f64 {
    if half_crack_length <= 0.0 {
        return 0.0;
    }
    stress * (core::f64::consts::PI * half_crack_length).sqrt()
}

/// Stress intensity factor for a single-edge crack of length `a`
/// in a semi-infinite plate under remote tension `σ`:
///
/// `K_I = 1.12 σ √(π a)`.
pub fn k_i_edge_crack(stress: f64, crack_length: f64) -> f64 {
    if crack_length <= 0.0 {
        return 0.0;
    }
    1.12 * stress * (core::f64::consts::PI * crack_length).sqrt()
}

/// Mode-II stress intensity factor for a centre crack under remote
/// in-plane shear `τ`:
///
/// `K_II = τ √(π a)`.
pub fn k_ii_centre_crack(shear_stress: f64, half_crack_length: f64) -> f64 {
    if half_crack_length <= 0.0 {
        return 0.0;
    }
    shear_stress * (core::f64::consts::PI * half_crack_length).sqrt()
}

/// Mode-III stress intensity factor for a penny-shaped crack under
/// remote anti-plane shear `τ`:
///
/// `K_III = (2/π) τ √(π a)`.
pub fn k_iii_centre_crack(shear_stress: f64, half_crack_length: f64) -> f64 {
    if half_crack_length <= 0.0 {
        return 0.0;
    }
    (2.0 / core::f64::consts::PI)
        * shear_stress
        * (core::f64::consts::PI * half_crack_length).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn centre_crack_k_grows_with_crack_length() {
        let k1 = k_i_centre_crack(100.0e6, 0.01);
        let k2 = k_i_centre_crack(100.0e6, 0.04);
        assert!(k2 > k1);
    }

    #[test]
    fn centre_crack_k_scales_with_stress() {
        let k1 = k_i_centre_crack(100.0e6, 0.01);
        let k2 = k_i_centre_crack(200.0e6, 0.01);
        assert!(approx(k2, 2.0 * k1, 1.0e-6));
    }

    #[test]
    fn zero_crack_returns_zero() {
        assert_eq!(k_i_centre_crack(100.0e6, 0.0), 0.0);
        assert_eq!(k_i_centre_crack(100.0e6, -0.01), 0.0);
        assert_eq!(k_i_edge_crack(100.0e6, -0.01), 0.0);
    }

    #[test]
    fn edge_crack_higher_than_centre_crack() {
        let stress = 100.0e6;
        let a = 0.01;
        assert!(k_i_edge_crack(stress, a) > k_i_centre_crack(stress, a));
    }
}