//! Poroelastic Fickian uptake on a slab.

use serde::{Deserialize, Serialize};

/// Compute the Fickian solvent uptake on a slab of half-thickness
/// `L` (m) at time `t` (s).
///
/// `M(t)/M_∞ = 1 − (8/π²) Σ_{n=0}^∞ 1/((2n+1)²) exp(−D ((2n+1)π / 2L)² t)`
///
/// truncated to `n_terms`.
pub fn slab_fickian_uptake(
    diffusivity: f64,
    half_thickness: f64,
    time: f64,
    n_terms: usize,
) -> f64 {
    if half_thickness <= 0.0 {
        return 1.0;
    }
    let mut sum = 0.0;
    let factor = core::f64::consts::PI / (2.0 * half_thickness);
    for n in 0..n_terms {
        let k = (2 * n + 1) as f64;
        let exponent = -diffusivity * factor * factor * k * k * time;
        sum += (1.0 / (k * k)) * exponent.exp();
    }
    (1.0 - (8.0 / (core::f64::consts::PI * core::f64::consts::PI)) * sum).clamp(0.0, 1.0)
}

/// Convenience wrapper that returns the uptake curve as a
/// `Vec<(time, uptake)>` sampled at `n_samples` equally-spaced
/// points from `t = 0` to `t = t_final`.
pub fn poro_uptake_curve(
    diffusivity: f64,
    half_thickness: f64,
    t_final: f64,
    n_samples: usize,
) -> Vec<(f64, f64)> {
    if n_samples < 2 {
        return Vec::new();
    }
    let dt = t_final / ((n_samples - 1) as f64);
    (0..n_samples)
        .map(|i| {
            let t = (i as f64) * dt;
            let u = slab_fickian_uptake(diffusivity, half_thickness, t, 30);
            (t, u)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn uptake_zero_at_t_zero() {
        // At t = 0 the series collapses to π²/8 so the
        // truncated sum leaves a small positive residual.
        let u = slab_fickian_uptake(1.0e-9, 0.001, 0.0, 5000);
        assert!(u.abs() < 1.0e-4, "u = {u}");
    }

    #[test]
    fn uptake_approaches_one_at_long_times() {
        let u = slab_fickian_uptake(1.0e-9, 0.001, 1.0e5, 30);
        assert!(u > 0.99);
    }

    #[test]
    fn uptake_monotonically_increases() {
        let mut last = 0.0;
        for t in [0.0, 1.0e3, 1.0e4, 1.0e5, 1.0e6] {
            let u = slab_fickian_uptake(1.0e-9, 0.001, t, 20);
            assert!(u >= last);
            last = u;
        }
    }

    #[test]
    fn curve_has_correct_length() {
        let c = poro_uptake_curve(1.0e-9, 0.001, 1.0e5, 50);
        assert_eq!(c.len(), 50);
    }
}