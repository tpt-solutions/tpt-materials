//! Nye tensor → scalar GND density.

/// Approximate scalar GND density from the Frobenius norm of the
/// Nye tensor `α_ij`:
// `ρ_GND = ||α|| / b`.
///
/// `nye_norm` is `√(Σ α_ij²)`.
pub fn gnd_from_curvature(nye_norm: f64, burgers_vector: f64) -> f64 {
    if burgers_vector.abs() < 1.0e-30 {
        return 0.0;
    }
    nye_norm.max(0.0) / burgers_vector
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn gnd_proportional_to_nye_norm() {
        let b = 2.5e-10;
        let r = gnd_from_curvature(1.0e-5, b);
        assert!(approx(r, 1.0e-5 / 2.5e-10, 1.0e-6));
    }

    #[test]
    fn zero_nye_yields_zero_gnd() {
        assert_eq!(gnd_from_curvature(0.0, 2.5e-10), 0.0);
    }

    #[test]
    fn zero_burgers_yields_zero() {
        assert_eq!(gnd_from_curvature(1.0, 0.0), 0.0);
    }
}
