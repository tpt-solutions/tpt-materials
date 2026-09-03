//! Smith–Watson–Topper (SWT, 1970) fatigue-indicator parameter.
//!
//! For uniaxial loading
//!
//! ```text
//! F = σ_max · (Δε / 2)
//! ```
//!
//! where `σ_max = (1 − R) / (1 + R) · σ_a · 2` for a sinusoidal
//! history, but the helper accepts the raw peak stress and total
//! strain range so any load history can be supplied.
//!
//! # Reference
//!
//! - Smith, K. N., Watson, P., & Topper, T. H. (1970).  "A
//!   stress-strain function for the fatigue of metals."  J. Mater.
//!   5(4), 767–778.

/// Smith–Watson–Topper FIP: `F = σ_max · Δε / 2`.
pub fn swt_fip(sigma_max: f64, sigma_min: f64, delta_eps_total: f64) -> f64 {
    let peak = sigma_max.abs().max(sigma_min.abs());
    0.5 * peak * delta_eps_total
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn swt_fully_reversed_zero_mean() {
        let f = swt_fip(100.0, -100.0, 0.01);
        assert!(approx(f, 0.5, 1.0e-9));
    }

    #[test]
    fn swt_tensile_mean_higher() {
        let f_zero = swt_fip(100.0, -100.0, 0.01);
        let f_tens = swt_fip(150.0, -50.0, 0.01);
        assert!(f_tens > f_zero);
    }

    #[test]
    fn swt_compressive_mean_higher_peak() {
        // SWT uses the peak |σ|, so a compressive mean with higher
        // absolute amplitudes raises F.
        let f_zero = swt_fip(100.0, -100.0, 0.01);
        let f_comp = swt_fip(50.0, -150.0, 0.01);
        assert!(f_comp > f_zero);
    }
}
