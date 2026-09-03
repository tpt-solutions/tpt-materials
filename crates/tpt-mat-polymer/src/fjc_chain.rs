//! Freely jointed chain force-extension.

use super::fjc::inverse_langevin_approx;

/// Force per chain `f = (1 / b) · L^{-1}(λ / (N b))`,
/// where `L^{-1}` is the inverse Langevin function.
pub fn fjc_force_extension(segment_length: f64, num_segments: u32, stretch: f64) -> f64 {
    if stretch <= 0.0 {
        return 0.0;
    }
    let n = num_segments as f64;
    let contour_length = n * segment_length;
    let x = (stretch * segment_length) / contour_length;
    let x = x.clamp(0.0, 0.999);
    inverse_langevin_approx(x) / segment_length
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn fjc_force_grows_with_stretch() {
        let f1 = fjc_force_extension(1.0e-9, 100, 1.1);
        let f2 = fjc_force_extension(1.0e-9, 100, 1.5);
        assert!(f2 > f1);
    }

    #[test]
    fn fjc_force_zero_at_zero_stretch() {
        assert!(approx(fjc_force_extension(1.0e-9, 100, 0.0), 0.0, 1.0e-12));
    }
}
