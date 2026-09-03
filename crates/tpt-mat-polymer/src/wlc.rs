//! Worm-like chain (Marko–Siggia) force-extension.

/// WLC force per chain in the Marko–Siggia interpolation:
/// `f = (k_B T / l_p) [1/(4 (1 − λ/L_c)²) − 1/4 + λ/L_c]`
/// for `λ ≤ L_c`.  Beyond `L_c`, force diverges.
pub fn wlc_force_extension(persistence_length: f64, contour_length: f64, stretch: f64) -> f64 {
    if stretch <= 0.0 {
        return 0.0;
    }
    if stretch >= contour_length {
        return f64::INFINITY;
    }
    let zeta = stretch / contour_length;
    let one_minus = 1.0 - zeta;
    if one_minus <= 0.0 {
        return f64::INFINITY;
    }
    1.0 / (4.0 * persistence_length * one_minus * one_minus) - 1.0 / (4.0 * persistence_length)
        + stretch / (persistence_length * contour_length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wlc_force_diverges_at_full_extension() {
        let f = wlc_force_extension(0.5e-9, 5.0e-8, 5.0e-8);
        assert!(f.is_infinite());
    }

    #[test]
    fn wlc_force_grows_with_stretch() {
        let f1 = wlc_force_extension(0.5e-9, 5.0e-8, 0.5e-8);
        let f2 = wlc_force_extension(0.5e-9, 5.0e-8, 1.0e-8);
        assert!(f2 > f1);
    }
}
