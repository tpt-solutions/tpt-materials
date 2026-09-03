//! Lifshitz–Slyozov–Wagner (LSW) Ostwald-ripening coarsening law.

/// LSW coarsening rate constant
/// `K = 8 γ D c_eq / (9 V_m)` (for 3-D, interface-reaction
/// control ignored).
pub fn lsw_coarsening_rate(
    interfacial_energy: f64,
    diffusivity: f64,
    equilibrium_concentration: f64,
    molar_volume: f64,
) -> f64 {
    if molar_volume.abs() < 1.0e-30 {
        return 0.0;
    }
    8.0 * interfacial_energy * diffusivity * equilibrium_concentration
        / (9.0 * molar_volume)
}

/// Mean-radius-cubed growth `R̄³(t) − R̄³(0) = K t`.
///
/// Returns the increment `ΔR³` for time step `dt`.
pub fn lsw_radius_cubed_growth(k: f64, dt: f64) -> f64 {
    k * dt
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn rate_constant_grows_with_diffusivity() {
        let k1 = lsw_coarsening_rate(0.5, 1.0e-15, 0.01, 1.0e-5);
        let k2 = lsw_coarsening_rate(0.5, 1.0e-14, 0.01, 1.0e-5);
        assert!(k2 > k1);
    }

    #[test]
    fn rate_constant_zero_when_molar_volume_zero() {
        assert!(approx(lsw_coarsening_rate(0.5, 1.0e-15, 0.01, 0.0), 0.0, 1.0e-30));
    }

    #[test]
    fn radius_cubed_growth_is_linear_in_time() {
        let k = 1.0e-25;
        let dr3_1 = lsw_radius_cubed_growth(k, 1.0);
        let dr3_2 = lsw_radius_cubed_growth(k, 2.0);
        assert!(approx(dr3_2, 2.0 * dr3_1, 1.0e-30));
    }
}