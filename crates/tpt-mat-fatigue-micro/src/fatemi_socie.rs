//! Fatemi–Socie (1988) critical-plane fatigue-indicator parameter.
//!
//! For each candidate plane the maximum shear strain amplitude
//! `Δγ_max / 2` is combined with the peak normal stress `σ_n_max`
//! via
//!
//! ```text
//! F = (Δγ_max / 2) · (1 + k · ⟨σ_n_max⟩ / σ_y)
//! ```
//!
//! where `⟨x⟩ = max(x, 0)` is the McCauley bracket (only tensile
//! normal stresses count).  For uniaxial tension the worst plane is
//! at 45° so `Δγ_max / 2 = Δε / 2` and the FIP reduces to
//! `(Δε/2)(1 + k σ_max/σ_y)`.
//!
//! # Reference
//!
//! - Fatemi, A., & Socie, D. F. (1988).  "A critical plane approach
//!   to multiaxial fatigue damage including out-of-phase loading."
//!   Fatigue Fract. Eng. Mater. Struct. 11(3), 149–165.

/// Fatemi–Socie FIP for uniaxial loading.
///
/// The first argument is the total cyclic strain *range* `Δε`; the
/// formula is
///
/// ```text
/// F = (Δε / 2) · (1 + k · max(σ_max, 0) / σ_y)
/// ```
pub fn fatemi_socie_fip(delta_eps_total: f64, sigma_max: f64, sigma_y: f64, k: f64) -> f64 {
    let s_n_pos = if sigma_max > 0.0 { sigma_max } else { 0.0 };
    0.5 * delta_eps_total * (1.0 + k * s_n_pos / sigma_y)
}

/// Fatemi–Socie FIP for a multiaxial plane-stress state.
///
/// Returns the maximum FIP across 360 candidate plane angles.
pub fn fatemi_socie_search(
    delta_eps_x: f64,
    delta_eps_y: f64,
    delta_eps_xy: f64,
    sigma_x_max: f64,
    sigma_y_max: f64,
    tau_xy_max: f64,
    sigma_y: f64,
    k: f64,
) -> f64 {
    let mut best = 0.0_f64;
    let n = 360;
    for i in 0..=n {
        let theta = core::f64::consts::PI * (i as f64) / (n as f64);
        let c = theta.cos();
        let s = theta.sin();
        let tg_amp = (s * c * (delta_eps_y - delta_eps_x)) + (c * c - s * s) * delta_eps_xy;
        let gamma_amp = tg_amp.abs();
        let sn_normal = c * c * sigma_x_max + s * s * sigma_y_max + 2.0 * c * s * tau_xy_max;
        let s_n_pos = if sn_normal > 0.0 { sn_normal } else { 0.0 };
        let f = gamma_amp * (1.0 + k * s_n_pos / sigma_y);
        if f > best {
            best = f;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn fatemi_socie_zero_mean_stress_matches_shear_amplitude() {
        let f = fatemi_socie_fip(0.01, 0.0, 250.0, 0.5);
        assert!(approx(f, 0.005, 1.0e-9));
    }

    #[test]
    fn fatemi_socie_tensile_mean_amplifies_fip() {
        let f_zero = fatemi_socie_fip(0.01, 0.0, 250.0, 0.5);
        let f_tens = fatemi_socie_fip(0.01, 100.0, 250.0, 0.5);
        assert!(f_tens > f_zero);
    }

    #[test]
    fn fatemi_socie_k_zero_is_pure_shear_amplitude() {
        let f = fatemi_socie_fip(0.005, 50.0, 250.0, 0.0);
        assert!(approx(f, 0.0025, 1.0e-9));
    }
}
