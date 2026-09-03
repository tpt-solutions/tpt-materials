//! Findley (1957) critical-plane fatigue-indicator parameter.
//!
//! On every candidate plane of normal direction `n = (cos θ, sin θ, 0)`
//! the shear stress amplitude `Δτ_max / 2` and the peak normal stress
//! `σ_n_max` are evaluated, and the Findley FIP is the maximum of
//! `τ_a + k σ_n_max` over the plane angle.
//!
//! For uniaxial loading `σ(t) = σ_max · s(t)` with load ratio `R`,
//! `σ_max = σ_a · (1 − R)`, the shear stress on a plane at angle `θ`
//! is `τ(θ) = (σ(t)/2) sin(2 θ)` so `Δτ/2 = |σ_max / 2 sin(2 θ)|`
//! and the normal stress on the plane is `σ_n(θ) = σ(t) cos²(θ)` so
//! `σ_n_max = |σ_max| cos²(θ)`.
//!
//! # Reference
//!
//! - Findley, W. N. (1957).  "Fatigue of metals under combined
//!   bending and torsion."  Proc. ASTM 57, 880–886.

/// Findley FIP for uniaxial loading.
///
/// `F = max_θ [ |σ_max / 2| · |sin 2θ| + k |σ_max| cos² θ ]`.
pub fn findley_fip(sigma_max: f64, sigma_min: f64, k: f64) -> f64 {
    let diff = (sigma_max - sigma_min).abs();
    let amp = 0.5 * diff;
    let peak = sigma_max.abs().max(sigma_min.abs());
    // Maximise f(θ) = amp · |sin 2θ| + k · peak · cos²θ over θ ∈ [0, π/2].
    // Sample 360 points to find the worst plane.
    let mut best = 0.0_f64;
    let n = 360;
    for i in 0..=n {
        let theta = core::f64::consts::FRAC_PI_2 * (i as f64) / (n as f64);
        let s2 = (2.0 * theta).sin().abs();
        let c = theta.cos();
        let f = amp * s2 + k * peak * c * c;
        if f > best {
            best = f;
        }
    }
    best
}

/// Findley FIP for a general (multiaxial) plane-stress state.
///
/// `σ_x`, `σ_y`, `τ_xy` are the in-plane principal stress components
/// in the plane of the RVE surface.
pub fn findley_from_multiaxial_stress(sigma_x: f64, sigma_y: f64, tau_xy: f64, k: f64) -> f64 {
    findley_search(sigma_x, sigma_y, tau_xy, k).fip
}

/// Search across candidate planes and return the maximum Findley FIP
/// along with the plane angle that achieves it.
pub fn findley_search(sigma_x: f64, sigma_y: f64, tau_xy: f64, k: f64) -> FindleyResult {
    let mut best = 0.0_f64;
    let mut best_theta = 0.0_f64;
    let n = 720;
    for i in 0..=n {
        let theta = core::f64::consts::PI * (i as f64) / (n as f64);
        let c = theta.cos();
        let s = theta.sin();
        let sn = c * c * sigma_x + s * s * sigma_y + 2.0 * c * s * tau_xy;
        let tn = (s * c * (sigma_y - sigma_x)) + (c * c - s * s) * tau_xy;
        let amp_t = tn.abs();
        let peak_sn = sn.abs();
        let f = amp_t + k * peak_sn;
        if f > best {
            best = f;
            best_theta = theta;
        }
    }
    FindleyResult {
        fip: best,
        theta: best_theta,
    }
}

/// Result of the Findley search over candidate planes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FindleyResult {
    /// The maximum FIP value.
    pub fip: f64,
    /// Plane angle (radians) that achieved the maximum.
    pub theta: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn findley_uniaxial_pure_tension() {
        // σ_max=100, σ_min=0, k=0.5.
        // F(θ) = 50|sin 2θ| + 0.5·100·cos²θ.
        // dF/dθ = 100 cos 2θ − 100 sin 2θ cos 2θ / ... = 0
        // ⇒ tan 2θ = 2 ⇒ θ ≈ 31.7° ⇒ F ≈ 80.90.
        let f = findley_fip(100.0, 0.0, 0.5);
        assert!(approx(f, 80.90, 1.0e-1));
    }

    #[test]
    fn findley_fully_reversed_higher_fip() {
        let f1 = findley_fip(100.0, 0.0, 0.3);
        let f2 = findley_fip(100.0, -100.0, 0.3);
        // Fully reversed ⇒ shear amplitude doubled ⇒ F higher.
        assert!(f2 > f1);
    }

    #[test]
    fn findley_search_max_shear_plane() {
        let fip = findley_from_multiaxial_stress(120.0, 50.0, -30.0, 0.5);
        assert!(fip > 0.0);
    }

    #[test]
    fn findley_k_zero_pure_shear() {
        // k=0 ⇒ F = amp · max|sin 2θ| = amp = 50.
        let f = findley_fip(100.0, 0.0, 0.0);
        assert!(approx(f, 50.0, 1.0e-3));
    }
}
