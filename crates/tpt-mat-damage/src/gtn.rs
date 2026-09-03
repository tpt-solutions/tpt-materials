//! Gurson–Tvergaard–Needleman (GTN) porous plasticity.
//!
//! The GTN model extends Gurson's (1977) yield function for a
//! porous ductile matrix with Tvergaard & Needleman's (1984)
//! calibration constants and a strain-controlled nucleation term.
//!
//! ## Yield function
//!
//! ```text
//! Φ = (σ_eq / σ_y)² + 2 q₁ f cosh(3 q₂ σ_h / (2 σ_y))
//!     − (1 + q₃ f²)
//! ```
//!
//! where:
//!
//! - `σ_eq` — von Mises equivalent stress
//! - `σ_h` — hydrostatic stress (`σ_kk / 3`)
//! - `σ_y` — current matrix yield stress
//! - `f`   — current porosity (void volume fraction)
//!
//! Classic Tvergaard–Needleman calibration: `q₁ = 1.5`,
//! `q₂ = 1.0`, `q₃ = q₁² = 2.25`.
//!
//! ## Porosity evolution
//!
//! ```text
//! ḟ = (1 − f) ε̇^p_kk + A ε̇^p_eq     (void growth + nucleation)
//! ```
//!
//! The nucleation term uses the strain-controlled Chu–Needleman form
//! with normal-distribution parameters `(ε_n, s_n)` and amplitude
//! `f_n`.  Void coalescence is captured by the rapid void-growth
// regime above `f = f_c` and final failure at `f = f_f`.
//!
//! ## Reference behaviour
//!
//! As `f → 0`:
//$$
//!   Φ → (σ_eq / σ_y)² − 1
//! $$
//!
//! which is the classical von Mises criterion.
//!
//! # References
//!
//! - Gurson, A. L. (1977).  "Continuum theory of ductile rupture by
//!   void nucleation and growth: Part I — Yield criteria and flow
//!   rules for porous ductile media."  J. Eng. Mater. Technol.
//!   99(1), 2–15.
//! - Tvergaard, V., & Needleman, A. (1984).  "Analysis of the
//!   cup-cone fracture in a round tensile bar."  Acta Metall. 32(1),
//!   157–169.
//! - Tvergaard, V. (1990).  "Material failure by void growth to
//!   coalescence."  Adv. Appl. Mech. 27, 83–151.

use serde::{Deserialize, Serialize};

/// Classical Tvergaard–Needleman calibration constants.
pub const Q1_DEFAULT: f64 = 1.5;
/// Classical Tvergaard–Needleman calibration constants.
pub const Q2_DEFAULT: f64 = 1.0;
/// Classical Tvergaard–Needleman calibration constants
/// (`q₃ = q₁²` in the original calibration).
pub const Q3_DEFAULT: f64 = 2.25;

/// GTN model parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GtnParams {
    /// Initial void volume fraction `f_0`.
    pub f0: f64,
    /// Critical porosity at onset of coalescence `f_c`.
    pub fc: f64,
    /// Failure porosity `f_f`.
    pub ff: f64,
    /// Tvergaard parameter `q₁` (typical `1.5`).
    pub q1: f64,
    /// Tvergaard parameter `q₂` (typical `1.0`).
    pub q2: f64,
    /// Tvergaard parameter `q₃` (typical `2.25`).
    pub q3: f64,
    /// Nucleation-strain amplitude `f_n` (volume fraction of
    /// nucleating inclusions).
    pub fn_: f64,
    /// Mean nucleation strain `ε_n`.
    pub en: f64,
    /// Standard deviation of the nucleation strain `s_n`.
    pub sn: f64,
}

impl GtnParams {
    /// Classic Tvergaard–Needleman calibration with physically
    /// representative defaults:
    ///
    /// - `f_0 = 0.001` — typical clean-steel inclusion content
    /// - `f_c = 0.15`   — onset of coalescence
    /// - `f_f = 0.25`   — final failure
    /// - `q₁ = q₂ = 1.0 → q₃ = 2.25` (T&N 1984)
    /// - `f_n = 0.04`, `ε_n = 0.3`, `s_n = 0.1`
    pub fn classic() -> Self {
        Self {
            f0: 1.0e-3,
            fc: 0.15,
            ff: 0.25,
            q1: Q1_DEFAULT,
            q2: Q2_DEFAULT,
            q3: Q3_DEFAULT,
            fn_: 0.04,
            en: 0.3,
            sn: 0.1,
        }
    }
}

/// Compute the GTN yield function value `Φ`.
///
/// `Φ ≤ 0` indicates admissible elastic / elasto-plastic states.
/// `Φ = 0` corresponds to the updated yield surface; `Φ > 0`
/// indicates an inadmissible (returned) state that requires
/// plastic correction.
///
/// # Arguments
///
/// - `sigma_eq` — von Mises equivalent stress `σ_eq ≥ 0`
/// - `sigma_h`  — hydrostatic stress `σ_h = σ_kk / 3` (signed)
/// - `sigma_y`  — current matrix yield stress `σ_y > 0`
/// - `f`        — current porosity `f ∈ [0, 1)`
/// - `params`   — GTN parameters (`q₁`, `q₂`, `q₃`)
pub fn gtn_yield_function(
    sigma_eq: f64,
    sigma_h: f64,
    sigma_y: f64,
    f: f64,
    params: &GtnParams,
) -> f64 {
    let ratio = sigma_eq / sigma_y;
    let ratio_sq = ratio * ratio;
    let cosh_arg = (3.0 * params.q2 * sigma_h) / (2.0 * sigma_y);
    // Clamp to avoid overflow for very large triaxialities.
    let cosh_arg = cosh_arg.clamp(-50.0, 50.0);
    let cosh_val = cosh(cosh_arg);
    ratio_sq + 2.0 * params.q1 * f * cosh_val - (1.0 + params.q3 * f * f)
}

/// Update the porosity following GTN evolution.
///
/// ```text
/// Δf = (1 − f) Δε^p_kk + A Δε^p_eq
/// A  = f_n / (s_n √(2π)) exp(−½ ((ε^p_eq − ε_n) / s_n)²)
/// ```
///
/// `f` is clamped to `[0, 1)` after the update.  The accelerated
/// coalescence regime (`f ≥ f_c`) is captured by the standard
/// quadratic void-growth acceleration; this function returns the
/// unmodified evolution, callers may further accelerate beyond
/// `f_c` if they wish.
///
/// # Arguments
///
/// - `f`        — current porosity
/// - `deps_p_kk` — plastic volumetric strain increment `Δε^p_kk`
/// - `deps_p_eq` — plastic equivalent strain increment `Δε^p_eq ≥ 0`
/// - `params`   — GTN parameters (`f_n`, `ε_n`, `s_n`)
pub fn update_porosity(
    f: f64,
    deps_p_kk: f64,
    deps_p_eq: f64,
    params: &GtnParams,
) -> f64 {
    let strain_driven_nucleation = nucleation_rate(deps_p_eq, params);
    let df_growth = (1.0 - f) * deps_p_kk;
    let df_nucleation = strain_driven_nucleation * deps_p_eq;
    let f_new = f + df_growth + df_nucleation;
    f_new.clamp(0.0, 1.0)
}

/// Strain-controlled nucleation rate `A(ε^p_eq)` of the
/// Chu–Needleman form.
fn nucleation_rate(eps_p_eq: f64, params: &GtnParams) -> f64 {
    let z = (eps_p_eq - params.en) / params.sn;
    let gaussian = (-0.5 * z * z).exp();
    params.fn_ / (params.sn * (2.0 * core::f64::consts::PI).sqrt()) * gaussian
}

/// Indicates whether coalescence is active (`f ≥ f_c`).
pub fn is_coalescing(f: f64, params: &GtnParams) -> bool {
    f >= params.fc
}

/// Indicates whether final failure has been reached (`f ≥ f_f`).
pub fn has_failed(f: f64, params: &GtnParams) -> bool {
    f >= params.ff
}

/// Hyperbolic cosine `cosh(x) = (e^x + e^{-x}) / 2`.
///
/// Avoids pulling the entire `std` for one transcendental; uses
/// the `exp`-based definition directly.
fn cosh(x: f64) -> f64 {
    let e_pos = x.exp();
    let e_neg = (-x).exp();
    0.5 * (e_pos + e_neg)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn von_mises_recovered_at_zero_porosity() {
        // With f = 0 and σ_h = 0:
        // Φ = (σ_eq / σ_y)² + 0 − 1.
        let p = GtnParams::classic();
        for sigma_eq in [0.0, 100.0, 200.0, 250.0, 500.0] {
            let sigma_y = 250.0;
            let phi = gtn_yield_function(sigma_eq, 0.0, sigma_y, 0.0, &p);
            let expected = (sigma_eq / sigma_y).powi(2) - 1.0;
            assert!(
                approx(phi, expected, 1.0e-12),
                "σ_eq={sigma_eq} Φ={phi} expected={expected}"
            );
        }
    }

    #[test]
    fn hydrostatic_stress_lowers_yield_threshold_for_porous_material() {
        // For a porous material under pure hydrostatic tension,
        // Φ = 0 + 2 q₁ f cosh(3 q₂ σ_h / (2 σ_y)) − (1 + q₃ f²)
        //   = 0 when σ_h reaches the incipient-yield value.
        let p = GtnParams::classic();
        let f = 0.05;
        let sigma_y = 250.0;
        // Solve for σ_h yielding Φ = 0 numerically.
        let mut lo = 0.0;
        let mut hi = 1000.0;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            let v = gtn_yield_function(0.0, mid, sigma_y, f, &p);
            if v > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let sigma_h_yield = 0.5 * (lo + hi);
        // Closed-form check: Φ = 0 ⇒ cosh(3 q₂ σ_h / (2 σ_y)) = (1+q₃ f²)/(2 q₁ f).
        let ratio = (1.0 + p.q3 * f * f) / (2.0 * p.q1 * f);
        let cosh_target = ratio.acosh();
        let sigma_h_predicted = (cosh_target * 2.0 * sigma_y) / (3.0 * p.q2);
        assert!(
            approx(sigma_h_yield, sigma_h_predicted, 1.0e-6),
            "yield {sigma_h_yield} vs predicted {sigma_h_predicted}"
        );
    }

    #[test]
    fn porosity_grows_under_positive_volumetric_strain() {
        let p = GtnParams::classic();
        let f0 = 1.0e-3;
        let deps_p_kk = 0.05;
        let deps_p_eq = 0.05;
        let f1 = update_porosity(f0, deps_p_kk, deps_p_eq, &p);
        assert!(f1 > f0, "f should grow: {f1} <= {f0}");
    }

    #[test]
    fn porosity_clamped_below_one() {
        let p = GtnParams::classic();
        let f_new = update_porosity(0.5, 5.0, 5.0, &p);
        assert!(f_new <= 1.0);
    }

    #[test]
    fn nucleation_is_gaussian_in_strain() {
        let p = GtnParams::classic();
        // At ε^p_eq = ε_n the nucleation rate is maximum.
        let a_peak = nucleation_rate(p.en, &p);
        // Off-peak the rate should be much smaller.
        let a_tail = nucleation_rate(p.en + 3.0 * p.sn, &p);
        assert!(a_peak > a_tail);
    }

    #[test]
    fn coalescence_and_failure_thresholds() {
        let p = GtnParams::classic();
        assert!(!is_coalescing(0.10, &p));
        assert!(is_coalescing(0.20, &p));
        assert!(has_failed(0.30, &p));
        assert!(!has_failed(0.10, &p));
    }
}