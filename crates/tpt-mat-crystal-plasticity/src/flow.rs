//! Viscoplastic flow rule and plastic velocity gradient.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_math_linalg_fixed::{Mat3, Vec3, Vec6};
use tpt_mat_crystallography::{SlipSystem, SlipSystemError};

use crate::model::RateSensitivity;

/// A single increment of plastic flow: the velocity gradient `L^p` and
/// the slip rates `γ̇^α` that produced it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlasticIncrement {
    /// Plastic part of the velocity gradient `L^p = Σ_α γ̇^α s^α ⊗ n^α`.
    pub velocity_gradient: Mat3,
    /// Slip rate per slip system (1/s).
    pub slip_rates: Vec<f64>,
}

impl PlasticIncrement {
    /// Construct from raw data.
    pub fn new(velocity_gradient: Mat3, slip_rates: Vec<f64>) -> Self {
        Self {
            velocity_gradient,
            slip_rates,
        }
    }
}

/// Errors raised by [`power_law_slip_rate`].
#[derive(Debug, Error, PartialEq)]
pub enum FlowRuleError {
    /// CRSS was non-positive.
    #[error("critical resolved shear stress must be positive")]
    NonPositiveCrss,
    /// Reference slip rate was non-positive.
    #[error("reference slip rate must be positive")]
    NonPositiveReferenceRate,
    /// Rate-sensitivity exponent was non-positive.
    #[error("rate-sensitivity exponent must be positive")]
    NonPositiveExponent,
}

/// Power-law slip rate `γ̇^α = γ̇_0 |τ^α / τ_c^α|^n sign(τ^α)`.
pub fn power_law_slip_rate(
    rss: f64,
    crss: f64,
    rate: &RateSensitivity,
) -> Result<f64, FlowRuleError> {
    if crss <= 0.0 {
        return Err(FlowRuleError::NonPositiveCrss);
    }
    if rate.reference_strain_rate <= 0.0 {
        return Err(FlowRuleError::NonPositiveReferenceRate);
    }
    if rate.exponent <= 0.0 {
        return Err(FlowRuleError::NonPositiveExponent);
    }
    let ratio = (rss / crss).abs();
    let mag = rate.reference_strain_rate * ratio.powf(rate.exponent);
    Ok(mag * rss.signum())
}

/// Compute the plastic velocity gradient `L^p = Σ_α γ̇^α s^α ⊗ n^α`.
///
/// Slip systems are supplied as `(slip_direction, plane_normal)` pairs
/// in the sample frame; if they were originally given in the crystal
/// frame, the caller must rotate them with the lattice rotation first.
pub fn viscoplastic_velocity_gradient(
    slip_directions: &[Vec3],
    plane_normals: &[Vec3],
    slip_rates: &[f64],
) -> Result<Mat3, SlipSystemError> {
    let n = slip_directions.len();
    if plane_normals.len() != n || slip_rates.len() != n {
        return Err(SlipSystemError::ZeroDirection);
    }
    let mut acc = Mat3::ZERO;
    for ((&s, &n_), &g) in slip_directions.iter().zip(plane_normals).zip(slip_rates) {
        // s ⊗ n stored column-major: column j = s * n_j.
        let mut m = Mat3::ZERO;
        m.set(0, 0, s[0] * n_[0]);
        m.set(1, 0, s[1] * n_[0]);
        m.set(2, 0, s[2] * n_[0]);
        m.set(0, 1, s[0] * n_[1]);
        m.set(1, 1, s[1] * n_[1]);
        m.set(2, 1, s[2] * n_[1]);
        m.set(0, 2, s[0] * n_[2]);
        m.set(1, 2, s[1] * n_[2]);
        m.set(2, 2, s[2] * n_[2]);
        acc = acc + m.scale(g);
    }
    Ok(acc)
}

/// Compute slip-system resolved shear stress `τ^α = σ_ij (s_i n_j)`
/// for every slip system.  Slip directions / normals are in the same
/// frame as the supplied stress (tensor convention; the supplied
/// [`Vec6`] is the Cauchy stress tensor in Voigt order without the
/// engineering-shear factor of 2 for the off-diagonal components).
pub fn resolved_shear_stresses(
    sigma: Vec6,
    slip_systems: &[SlipSystem],
) -> Vec<f64> {
    let s = sigma.data;
    slip_systems
        .iter()
        .map(|slip| {
            let sd = slip.slip_direction;
            let nd = slip.plane_normal;
            // Direct tensor contraction σ_ij · s_i · n_j.
            //   τ = σ_xx s_x n_x + σ_yy s_y n_y + σ_zz s_z n_z
            //     + σ_xy (s_x n_y + s_y n_x) + σ_yz (s_y n_z + s_z n_y)
            //     + σ_xz (s_x n_z + s_z n_x).
            s[0] * sd[0] * nd[0]
                + s[1] * sd[1] * nd[1]
                + s[2] * sd[2] * nd[2]
                + s[3] * (sd[0] * nd[1] + sd[1] * nd[0])
                + s[4] * (sd[1] * nd[2] + sd[2] * nd[1])
                + s[5] * (sd[0] * nd[2] + sd[2] * nd[0])
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use tpt_mat_crystallography::CrystalStructure;

    #[test]
    fn power_law_zero_rss_yields_zero_rate() {
        let r = power_law_slip_rate(0.0, 10.0, &RateSensitivity::default()).unwrap();
        assert_eq!(r, 0.0);
    }

    #[test]
    fn power_law_sign_matches_rss() {
        let rate = RateSensitivity::default();
        let pos = power_law_slip_rate(15.0, 10.0, &rate).unwrap();
        let neg = power_law_slip_rate(-15.0, 10.0, &rate).unwrap();
        assert!(pos > 0.0 && neg < 0.0);
        assert!((pos + neg).abs() < 1e-12);
    }

    #[test]
    fn velocity_gradient_is_trace_free_for_pure_slip() {
        // For FCC {111}<110>, s ⊥ n ⇒ s ⊗ n is traceless ⇒ L^p traceless.
        let slips = CrystalStructure::FCC.slip_systems();
        let s_dirs: Vec<_> = slips.iter().map(|s| s.slip_direction).collect();
        let n_dirs: Vec<_> = slips.iter().map(|s| s.plane_normal).collect();
        let rates = vec![1.0; slips.len()];
        let l = viscoplastic_velocity_gradient(&s_dirs, &n_dirs, &rates).unwrap();
        assert_relative_eq!(l.trace(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn rss_matches_definition_for_simple_stress() {
        // For shear stress σ_xz on the FCC system (s=[1,-1,0]/√2,
        // n=[1,1,-1]/√3):
        //   τ = σ_ij (s_i n_j) = σ_xz (s_x n_z + s_z n_x).
        let s = Vec3::new(1.0, -1.0, 0.0).normalized();
        let n = Vec3::new(1.0, 1.0, -1.0).normalized();
        let sigma = Vec6::new(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let expected = sigma[5] * (s[0] * n[2] + s[2] * n[0]);
        let slip = tpt_mat_crystallography::SlipSystem::from_vectors(
            s,
            n,
            1.0,
            tpt_mat_crystallography::SlipFamily::Fcc110,
        )
        .unwrap();
        let rss = resolved_shear_stresses(sigma, &[slip]);
        assert!((rss[0] - expected).abs() < 1e-9, "got {} expected {}", rss[0], expected);
    }
}