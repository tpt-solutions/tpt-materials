//! Precipitation-strengthening increments: Orowan bypass and
//! particle shearing.
//!
//! These are returned as `StrengtheningIncrement`s ready for
//! hand-off to the `tpt_mat_hardening` crate as a Hall-Petch
//! style additive yield-strength bump.

use serde::{Deserialize, Serialize};

/// Strengthening mechanism for a precipitate population.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StrengtheningMechanism {
    /// Dislocations bow around impenetrable particles (Orowan).
    OrowanBypass,
    /// Dislocations shear through coherent / ordered particles.
    Shearing,
}

/// A strengthening increment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StrengtheningIncrement {
    /// Strengthening mechanism that produced it.
    pub mechanism: StrengtheningMechanism,
    /// Critical resolved shear-stress increment `Δτ` (Pa).
    pub delta_tau: f64,
}

/// Orowan bypass strengthening
/// `Δτ = 0.4 μ b / (π √(1 − ν) L)`
/// where `L = r (√(2π / (3 f)) − π / 2)` is the inter-particle
/// spacing for a 2-D random array of obstacles.
///
/// - `shear_modulus` `μ` (Pa)
/// - `burgers_vector` `b` (m)
/// - `poissons_ratio` `ν`
/// - `precipitate_radius` `r` (m)
/// - `volume_fraction` `f`
pub fn orowan_bypass_strengthening(
    shear_modulus: f64,
    burgers_vector: f64,
    poissons_ratio: f64,
    precipitate_radius: f64,
    volume_fraction: f64,
) -> f64 {
    if volume_fraction <= 0.0 || precipitate_radius <= 0.0 {
        return 0.0;
    }
    let denom = 1.0 - poissons_ratio;
    if denom <= 0.0 {
        return 0.0;
    }
    let one_minus_nu = denom.sqrt();
    let spacing = precipitate_radius * ((2.0 * core::f64::consts::PI / (3.0 * volume_fraction)).sqrt()
        - core::f64::consts::PI / 2.0);
    if spacing <= 0.0 {
        return 0.0;
    }
    0.4 * shear_modulus * burgers_vector / (core::f64::consts::PI * one_minus_nu * spacing)
}

/// Shear-strengthening increment (simplified Friedel form)
/// `Δτ = 0.13 μ b / L` for a sheared coherent particle.
pub fn shearing_strengthening(
    shear_modulus: f64,
    burgers_vector: f64,
    precipitate_radius: f64,
    volume_fraction: f64,
) -> f64 {
    if volume_fraction <= 0.0 || precipitate_radius <= 0.0 {
        return 0.0;
    }
    let spacing = precipitate_radius * ((2.0 * core::f64::consts::PI / (3.0 * volume_fraction)).sqrt()
        - core::f64::consts::PI / 2.0);
    if spacing <= 0.0 {
        return 0.0;
    }
    0.13 * shear_modulus * burgers_vector / spacing
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn orowan_zero_at_zero_volume_fraction() {
        let dt = orowan_bypass_strengthening(80.0e9, 2.5e-10, 0.3, 5.0e-9, 0.0);
        assert!(approx(dt, 0.0, 1.0e-12));
    }

    #[test]
    fn orowan_grows_with_volume_fraction() {
        let dt1 = orowan_bypass_strengthening(80.0e9, 2.5e-10, 0.3, 5.0e-9, 0.01);
        let dt2 = orowan_bypass_strengthening(80.0e9, 2.5e-10, 0.3, 5.0e-9, 0.05);
        assert!(dt2 > dt1);
    }

    #[test]
    fn shearing_zero_at_zero_volume_fraction() {
        let dt = shearing_strengthening(80.0e9, 2.5e-10, 5.0e-9, 0.0);
        assert!(approx(dt, 0.0, 1.0e-12));
    }
}