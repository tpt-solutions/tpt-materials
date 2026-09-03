//! Residual-stress field from thermal-contraction eigenstrain.
//!
//! During cooling from the solidus, the substrate constrains the
//! contraction of the deposit.  The elastic eigenstrain is
//!
//! ```text
//! ε* = α · (T_solidus − T_room)
//! ```
//!
//! and the resulting stress in a 1-D bar with substrate constraint is
//!
//! ```text
//! σ = E · ε* · (1 − f) / (1 − ν)      (plane strain)
//! ```
//!
//! with `f` the fraction of strain relaxed by plastic flow.

use serde::{Deserialize, Serialize};

/// Thermal-contraction descriptor for a material.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThermalContraction {
    /// Coefficient of thermal expansion (1/K).
    pub cte: f64,
    /// Solidus temperature (K).
    pub solidus: f64,
    /// Final (room) temperature (K).
    pub room_temperature: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poissons_ratio: f64,
    /// Fraction of strain relaxed by plastic flow `f` (0–1).
    pub plastic_relaxation: f64,
}

impl Default for ThermalContraction {
    fn default() -> Self {
        Self {
            cte: 12.0e-6,
            solidus: 1878.0, // Ti-6Al-4V
            room_temperature: 298.0,
            youngs_modulus: 110.0e9,
            poissons_ratio: 0.32,
            plastic_relaxation: 0.7,
        }
    }
}

/// Residual-stress field sampled at a grid of positions (1-D along
/// build direction).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResidualStressField {
    /// Position along the build (m).
    pub positions: Vec<f64>,
    /// Residual stress at each position (Pa); tensile is positive.
    pub stresses: Vec<f64>,
}

/// Compute the residual stress from a single thermal-contraction
/// descriptor (point estimate).
pub fn residual_stress(tc: &ThermalContraction) -> f64 {
    let eigen = tc.cte * (tc.solidus - tc.room_temperature);
    let factor = (1.0 - tc.plastic_relaxation) / (1.0 - tc.poissons_ratio);
    tc.youngs_modulus * eigen * factor
}

/// Build a 1-D residual-stress field over `n` samples in
/// `[0, height]`.
pub fn residual_stress_field(
    tc: &ThermalContraction,
    height: f64,
    n: usize,
) -> ResidualStressField {
    let mut positions = Vec::with_capacity(n);
    let mut stresses = Vec::with_capacity(n);
    let sigma_0 = residual_stress(tc);
    for i in 0..n {
        let x = if n > 1 {
            height * (i as f64) / (n as f64 - 1.0)
        } else {
            0.0
        };
        // Linear stress profile (max at substrate, zero at top).
        let s = sigma_0 * (1.0 - x / height.max(1.0e-12));
        positions.push(x);
        stresses.push(s);
    }
    ResidualStressField {
        positions,
        stresses,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_stress_positive() {
        let tc = ThermalContraction::default();
        let s = residual_stress(&tc);
        assert!(s > 0.0);
    }

    #[test]
    fn residual_stress_zero_at_top_of_field() {
        let tc = ThermalContraction::default();
        let f = residual_stress_field(&tc, 1.0e-2, 11);
        let last = *f.stresses.last().unwrap();
        assert!(last.abs() < 1.0e-3);
    }

    #[test]
    fn residual_stress_max_at_substrate() {
        let tc = ThermalContraction::default();
        let f = residual_stress_field(&tc, 1.0e-2, 11);
        let first = f.stresses[0];
        for &s in &f.stresses {
            assert!(first >= s);
        }
    }

    #[test]
    fn plastic_relaxation_lowers_stress() {
        let mut tc = ThermalContraction::default();
        let s1 = residual_stress(&tc);
        tc.plastic_relaxation = 0.95;
        let s2 = residual_stress(&tc);
        assert!(s2 < s1);
    }
}
