//! Grain-boundary mobility, including misorientation dependence.

use serde::{Deserialize, Serialize};

/// Mobility parameters for grain boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MobilityParams {
    /// Base mobility `M_0` (units of `mobility` in the solver).
    pub base_mobility: f64,
    /// Activation energy misorientation threshold (radians).  Below
    /// this angle the boundary is treated as a low-angle boundary
    /// (Read–Shockley, mobility ∝ θ).
    pub low_angle_threshold: f64,
    /// High-angle mobility multiplier (>1 for clean HAGBs).
    pub high_angle_multiplier: f64,
}

impl Default for MobilityParams {
    fn default() -> Self {
        Self {
            base_mobility: 1.0,
            low_angle_threshold: 10.0_f64.to_radians(),
            high_angle_multiplier: 5.0,
        }
    }
}

/// Misorientation-dependent grain-boundary mobility.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GrainBoundaryMobility {
    /// Mobility parameters.
    pub params: MobilityParams,
}

impl Default for GrainBoundaryMobility {
    fn default() -> Self {
        Self {
            params: MobilityParams::default(),
        }
    }
}

impl GrainBoundaryMobility {
    /// Mobility for a boundary with misorientation `theta` (radians).
    ///
    /// Low-angle boundaries (`θ < θ_c`) follow `M = M_0 · θ / θ_c`
    /// (Read–Shockley linearised mobility).  High-angle boundaries
    /// get a fixed `M = M_0 · high_angle_multiplier`.
    pub fn at_misorientation(&self, theta: f64) -> f64 {
        let p = self.params;
        if theta.abs() < p.low_angle_threshold {
            p.base_mobility * (theta.abs() / p.low_angle_threshold)
        } else {
            p.base_mobility * p.high_angle_multiplier
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_angle_is_linear() {
        let m = GrainBoundaryMobility::default();
        let m1 = m.at_misorientation(0.05_f64.to_radians());
        let m2 = m.at_misorientation(0.10_f64.to_radians());
        assert!((m2 - 2.0 * m1).abs() < 1e-12);
    }

    #[test]
    fn high_angle_saturates() {
        let m = GrainBoundaryMobility::default();
        let m_high = m.at_misorientation(0.5);
        let m_higher = m.at_misorientation(1.0);
        assert_eq!(m_high, m_higher);
    }
}
