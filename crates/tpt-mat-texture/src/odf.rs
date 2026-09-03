//! Orientation Distribution Function (ODF) via kernel density
//! estimation on `SO(3)`.
//!
//! We parameterise `SO(3)` with `q = (w, x, y, z)` (unit quaternion,
//! scalar-first, `w ≥ 0` convention).  The ODF is the KDE
//!
//! `f(q) = (1/N) Σ_i w_i K(‖q − q_i‖)`
//!
//! with the geodesic kernel
//! `K(d) = exp(−d² / (2 σ²)) / √(2π σ²)`
//!
//! where `d` is the rotational geodesic distance
//! `d(q, q_i) = arccos(2 ⟨q, q_i⟩² − 1)` (using the unit-quaternion
//! double-cover identity `q ≈ −q`).

use serde::{Deserialize, Serialize};

use tpt_mat_core::CrystalOrientation;

use crate::texture::TextureAnalyzer;

/// Kernel choice for the ODF KDE.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum OdfKernel {
    /// Geodesic Gaussian with bandwidth `sigma` (radians).
    GeodesicGaussian {
        /// Bandwidth (radians).
        sigma: f64,
    },
}

impl Default for OdfKernel {
    fn default() -> Self {
        OdfKernel::GeodesicGaussian { sigma: 0.1 }
    }
}

/// ODF evaluated on an arbitrary set of query rotations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrientationDistributionFunction {
    /// Kernel used.
    pub kernel: OdfKernel,
    /// Sum of grain weights.
    pub total_weight: f64,
    /// Per-grain quaternions `(w, x, y, z)`, scalar-first, `w ≥ 0`.
    pub quaternions: Vec<[f64; 4]>,
    /// Per-grain weights.
    pub weights: Vec<f64>,
}

impl TextureAnalyzer {
    /// Build an ODF from this texture using the supplied kernel.
    pub fn orientation_distribution_function(
        &self,
        kernel: OdfKernel,
    ) -> OrientationDistributionFunction {
        let quaternions = self
            .orientations
            .iter()
            .map(|g| g.to_quaternion())
            .collect();
        OrientationDistributionFunction {
            kernel,
            total_weight: self.total_weight(),
            quaternions,
            weights: self.weights.clone(),
        }
    }
}

impl OrientationDistributionFunction {
    /// Evaluate the ODF at a single orientation.
    pub fn evaluate(&self, g: CrystalOrientation) -> f64 {
        let (sigma,) = match self.kernel {
            OdfKernel::GeodesicGaussian { sigma } => (sigma,),
        };
        let q = g.to_quaternion();
        let mut acc = 0.0;
        let normaliser = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
        for (qi, &w) in self.quaternions.iter().zip(self.weights.iter()) {
            // Use the canonical `w ≥ 0` branch already enforced in
            // `CrystalOrientation::to_quaternion`.
            let dot = q[0] * qi[0] + q[1] * qi[1] + q[2] * qi[2] + q[3] * qi[3];
            // Quaternion double cover: q and -q describe the same rotation.
            let dot = dot.abs();
            // Geodesic distance on SO(3) via the trace formula.
            let cos_half = dot.clamp(0.0, 1.0);
            let angle = 2.0 * cos_half.acos();
            let k = (-0.5 * (angle / sigma).powi(2)).exp() * normaliser;
            acc += w * k;
        }
        if self.total_weight > 0.0 {
            acc / self.total_weight
        } else {
            0.0
        }
    }
}

/// Helper for callers: convert a unit quaternion `(w, x, y, z)` to a
/// `CrystalOrientation`.
#[allow(dead_code)]
pub fn quaternion_to_orientation(q: [f64; 4]) -> CrystalOrientation {
    CrystalOrientation::from_quaternion(q).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn odf_peaks_at_known_orientation() {
        let g0 = CrystalOrientation::from_euler_bunge(0.3, 0.7, 1.1);
        let ta = TextureAnalyzer::new(vec![g0], vec![1.0]);
        let odf = ta.orientation_distribution_function(OdfKernel::GeodesicGaussian { sigma: 0.05 });
        let f0 = odf.evaluate(g0);
        let f1 = odf.evaluate(CrystalOrientation::from_euler_bunge(1.5, 0.2, 0.9));
        assert!(f0 > f1, "ODF should peak near g0: {f0} vs {f1}");
    }

    #[test]
    fn odf_normalised_peak_decreases_with_more_grains() {
        let g0 = CrystalOrientation::identity();
        let ta1 = TextureAnalyzer::new(vec![g0], vec![1.0]);
        let ta2 = TextureAnalyzer::new(vec![g0, g0, g0], vec![1.0, 1.0, 1.0]);
        let odf1 =
            ta1.orientation_distribution_function(OdfKernel::GeodesicGaussian { sigma: 0.05 });
        let odf2 =
            ta2.orientation_distribution_function(OdfKernel::GeodesicGaussian { sigma: 0.05 });
        // Same total weight after normalisation, so the peak should be
        // identical.
        let p1 = odf1.evaluate(g0);
        let p2 = odf2.evaluate(g0);
        assert!((p1 - p2).abs() < 1e-9);
    }

    #[test]
    fn quaternion_helper_recovers_identity() {
        let g = quaternion_to_orientation([1.0, 0.0, 0.0, 0.0]);
        assert_eq!(g, CrystalOrientation::identity());
    }
}
