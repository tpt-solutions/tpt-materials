//! Texture analyzer: a collection of grain orientations + weights.

use serde::{Deserialize, Serialize};

use tpt_mat_core::{CrystalOrientation, Grain, MaterialMicrostructure};

/// Collection of grain orientations with associated weights (typically
/// volume fractions).  Built from either an explicit list or from the
/// grains of a [`MaterialMicrostructure`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TextureAnalyzer {
    /// Orientations (sample → crystal frame).
    pub orientations: Vec<CrystalOrientation>,
    /// Per-grain weight (default 1.0 per grain, i.e. equal volume).
    pub weights: Vec<f64>,
}

impl TextureAnalyzer {
    /// Construct an explicit analyzer.
    pub fn new(orientations: Vec<CrystalOrientation>, weights: Vec<f64>) -> Self {
        debug_assert_eq!(orientations.len(), weights.len());
        Self {
            orientations,
            weights,
        }
    }

    /// Build an analyzer from the grains of a [`MaterialMicrostructure`].
    /// Each grain gets unit weight; callers who need volume-weighted
    /// statistics should call [`Self::with_weights`].
    pub fn from_ebsd(material: &MaterialMicrostructure) -> Self {
        let mut orientations = Vec::with_capacity(material.grains.len());
        for grain in material.grains.values() {
            orientations.push(grain.orientation);
        }
        let weights = vec![1.0; orientations.len()];
        Self::new(orientations, weights)
    }

    /// Replace the weights (e.g. with grain volume fractions).
    pub fn with_weights(mut self, weights: Vec<f64>) -> Self {
        debug_assert_eq!(self.orientations.len(), weights.len());
        self.weights = weights;
        self
    }

    /// Number of grains.
    pub fn len(&self) -> usize {
        self.orientations.len()
    }

    /// Whether the analyzer is empty.
    pub fn is_empty(&self) -> bool {
        self.orientations.is_empty()
    }

    /// Sum of the weights (for normalisation).
    pub fn total_weight(&self) -> f64 {
        self.weights.iter().sum()
    }
}

/// Convenience: collect orientations from a slice of grains.
#[allow(dead_code)]
pub fn orientations_from_grains(grains: &[Grain]) -> Vec<CrystalOrientation> {
    grains.iter().map(|g| g.orientation).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_mat_core::{Grain, GrainId, Phase, PhaseId};

    #[test]
    fn empty_has_zero_total_weight() {
        let ta = TextureAnalyzer::new(Vec::new(), Vec::new());
        assert_eq!(ta.total_weight(), 0.0);
        assert!(ta.is_empty());
    }

    #[test]
    fn from_ebsd_assigns_unit_weights() {
        let mut mat = MaterialMicrostructure::new("m", "test");
        let c = tpt_mat_core::Composition::new(
            tpt_mat_core::CompositionBasis::Atomic,
            std::iter::once(("Fe".to_string(), 1.0)).collect(),
            1e-9,
        )
        .unwrap();
        mat.add_phase(Phase::new(PhaseId(0), "alpha-Fe", c, 1.0));
        for i in 0..3 {
            mat.add_grain(Grain::new(
                GrainId(i),
                PhaseId(0),
                CrystalOrientation::identity(),
                tpt_math_linalg_fixed::Vec3::ZERO,
                1e-5,
            ));
        }
        let ta = TextureAnalyzer::from_ebsd(&mat);
        assert_eq!(ta.len(), 3);
        assert!((ta.total_weight() - 3.0).abs() < 1e-12);
    }
}