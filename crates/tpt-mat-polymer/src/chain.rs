//! Chain-model definitions.

use serde::{Deserialize, Serialize};

/// Single-chain constitutive model selection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ChainModel {
    /// Arruda–Boyce 8-chain model.
    ArrudaBoyce {
        /// Number of segments per chain `N`.
        n_segments: u32,
        /// Shear modulus of the crosslinked network (Pa).
        shear_modulus: f64,
    },
    /// Worm-like chain (Marko–Siggia interpolation).
    WormLikeChain {
        /// Persistence length `l_p` (m).
        persistence_length: f64,
        /// Contour length `L_c` (m).
        contour_length: f64,
        /// Chain number density per unit volume (m⁻³).
        k_b: f64,
    },
    /// Freely jointed chain (inverse-Langevin exact).
    FreelyJointedChain {
        /// Segment length `b` (m).
        segment_length: f64,
        /// Number of segments `N`.
        num_segments: u32,
        /// Chain number density per unit volume (m⁻³).
        k_b: f64,
    },
}

/// A polymer material model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolymerModel {
    /// Chain model.
    pub chain: ChainModel,
    /// Crosslink density (mol/m³) — currently only used as a
    /// metadata field; future extensions may use it to scale the
    /// network shear modulus.
    pub crosslink_density: f64,
}
