//! Phase-field model selection.

use serde::{Deserialize, Serialize};

/// Phase-field model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum PhaseFieldModel {
    /// Allen-Cahn (non-conserved order parameter, `η`):
    /// `∂η/∂t = −L (δF/δη) = −L (f_bulk′(η) − κ ∇²η)`.
    AllenCahn,
    /// Cahn-Hilliard (conserved concentration, `c`):
    /// `∂c/∂t = ∇ · (M ∇(μ))`, `μ = f_bulk′(c) − κ ∇²c`.
    CahnHilliard,
    /// Kobayashi (Allen-Cahn with latent-heat thermal coupling):
    /// `∂η/∂t = L (κ ∇²η − f_bulk′(η) + λ u)`,
    /// `∂u/∂t = D ∇²u + L/2 ∂η/∂t`,
    /// where `u` is dimensionless temperature.
    Kobayashi,
    /// Multi-phase Allen-Cahn: `num_grains` non-conserved order
    /// parameters coupled through `Σ η_i = 1` at every point.
    MultiPhase {
        /// Number of grains (one order parameter per grain).
        num_grains: usize,
    },
}
