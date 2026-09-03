//! Anisotropy model for the interface energy / mobility.

use serde::{Deserialize, Serialize};

/// Anisotropy mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum AnisotropyMode {
    /// 4-fold cubic: `a(θ) = α (1 + ε cos(4 θ))`.
    Cubic4Fold,
    /// 6-fold hexagonal: `a(θ) = α (1 + ε cos(6 θ))`.
    Hexagonal6Fold,
    /// Isotropic (no angular dependence).
    Isotropic,
}

/// Anisotropy parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AnisotropyModel {
    /// Strength `α ∈ [0, 1]` (fraction of the interface energy that is anisotropic).
    pub strength: f64,
    /// Anisotropy mode.
    pub mode: AnisotropyMode,
}

impl Default for AnisotropyModel {
    fn default() -> Self {
        Self {
            strength: 0.03,
            mode: AnisotropyMode::Cubic4Fold,
        }
    }
}

impl AnisotropyModel {
    /// Evaluate the anisotropy multiplier `a(θ)` where `θ` is the
    /// angle between the interface normal and the crystal x-axis.
    pub fn evaluate(&self, theta: f64) -> f64 {
        match self.mode {
            AnisotropyMode::Isotropic => 1.0,
            AnisotropyMode::Cubic4Fold => 1.0 + self.strength * (4.0 * theta).cos(),
            AnisotropyMode::Hexagonal6Fold => 1.0 + self.strength * (6.0 * theta).cos(),
        }
    }

    /// `da/dθ` (derivative of the anisotropy multiplier).
    pub fn derivative(&self, theta: f64) -> f64 {
        match self.mode {
            AnisotropyMode::Isotropic => 0.0,
            AnisotropyMode::Cubic4Fold => -4.0 * self.strength * (4.0 * theta).sin(),
            AnisotropyMode::Hexagonal6Fold => -6.0 * self.strength * (6.0 * theta).sin(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotropic_is_unity() {
        let a = AnisotropyModel {
            strength: 0.5,
            mode: AnisotropyMode::Isotropic,
        };
        assert_eq!(a.evaluate(0.5), 1.0);
        assert_eq!(a.derivative(0.5), 0.0);
    }

    #[test]
    fn cubic_has_maxima_at_zero_and_quarter_period() {
        let a = AnisotropyModel {
            strength: 0.05,
            mode: AnisotropyMode::Cubic4Fold,
        };
        // a(0) = 1 + ε; minima at θ = π/4.
        assert!(a.evaluate(0.0) > 1.0);
        assert!(a.evaluate(std::f64::consts::FRAC_PI_4) < 1.0);
    }
}
