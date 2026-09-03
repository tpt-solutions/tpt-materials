//! Grain-boundary diffusion (Fisher-type approximation).
//!
//! In the Fisher model the effective in-plane diffusivity through a
//! polycrystal with grain-boundary width `δ` and grain-boundary
//! diffusivity `D_gb` is `D_eff = D_l + (π δ / L) D_gb` (Suzuki et al.
//! form, Regime A/B transition not modelled here).  This is enough for
//! a sanity-check component; full Lea-Fox-Le Claire treatment is
//! deferred.

use serde::{Deserialize, Serialize};

/// Arrhenius parameters describing a lattice / grain-boundary
/// diffusivity pair.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GBDiffusivity {
    /// Lattice diffusivity `D_l`.
    pub d_lattice: f64,
    /// Grain-boundary diffusivity `D_gb`.
    pub d_grain_boundary: f64,
    /// Grain-boundary thickness `δ` (same units as grain size).
    pub gb_width: f64,
    /// Mean grain size `L`.
    pub grain_size: f64,
}

impl GBDiffusivity {
    /// `D_lattice` (lattice-only contribution).
    pub fn lattice(&self) -> f64 {
        self.d_lattice
    }

    /// Effective 1D diffusivity (Fisher / Harrison Regime A):
    /// `D_eff = D_l + (π δ / L) D_gb`.
    pub fn effective(&self) -> f64 {
        if self.grain_size <= 0.0 {
            return self.d_lattice;
        }
        self.d_lattice
            + std::f64::consts::PI * (self.gb_width / self.grain_size) * self.d_grain_boundary
    }
}

/// Helper bundling a grain-boundary diffusivity.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GrainBoundaryDiffusion {
    /// Parameters.
    pub params: GBDiffusivity,
}

impl Default for GrainBoundaryDiffusion {
    fn default() -> Self {
        Self {
            params: GBDiffusivity {
                d_lattice: 1.0e-12,
                d_grain_boundary: 1.0e-9,
                gb_width: 5.0e-10,
                grain_size: 1.0e-5,
            },
        }
    }
}

impl GrainBoundaryDiffusion {
    /// Convenience: returns the Fisher effective diffusivity.
    pub fn effective_diffusivity(&self) -> f64 {
        self.params.effective()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gb_enhances_lattice() {
        let g = GrainBoundaryDiffusion::default();
        assert!(g.effective_diffusivity() > g.params.d_lattice);
    }

    #[test]
    fn large_grains_recover_lattice() {
        // GB enhancement = π δ D_gb / L → 0 as L → ∞, so the
        // effective diffusivity tends to D_lattice.
        let g = GrainBoundaryDiffusion {
            params: GBDiffusivity {
                d_lattice: 1.0,
                d_grain_boundary: 100.0,
                gb_width: 0.01,
                grain_size: 1.0e3,
            },
            ..GrainBoundaryDiffusion::default()
        };
        let eff = g.effective_diffusivity();
        let enhancement = (eff - g.params.d_lattice) / g.params.d_lattice;
        assert!(
            enhancement < 1e-2,
            "expected small enhancement, got {enhancement}"
        );
    }
}
