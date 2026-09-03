//! Per-slip-system dislocation-density state.

use serde::{Deserialize, Serialize};

/// Statistically-stored dislocation (SSD) and geometrically-
/// necessary dislocation (GND) density state for a single slip
/// system, plus a scalar "forest" density that the storage
/// term `k₁ √ρ_f` uses.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DislocationDensityState {
    /// SSD density (m⁻²).
    pub rho_ssd: f64,
    /// GND density (m⁻²), `≥ 0`.
    pub rho_gnd: f64,
    /// Forest dislocation density (m⁻²), `≥ 0`.
    pub rho_forest: f64,
}

impl Default for DislocationDensityState {
    fn default() -> Self {
        Self {
            rho_ssd: 1.0e10,
            rho_gnd: 0.0,
            rho_forest: 1.0e10,
        }
    }
}

impl DislocationDensityState {
    /// Total dislocation density `ρ_tot = ρ_SSD + ρ_GND`.
    pub fn total(&self) -> f64 {
        self.rho_ssd + self.rho_gnd
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_density_is_sum() {
        let s = DislocationDensityState {
            rho_ssd: 1.0e14,
            rho_gnd: 5.0e13,
            rho_forest: 1.5e14,
        };
        assert!((s.total() - 1.5e14).abs() < 1.0e-6);
    }
}
