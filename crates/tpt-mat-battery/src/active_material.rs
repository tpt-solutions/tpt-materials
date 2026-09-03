//! Battery active-material descriptors.

use serde::{Deserialize, Serialize};

use crate::chemistry::BatteryChemistry;

/// A battery active material (positive or negative electrode).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ActiveMaterial {
    /// Human-readable label.
    pub label: String,
    /// Chemistry identifier.
    pub chemistry: BatteryChemistry,
    /// Spherical particle radius (m).
    pub particle_radius: f64,
    /// Solid-state Li diffusion coefficient `D` at the reference
    /// temperature (m²/s).
    pub diffusion_coefficient: f64,
    /// Activation energy for diffusion `Q` (J/mol).
    pub activation_energy_diffusion: f64,
    /// Partial molar volume of Li `Ω` (m³/mol).
    pub partial_molar_volume: f64,
    /// Young's modulus `E` (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio `ν`.
    pub poissons_ratio: f64,
    /// Critical stress for particle cracking `σ_c` (Pa).
    pub fracture_stress: f64,
    /// Maximum Li concentration `c_max` (mol/m³).
    pub max_concentration: f64,
}

impl ActiveMaterial {
    /// Diffusion coefficient at temperature `T` (K).
    ///
    /// `D(T) = D_0 exp(−Q / (R T))`.
    pub fn diffusion_at(&self, temperature: f64) -> f64 {
        const R: f64 = 8.314_462_618;
        self.diffusion_coefficient * (-self.activation_energy_diffusion / (R * temperature)).exp()
    }
}
