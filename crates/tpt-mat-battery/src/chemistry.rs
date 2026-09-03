//! Battery chemistry descriptors and defaults.
//!
//! The diffusion coefficient `D` and partial molar volume `Ω` for a
//! range of common Li-ion chemistries are bundled here.  The values
//! are textbook ranges (Doyle, Fuller, Newman 1993; Zhang et al.
//! 2014); the caller should validate against their own
//! characterisation for predictive work.

use serde::{Deserialize, Serialize};

use crate::active_material::ActiveMaterial;

/// Battery chemistry (positive-electrode active material).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BatteryChemistry {
    /// NMC 811 (`LiNi₀.₈Mn₀.₁Co₀.₁O₂`).
    Nmc811,
    /// NMC 622 (`LiNi₀.₆Mn₀.₂Co₀.₂O₂`).
    Nmc622,
    /// LFP (`LiFePO₄`).
    Lfp,
    /// NCA (`LiNi₀.₈Co₀.₁₅Al₀.₀₅O₂`).
    Nca,
    /// Graphite negative electrode.
    Graphite,
    /// Silicon negative electrode.
    Silicon,
}

/// Return a sensible-default [`ActiveMaterial`] for a chemistry.
pub fn default_active_material(chemistry: BatteryChemistry) -> ActiveMaterial {
    match chemistry {
        BatteryChemistry::Nmc811 => ActiveMaterial {
            label: "NMC811".to_string(),
            chemistry,
            particle_radius: 5.0e-6,
            diffusion_coefficient: 1.0e-14,
            activation_energy_diffusion: 30_000.0,
            partial_molar_volume: 1.0e-5,
            youngs_modulus: 200.0e9,
            poissons_ratio: 0.25,
            fracture_stress: 700.0e6,
            max_concentration: 50_000.0,
        },
        BatteryChemistry::Nmc622 => ActiveMaterial {
            label: "NMC622".to_string(),
            chemistry,
            particle_radius: 5.0e-6,
            diffusion_coefficient: 2.0e-14,
            activation_energy_diffusion: 28_000.0,
            partial_molar_volume: 1.0e-5,
            youngs_modulus: 200.0e9,
            poissons_ratio: 0.25,
            fracture_stress: 700.0e6,
            max_concentration: 49_000.0,
        },
        BatteryChemistry::Lfp => ActiveMaterial {
            label: "LFP".to_string(),
            chemistry,
            particle_radius: 1.0e-7,
            diffusion_coefficient: 1.0e-13,
            activation_energy_diffusion: 35_000.0,
            partial_molar_volume: 3.0e-5,
            youngs_modulus: 120.0e9,
            poissons_ratio: 0.28,
            fracture_stress: 500.0e6,
            max_concentration: 23_000.0,
        },
        BatteryChemistry::Nca => ActiveMaterial {
            label: "NCA".to_string(),
            chemistry,
            particle_radius: 5.0e-6,
            diffusion_coefficient: 1.5e-14,
            activation_energy_diffusion: 30_000.0,
            partial_molar_volume: 1.0e-5,
            youngs_modulus: 200.0e9,
            poissons_ratio: 0.25,
            fracture_stress: 700.0e6,
            max_concentration: 50_000.0,
        },
        BatteryChemistry::Graphite => ActiveMaterial {
            label: "Graphite".to_string(),
            chemistry,
            particle_radius: 5.0e-6,
            diffusion_coefficient: 3.0e-14,
            activation_energy_diffusion: 35_000.0,
            partial_molar_volume: 4.0e-6,
            youngs_modulus: 10.0e9,
            poissons_ratio: 0.30,
            fracture_stress: 100.0e6,
            max_concentration: 30_000.0,
        },
        BatteryChemistry::Silicon => ActiveMaterial {
            label: "Silicon".to_string(),
            chemistry,
            particle_radius: 50.0e-9,
            diffusion_coefficient: 1.0e-16,
            activation_energy_diffusion: 40_000.0,
            partial_molar_volume: 9.0e-5,
            youngs_modulus: 90.0e9,
            poissons_ratio: 0.28,
            fracture_stress: 1.0e9,
            max_concentration: 220_000.0,
        },
    }
}
