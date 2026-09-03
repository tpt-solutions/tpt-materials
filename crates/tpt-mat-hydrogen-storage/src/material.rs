//! Hydrogen storage material types.

use serde::{Deserialize, Serialize};

/// Hydride type classification.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum HydrideType {
    /// Interstitial metal hydride (e.g. `LaNi_5H_6`).
    MetalHydride {
        /// Alloy designation (e.g. "LaNi5").
        alloy: String,
    },
    /// Chemical hydride (e.g. `NaBH_4`, `LiBH_4`).
    ChemicalHydride {
        /// Compound formula.
        compound: String,
    },
    /// Porous adsorbent (MOF, activated carbon, …).
    PorousMaterial(PorousAdsorbent),
}

/// Porous-adsorbent parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PorousAdsorbent {
    /// BET surface area (m²/g).
    pub surface_area_m2_per_g: f64,
    /// Pore volume (cm³/g).
    pub pore_volume_cm3_per_g: f64,
}

/// A hydrogen storage material.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HydrogenStorageMaterial {
    /// Hydride type.
    pub hydride_type: HydrideType,
    /// Reversible storage capacity in `H/M` (hydrogen atoms per
    /// metal atom) or in `wt%`, depending on the material.
    pub storage_capacity: f64,
    /// Enthalpy of desorption `ΔH_des` (J/mol H₂).
    pub desorption_enthalpy_j_per_mol: f64,
    /// Entropy of desorption `ΔS_des` (J/(mol H₂·K)).
    pub desorption_entropy_j_per_mol_k: f64,
}