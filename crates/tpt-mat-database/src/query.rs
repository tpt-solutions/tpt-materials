//! Property-range queries.

use serde::{Deserialize, Serialize};

use super::record::MaterialRecord;

/// Single named property that can be queried on a record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Property {
    /// Young's modulus (GPa).
    YoungsModulusGPa,
    /// Yield strength (MPa).
    YieldStrengthMPa,
    /// Ultimate tensile strength (MPa).
    UltimateTensileStrengthMPa,
    /// Density (kg/m³).
    DensityKgPerM3,
    /// Poisson's ratio.
    PoissonsRatio,
    /// Thermal conductivity (W/(m·K)).
    ThermalConductivityWPerMK,
    /// Specific heat (J/(kg·K)).
    SpecificHeatJPerKgK,
    /// Coefficient of thermal expansion (1/K).
    CoefficientThermalExpansionPerK,
    /// Melting point (K).
    MeltingPointK,
    /// Electrical resistivity (Ω·m).
    ElectricalResistivityOhmM,
}

/// Property-range query: inclusive `[min, max]`.  Either bound may
/// be `None` for an open-ended query.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PropertyQuery {
    /// Property to query.
    pub property: Property,
    /// Lower bound (inclusive); `None` for `-∞`.
    pub min: Option<f64>,
    /// Upper bound (inclusive); `None` for `+∞`.
    pub max: Option<f64>,
}

impl PropertyQuery {
    /// Construct a query for `property`.
    pub fn new(property: Property) -> Self {
        Self {
            property,
            min: None,
            max: None,
        }
    }

    /// Set the lower bound (inclusive).
    pub fn min(mut self, v: f64) -> Self {
        self.min = Some(v);
        self
    }

    /// Set the upper bound (inclusive).
    pub fn max(mut self, v: f64) -> Self {
        self.max = Some(v);
        self
    }

    /// `true` iff the record's value lies in the closed range.
    pub fn matches(&self, record: &MaterialRecord) -> bool {
        let v = self.extract(record);
        if let Some(lo) = self.min {
            if v < lo {
                return false;
            }
        }
        if let Some(hi) = self.max {
            if v > hi {
                return false;
            }
        }
        true
    }

    fn extract(&self, record: &MaterialRecord) -> f64 {
        match self.property {
            Property::YoungsModulusGPa => record.mechanical.youngs_modulus_gpa,
            Property::YieldStrengthMPa => record.mechanical.yield_strength_mpa,
            Property::UltimateTensileStrengthMPa => {
                record.mechanical.ultimate_tensile_strength_mpa
            }
            Property::DensityKgPerM3 => record.mechanical.density_kg_per_m3,
            Property::PoissonsRatio => record.mechanical.poissons_ratio,
            Property::ThermalConductivityWPerMK => {
                record.thermal.thermal_conductivity_w_per_mk
            }
            Property::SpecificHeatJPerKgK => record.thermal.specific_heat_j_per_kgk,
            Property::CoefficientThermalExpansionPerK => {
                record.thermal.coefficient_thermal_expansion_per_k
            }
            Property::MeltingPointK => record.thermal.melting_point_k,
            Property::ElectricalResistivityOhmM => {
                record.electrical.electrical_resistivity_ohm_m
            }
        }
    }
}