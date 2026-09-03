//! Material record schema.

use serde::{Deserialize, Serialize};

use super::provenance::DataSource;

/// Composition expressed in mass fractions.  Sums do not have to
/// be `1.0` — fields default to `0.0` and `*_balance` is the
/// remainder (commonly Fe, Al, Ti, O, etc.).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Composition {
    /// Fe mass fraction.
    pub fe_balance: f64,
    /// Al mass fraction as element.
    pub al: f64,
    /// Al mass fraction when used as a balance element.
    pub al_balance: f64,
    /// Ti mass fraction.
    pub ti_balance: f64,
    /// O mass fraction.
    pub o_balance: f64,
    /// Cr mass fraction.
    pub cr: f64,
    /// Ni mass fraction.
    pub ni: f64,
    /// Mo mass fraction.
    pub mo: f64,
    /// Mg mass fraction.
    pub mg: f64,
    /// Si mass fraction.
    pub si: f64,
    /// V mass fraction.
    pub v: f64,
    /// Mn mass fraction.
    pub mn: f64,
    /// Co mass fraction.
    pub co: f64,
    /// Li mass fraction.
    pub li: f64,
    /// C mass fraction.
    pub c: f64,
}

/// Mechanical properties at room temperature.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Mechanical {
    /// Young's modulus (GPa).
    pub youngs_modulus_gpa: f64,
    /// Yield strength (MPa).
    pub yield_strength_mpa: f64,
    /// Ultimate tensile strength (MPa).
    pub ultimate_tensile_strength_mpa: f64,
    /// Density (kg/m³).
    pub density_kg_per_m3: f64,
    /// Poisson's ratio.
    pub poissons_ratio: f64,
}

/// Thermal properties.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Thermal {
    /// Thermal conductivity (W/(m·K)).
    pub thermal_conductivity_w_per_mk: f64,
    /// Specific heat at constant pressure (J/(kg·K)).
    pub specific_heat_j_per_kgk: f64,
    /// Linear coefficient of thermal expansion (1/K).
    pub coefficient_thermal_expansion_per_k: f64,
    /// Melting point (K).
    pub melting_point_k: f64,
}

/// Electrical properties.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Electrical {
    /// Electrical resistivity (Ω·m).
    pub electrical_resistivity_ohm_m: f64,
}

/// One row of the materials database.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MaterialRecord {
    /// Common name (e.g. `"AISI 304"`).
    pub name: String,
    /// Composition (mass fractions).
    pub composition: Composition,
    /// Mechanical properties.
    pub mechanical: Mechanical,
    /// Thermal properties.
    pub thermal: Thermal,
    /// Electrical properties.
    pub electrical: Electrical,
    /// Provenance for these values.
    pub sources: Vec<DataSource>,
}
