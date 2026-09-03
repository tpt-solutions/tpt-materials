//! Welding process models.
//!
//! - [`WeldModel`] — base metal + filler + process.
//! - [`heat_affected_zone`] — HAZ width and peak-temperature profile
//!   from the Rosenthal line-source solution.
//! - [`predict_haz_microstructure`] — grain coarsening + phase
//!   transformation at the HAZ peak-temperature point, with the
//!   cooling-rate coupling from the heat-input.
//!
//! # References
//!
//! - Rosenthal, D. (1941).  "Mathematical theory of heat
//!   distribution during welding and cutting."  Weld. J. 20,
//!   220s–234s.
//! - Kou, S. (2003).  *Welding Metallurgy.*  2nd ed., Wiley.

#![warn(missing_docs)]

mod haz;
mod microstructure;

use serde::{Deserialize, Serialize};

pub use haz::{
    heat_affected_zone, peak_temperature_at, peak_temperature_profile, HazPeakProfile, HazResult,
    RosenthalWeld,
};
pub use microstructure::{predict_haz_microstructure, HazMicrostructure, MicrostructurePhase};

/// Filler metal chemistry descriptor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FillerMetal {
    /// Filler label (e.g. `"ER70S"`).
    pub label: String,
    /// Carbon content (wt%).
    pub carbon_wt_pct: f64,
    /// M_s temperature (K).
    pub ms_temperature: f64,
    /// Reference Avrami rate `k₀` at `T = 873 K` (1/s).
    pub avrami_k0: f64,
    /// Avrami time exponent.
    pub avrami_n: f64,
}

/// Base-metal descriptor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BaseMetal {
    /// Label.
    pub label: String,
    /// Carbon content (wt%).
    pub carbon_wt_pct: f64,
    /// A₁ / A₃ transformation temperatures (K).
    pub a1_temperature: f64,
    pub a3_temperature: f64,
    /// M_s temperature (K).
    pub ms_temperature: f64,
    /// Thermal conductivity (W/m·K).
    pub thermal_conductivity: f64,
    /// Thermal diffusivity (m²/s).
    pub thermal_diffusivity: f64,
    /// Pre-heat temperature (K).
    pub preheat_temperature: f64,
}

impl BaseMetal {
    /// Generic AISI 4140-style steel descriptor.
    pub fn aisi_4140() -> Self {
        Self {
            label: "AISI 4140".to_string(),
            carbon_wt_pct: 0.4,
            a1_temperature: 996.0,
            a3_temperature: 1073.0,
            ms_temperature: 600.0,
            thermal_conductivity: 42.0,
            thermal_diffusivity: 1.0e-5,
            preheat_temperature: 298.0,
        }
    }
}

/// Welding process parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum WeldProcess {
    /// Gas-metal arc welding (GMAW / MIG).
    Gmaw,
    /// Gas-tungsten arc welding (GTAW / TIG).
    Gtaw,
    /// Submerged-arc welding (SAW).
    Saw,
    /// Laser welding.
    Laser,
}

/// Weld-model bundle.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WeldModel {
    /// Base metal.
    pub base_metal: BaseMetal,
    /// Filler metal (None for autogenous).
    pub filler_metal: Option<FillerMetal>,
    /// Process parameters.
    pub process: WeldProcess,
    /// Effective linear heat input `Q = η V I / v` (J/m).
    pub linear_heat_input: f64,
    /// Scan speed `v` (m/s).
    pub scan_speed: f64,
    /// Process efficiency `η` (0–1).
    pub efficiency: f64,
}

impl WeldModel {
    /// Construct a weld model.
    pub fn new(
        base_metal: BaseMetal,
        filler_metal: Option<FillerMetal>,
        process: WeldProcess,
        linear_heat_input: f64,
        scan_speed: f64,
        efficiency: f64,
    ) -> Self {
        Self {
            base_metal,
            filler_metal,
            process,
            linear_heat_input,
            scan_speed,
            efficiency,
        }
    }
}
