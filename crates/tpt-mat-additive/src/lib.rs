//! Additive-manufacturing process models.
//!
//! Analytical thermal-history, microstructure and residual-stress
//! predictors for laser powder-bed fusion (LPBF), directed-energy
//! deposition (DED) and electron-beam melting (EBM).
//!
//! # Models
//!
//! - [`AmProcess`] — process parameters (LPBF / DED / EBM).
//! - [`thermal_history`] — moving-source analytic solution
//!   (Rosenthal 1941 / Eagar–Tsai 1983 with finite-width correction
//!   for a Gaussian beam).
//! - [`predict_microstructure`] — columnar / equiaxed map from the
//!   G–R solidification map and resulting grain size (Hall–Petch).
//! - [`residual_stress`] — thermal-contraction eigenstrain field.
//!
//! # References
//!
//! - Rosenthal, D. (1941).  "Mathematical theory of heat
//!   distribution during welding and cutting."  Weld. J. 20,
//!   220s–234s.
//! - Eagar, T. W., & Tsai, N.-S. (1983).  "Temperature fields
//!   produced by traveling distributed heat sources."  Weld. J.
//!   62(12), 346s–355s.
//! - Kou, S. (2003).  *Welding Metallurgy.*  2nd ed., Wiley.

#![warn(missing_docs)]

mod history;
mod microstructure;
mod residual;

use serde::{Deserialize, Serialize};

pub use history::{
    cooling_rate, rosenthal_2d, thermal_history, CoolingRate, GaussianBeam, MovingPointSource,
    ThermalHistory,
};
pub use microstructure::{
    columnar_equiaxed_threshold, grain_size_from_cooling_rate, predict_microstructure,
    PredictedMicrostructure, SolidificationMap,
};
pub use residual::{
    residual_stress, residual_stress_field, ResidualStressField, ThermalContraction,
};

/// Additive-manufacturing process parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AmProcess {
    /// Laser powder-bed fusion (LPBF): thin layers (20–100 µm),
    /// small hatch spacing (~100 µm), moderate scan speed.
    LaserPowderBedFusion {
        /// Laser power `P` (W).
        laser_power: f64,
        /// Scan speed `v` (m/s).
        scan_speed: f64,
        /// Hatch spacing `h` (m).
        hatch_spacing: f64,
        /// Layer thickness `t` (m).
        layer_thickness: f64,
        /// Beam radius (m).
        beam_radius: f64,
    },
    /// Directed-energy deposition (DED): larger melt pool, lower
    /// scan speed.
    DirectedEnergyDeposition {
        /// Laser power `P` (W).
        laser_power: f64,
        /// Scan speed `v` (m/s).
        scan_speed: f64,
        /// Deposit width `w` (m).
        deposit_width: f64,
        /// Layer thickness `t` (m).
        layer_thickness: f64,
        /// Beam radius (m).
        beam_radius: f64,
    },
    /// Electron-beam melting (EBM): high scan speed, preheated bed.
    ElectronBeamMelting {
        /// Beam power `P` (W).
        beam_power: f64,
        /// Scan speed `v` (m/s).
        scan_speed: f64,
        /// Hatch spacing `h` (m).
        hatch_spacing: f64,
        /// Layer thickness `t` (m).
        layer_thickness: f64,
        /// Beam radius (m).
        beam_radius: f64,
    },
}

impl AmProcess {
    /// Effective heat input per unit length `Q = P / v` (J/m).
    pub fn linear_heat_input(&self) -> f64 {
        match self {
            Self::LaserPowderBedFusion {
                laser_power,
                scan_speed,
                ..
            } => *laser_power / *scan_speed,
            Self::DirectedEnergyDeposition {
                laser_power,
                scan_speed,
                ..
            } => *laser_power / *scan_speed,
            Self::ElectronBeamMelting {
                beam_power,
                scan_speed,
                ..
            } => *beam_power / *scan_speed,
        }
    }

    /// Volumetric energy density `E = P / (v h t)` (J/m³).
    pub fn volumetric_energy_density(&self) -> f64 {
        match self {
            Self::LaserPowderBedFusion {
                laser_power,
                scan_speed,
                hatch_spacing,
                layer_thickness,
                ..
            } => *laser_power / (*scan_speed * *hatch_spacing * *layer_thickness),
            Self::DirectedEnergyDeposition {
                laser_power,
                scan_speed,
                deposit_width,
                layer_thickness,
                ..
            } => *laser_power / (*scan_speed * *deposit_width * *layer_thickness),
            Self::ElectronBeamMelting {
                beam_power,
                scan_speed,
                hatch_spacing,
                layer_thickness,
                ..
            } => *beam_power / (*scan_speed * *hatch_spacing * *layer_thickness),
        }
    }
}

/// Pre-heated substrate / build-plate temperature (K) for EBM.
pub const EBM_PREHEAT_TEMPERATURE: f64 = 1023.0;
/// LPBF pre-heat (room temperature, K).
pub const LPBF_PREHEAT_TEMPERATURE: f64 = 298.0;
/// DED pre-heat (room temperature, K).
pub const DED_PREHEAT_TEMPERATURE: f64 = 298.0;
