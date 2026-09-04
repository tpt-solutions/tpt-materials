//! `tpt-materials` — umbrella / facade crate.
//!
//! Every spec §5 domain crate is re-exported behind an optional
//! feature so callers can pull in only what they need:
//!
//! ```toml
//! [dependencies]
//! tpt-materials = { version = "0.1", features = ["full"] }
//! ```
//!
//! or selectively:
//!
//! ```toml
//! [dependencies]
//! tpt-materials = { version = "0.1", features = ["crystal-plasticity", "homogenization"] }
//! ```
//!
//! # Spec §6 / §13 import snippets
//!
//! ```ignore
//! use tpt_materials::crystallography::CrystalStructure;
//! use tpt_materials::crystal_plasticity::{CrystalPlasticityModel, RateSensitivity};
//! use tpt_materials::homogenization::voigt;
//! use tpt_materials::hardening::{Hardening, VoceHardening, VoceParams};
//! ```
//!
//! Lightweight runtime smoke-test of the facade:
//!
//! ```
//! # #[cfg(feature = "crystallography")]
//! # {
//! use tpt_materials::crystallography::CrystalStructure;
//! let n_slip = CrystalStructure::FCC.slip_systems().len();
//! assert_eq!(n_slip, 12);
//! # }
//! ```

#![warn(missing_docs)]

#[cfg(feature = "additive")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_additive;
#[cfg(feature = "battery")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_battery;
#[cfg(feature = "calphad")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_calphad;
#[cfg(feature = "composite-micro")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_composite_micro;
#[cfg(feature = "corrosion")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_corrosion;
#[cfg(feature = "creep")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_creep;
#[cfg(feature = "crystal-plasticity")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_crystal_plasticity;
#[cfg(feature = "crystallography")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_crystallography;
#[cfg(feature = "damage")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_damage;
#[cfg(feature = "database")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_database;
#[cfg(feature = "diffusion")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_diffusion;
#[cfg(feature = "dislocation")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_dislocation;
#[cfg(feature = "fatigue")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_fatigue;
#[cfg(feature = "fatigue-micro")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_fatigue_micro;
#[cfg(feature = "fracture")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_fracture;
#[cfg(feature = "grain-growth")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_grain_growth;
#[cfg(feature = "hardening")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_hardening;
#[cfg(feature = "heat-treatment")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_heat_treatment;
#[cfg(feature = "homogenization")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_homogenization;
#[cfg(feature = "hydrogel")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_hydrogel;
#[cfg(feature = "hydrogen-embrittlement")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_hydrogen_embrittlement;
#[cfg(feature = "hydrogen-storage")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_hydrogen_storage;
#[cfg(feature = "machine-learning")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_machine_learning;
#[cfg(feature = "phase-field")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_phase_field;
#[cfg(feature = "phase-transform")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_phase_transform;
#[cfg(feature = "polymer")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_polymer;
#[cfg(feature = "precipitation")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_precipitation;
#[cfg(feature = "rve")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_rve;
#[cfg(feature = "solidification")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_solidification;
#[cfg(feature = "texture")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_texture;
#[cfg(feature = "thermal")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_thermal;
#[cfg(feature = "wasm")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_wasm;
#[cfg(feature = "welding")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_welding;
#[cfg(feature = "inverse")]
#[allow(unused_extern_crates)]
extern crate tpt_mat_inverse;

#[cfg(feature = "crystallography")]
pub mod crystallography {
    //! Re-export of `tpt-mat-crystallography`.
    pub use tpt_mat_crystallography::*;
}

#[cfg(feature = "crystal-plasticity")]
pub mod crystal_plasticity {
    //! Re-export of `tpt-mat-crystal-plasticity`.
    pub use tpt_mat_crystal_plasticity::*;
}

#[cfg(feature = "hardening")]
pub mod hardening {
    //! Re-export of `tpt-mat-hardening`.
    pub use tpt_mat_hardening::*;
}

#[cfg(feature = "texture")]
pub mod texture {
    //! Re-export of `tpt-mat-texture`.
    pub use tpt_mat_texture::*;
}

#[cfg(feature = "phase-field")]
pub mod phase_field {
    //! Re-export of `tpt-mat-phase-field`.
    pub use tpt_mat_phase_field::*;
}

#[cfg(feature = "grain-growth")]
pub mod grain_growth {
    //! Re-export of `tpt-mat-grain-growth`.
    pub use tpt_mat_grain_growth::*;
}

#[cfg(feature = "solidification")]
pub mod solidification {
    //! Re-export of `tpt-mat-solidification`.
    pub use tpt_mat_solidification::*;
}

#[cfg(feature = "diffusion")]
pub mod diffusion {
    //! Re-export of `tpt-mat-diffusion`.
    pub use tpt_mat_diffusion::*;
}

#[cfg(feature = "phase-transform")]
pub mod phase_transform {
    //! Re-export of `tpt-mat-phase-transform`.
    pub use tpt_mat_phase_transform::*;
}

#[cfg(feature = "calphad")]
pub mod calphad {
    //! Re-export of `tpt-mat-calphad`.
    pub use tpt_mat_calphad::*;
}

#[cfg(feature = "homogenization")]
pub mod homogenization {
    //! Re-export of `tpt-mat-homogenization`.
    pub use tpt_mat_homogenization::*;
}

#[cfg(feature = "rve")]
pub mod rve {
    //! Re-export of `tpt-mat-rve`.
    pub use tpt_mat_rve::*;
}

#[cfg(feature = "composite-micro")]
pub mod composite_micro {
    //! Re-export of `tpt-mat-composite-micro`.
    pub use tpt_mat_composite_micro::*;
}

#[cfg(feature = "damage")]
pub mod damage {
    //! Re-export of `tpt-mat-damage`.
    pub use tpt_mat_damage::*;
}

#[cfg(feature = "fatigue")]
pub mod fatigue {
    //! Re-export of `tpt-mat-fatigue`.
    pub use tpt_mat_fatigue::*;
}

#[cfg(feature = "fatigue-micro")]
pub mod fatigue_micro {
    //! Re-export of `tpt-mat-fatigue-micro`.
    pub use tpt_mat_fatigue_micro::*;
}

#[cfg(feature = "corrosion")]
pub mod corrosion {
    //! Re-export of `tpt-mat-corrosion`.
    pub use tpt_mat_corrosion::*;
}

#[cfg(feature = "creep")]
pub mod creep {
    //! Re-export of `tpt-mat-creep`.
    pub use tpt_mat_creep::*;
}

#[cfg(feature = "heat-treatment")]
pub mod heat_treatment {
    //! Re-export of `tpt-mat-heat-treatment`.
    pub use tpt_mat_heat_treatment::*;
}

#[cfg(feature = "welding")]
pub mod welding {
    //! Re-export of `tpt-mat-welding`.
    pub use tpt_mat_welding::*;
}

#[cfg(feature = "battery")]
pub mod battery {
    //! Re-export of `tpt-mat-battery`.
    pub use tpt_mat_battery::*;
}

#[cfg(feature = "additive")]
pub mod additive {
    //! Re-export of `tpt-mat-additive`.
    pub use tpt_mat_additive::*;
}

#[cfg(feature = "database")]
pub mod database {
    //! Re-export of `tpt-mat-database`.
    pub use tpt_mat_database::*;
}

#[cfg(feature = "machine-learning")]
pub mod machine_learning {
    //! Re-export of `tpt-mat-machine-learning`.
    pub use tpt_mat_machine_learning::*;
}

#[cfg(feature = "fracture")]
pub mod fracture {
    //! Re-export of `tpt-mat-fracture`.
    pub use tpt_mat_fracture::*;
}

#[cfg(feature = "thermal")]
pub mod thermal {
    //! Re-export of `tpt-mat-thermal`.
    pub use tpt_mat_thermal::*;
}

#[cfg(feature = "dislocation")]
pub mod dislocation {
    //! Re-export of `tpt-mat-dislocation`.
    pub use tpt_mat_dislocation::*;
}

#[cfg(feature = "precipitation")]
pub mod precipitation {
    //! Re-export of `tpt-mat-precipitation`.
    pub use tpt_mat_precipitation::*;
}

#[cfg(feature = "hydrogen-embrittlement")]
pub mod hydrogen_embrittlement {
    //! Re-export of `tpt-mat-hydrogen-embrittlement`.
    pub use tpt_mat_hydrogen_embrittlement::*;
}

#[cfg(feature = "polymer")]
pub mod polymer {
    //! Re-export of `tpt-mat-polymer`.
    pub use tpt_mat_polymer::*;
}

#[cfg(feature = "hydrogel")]
pub mod hydrogel {
    //! Re-export of `tpt-mat-hydrogel`.
    pub use tpt_mat_hydrogel::*;
}

#[cfg(feature = "hydrogen-storage")]
pub mod hydrogen_storage {
    //! Re-export of `tpt-mat-hydrogen-storage`.
    pub use tpt_mat_hydrogen_storage::*;
}

#[cfg(feature = "wasm")]
pub mod wasm_bindings {
    //! Re-export of `tpt-mat-wasm`.
    pub use tpt_mat_wasm::*;
}

#[cfg(feature = "inverse")]
pub mod inverse {
    //! Re-export of `tpt-mat-inverse` (closed-form calibration
    //! helpers + Levenberg–Marquardt driver).
    pub use tpt_mat_inverse::*;
}

/// Cross-repo output adapters (spec §6).
///
/// These functions convert in-repository material-model outputs into
/// the wire-format types consumed by the sibling `tpt-energy`,
/// `tpt-transport`, `tpt-electronics`, and `tpt-medical` crates.
/// They are deliberately lightweight — each adapter is a small
/// struct that bundles a few scalar/vector fields plus provenance
/// metadata, and is `serde::Serialize` so it can be JSON-encoded at
/// the cross-repo boundary.
#[cfg(feature = "adapters")]
pub mod adapters {
    use serde::{Deserialize, Serialize};

    /// Provenance metadata attached to every cross-repo adapter.
    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    pub struct AdapterProvenance {
        /// Source crate name (e.g. `tpt-mat-battery`).
        pub source_crate: String,
        /// Source model version (semver).
        pub source_version: String,
        /// ISO-8601 generation timestamp.
        pub generated_at: String,
    }

    /// `tpt-energy` adapter: battery capacity-fade curve.
    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    pub struct EnergyDegradationAdapter {
        /// Provenance.
        pub provenance: AdapterProvenance,
        /// Cell chemistry label (e.g. "NMC811").
        pub chemistry: String,
        /// Cycle indices (1..=n).
        pub cycles: Vec<usize>,
        /// Capacity retention (0..=1).
        pub capacity_retention: Vec<f64>,
        /// Internal-resistance growth (relative).
        pub resistance_growth: Vec<f64>,
        /// Reference C-rate used for the curve.
        pub c_rate: f64,
        /// Reference temperature in K.
        pub temperature_k: f64,
    }

    impl EnergyDegradationAdapter {
        /// Build from a [`tpt_mat_battery::DegradationCurve`].
        #[cfg(feature = "battery")]
        pub fn from_degradation_curve(
            curve: &tpt_mat_battery::DegradationCurve,
            chemistry: impl Into<String>,
            c_rate: f64,
            temperature_k: f64,
            source_version: impl Into<String>,
        ) -> Self {
            Self {
                provenance: AdapterProvenance {
                    source_crate: "tpt-mat-battery".to_string(),
                    source_version: source_version.into(),
                    generated_at: "1970-01-01T00:00:00Z".to_string(),
                },
                chemistry: chemistry.into(),
                cycles: curve.cycles.clone(),
                capacity_retention: curve.capacity_retention.clone(),
                resistance_growth: curve.resistance_growth.clone(),
                c_rate,
                temperature_k,
            }
        }
    }

    /// `tpt-transport` adapter: composite fatigue S-N curve.
    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    pub struct TransportFatigueAdapter {
        /// Provenance.
        pub provenance: AdapterProvenance,
        /// Material label.
        pub material: String,
        /// Stress amplitude in MPa.
        pub stress_amplitude_mpa: Vec<f64>,
        /// Cycles to failure.
        pub cycles_to_failure: Vec<f64>,
        /// Basquin exponent.
        pub basquin_b: f64,
        /// Fatigue strength coefficient (MPa).
        pub fatigue_strength_coeff_mpa: f64,
    }

    /// `tpt-transport` adapter: alloy creep strain rate.
    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    pub struct TransportCreepAdapter {
        /// Provenance.
        pub provenance: AdapterProvenance,
        /// Alloy name.
        pub alloy: String,
        /// Temperature in K.
        pub temperature_k: f64,
        /// Applied stress in MPa.
        pub stress_mpa: f64,
        /// Steady-state creep strain rate (1/s).
        pub strain_rate: f64,
        /// Norton-Bailey exponent.
        pub norton_n: f64,
    }

    /// `tpt-electronics` adapter: solder-joint fatigue cycles.
    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    pub struct ElectronicsSolderAdapter {
        /// Provenance.
        pub provenance: AdapterProvenance,
        /// Solder alloy (e.g. "SAC305").
        pub solder: String,
        /// Shear strain range.
        pub shear_strain_range: f64,
        /// Predicted cycles to failure (Coffin-Manson).
        pub cycles_to_failure: f64,
        /// Fatigue exponent.
        pub coffin_manson_c: f64,
    }

    /// `tpt-medical` adapter: implant corrosion / biocompatibility.
    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    pub struct MedicalCorrosionAdapter {
        /// Provenance.
        pub provenance: AdapterProvenance,
        /// Implant alloy (e.g. "Ti-6Al-4V").
        pub alloy: String,
        /// Penetration rate in mm/yr.
        pub penetration_rate_mm_per_yr: f64,
        /// Mass-loss rate in g/(m^2·day).
        pub mass_loss_rate_g_per_m2_day: f64,
        /// Body environment (e.g. "saline-37C").
        pub environment: String,
    }
}

/// Crate version string.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_non_empty() {
        assert!(!VERSION.is_empty());
    }
}
