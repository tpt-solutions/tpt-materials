//! Electrochemical corrosion models.
//!
//! This crate implements the mixed-potential theory of corrosion
//! following Wagner & Traud (1938) and the Butler–Volmer charge-
//! transfer kinetics.  The electrochemical models treat the
//! corroding metal as two coupled half-reactions (metal dissolution
//! and a cathodic process such as hydrogen evolution or oxygen
//! reduction) whose intersection on a polarization diagram gives the
//! corrosion current density `i_corr` and corrosion potential
//! `E_corr`.
//!
//! # Models
//!
//! - [`ElectrodeKinetics`]: a single half-reaction described by its
//!   Tafel slope `b` (V/dec), exchange current density `i_0`
//!   (A/m²) and equilibrium potential `E_eq` (V).
//! - [`CorrosionModel`]: a coupled anode + cathode on a common
//!   electrolyte.
//! - [`corrosion_rate`]: solve `i_a(E) = |i_c(E)|` for the corrosion
//!   current and derive a [`CorrosionRate`].
//! - [`polarization_curve`]: scan the potential axis and return the
//!   [`PolarizationCurve`] (anodic + cathodic branches).
//!
//! # References
//!
//! - Wagner, C., & Traud, W. (1938).  "Über die Deutung von
//!   Korrosionsvorgängen durch Überlagerung von
//!   elektrochemischen Teilvorgängen und über die Potentialbildung
//!   bei Mischelektroden."  Z. Elektrochem. 44(7), 391–402.
//! - Bard, A. J., & Faulkner, L. R. (2001).  *Electrochemical
//!   Methods: Fundamentals and Applications.*  2nd ed., Wiley.
//! - Jones, D. A. (1996).  *Principles and Prevention of
//!   Corrosion.*  2nd ed., Prentice-Hall.

#![warn(missing_docs)]

mod butler_volmer;
mod polarization;
mod rate;

use serde::{Deserialize, Serialize};

pub use butler_volmer::{
    butler_volmer_current, butler_volmer_current_density, tafel_anodic, tafel_cathodic,
    ElectrodeKinetics,
};
pub use polarization::{polarization_curve, PolarizationBranch, PolarizationCurve};
pub use rate::{
    corrosion_rate, mixed_potential, CorrosionModel, CorrosionRate, MixedPotentialResult,
};

/// Faraday constant (C/mol).
pub const FARADAY: f64 = 96_485.332_12;
/// Universal gas constant (J/mol·K).
pub const GAS_CONSTANT: f64 = 8.314_462_618;

/// Compute the molar mass of a metal `M` (kg/mol) given its atomic
/// mass and valence state.
pub fn molar_mass(atomic_mass: f64, valence: f64) -> f64 {
    atomic_mass / valence
}

/// Density of the corroding metal (kg/m³).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MetalDensity {
    /// Atomic mass (kg/mol, e.g. `0.055_845` for Fe).
    pub atomic_mass: f64,
    /// Metal valence for the dissolution reaction.
    pub valence: f64,
    /// Bulk density (kg/m³).
    pub density: f64,
}

impl MetalDensity {
    /// Iron `Fe → Fe²⁺ + 2 e⁻`.
    pub fn iron() -> Self {
        Self {
            atomic_mass: 0.055_845,
            valence: 2.0,
            density: 7874.0,
        }
    }
    /// Titanium `Ti → Ti⁴⁺ + 4 e⁻`.
    pub fn titanium() -> Self {
        Self {
            atomic_mass: 0.047_867,
            valence: 4.0,
            density: 4506.0,
        }
    }
    /// Aluminium `Al → Al³⁺ + 3 e⁻`.
    pub fn aluminium() -> Self {
        Self {
            atomic_mass: 0.026_982,
            valence: 3.0,
            density: 2700.0,
        }
    }
}
