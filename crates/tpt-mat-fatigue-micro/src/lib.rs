//! Fatigue-indicator parameters (FIPs) for microstructural fatigue.
//!
//! This crate sits at the boundary between the crystal-plasticity
//! constitutive model ([`tpt_mat_crystal_plasticity`]) and the
//! microstructural RVE ([`tpt_mat_rve`]).  Given the converged state
//! of a cycle-by-cycle CP-FEM simulation, it computes the
//! *fatigue-indicator parameter* (FIP) field — a scalar field whose
//! maximum over the RVE identifies the most likely crack-initiation
//! site.
//!
//! The four classical FIPs are implemented:
//!
//! - **Findley (1957)**: shear-stress amplitude on the maximum
//!   shear plane, `F = max_θ (Δτ_max(θ)/2 + k σ_n_max(θ))`.
//! - **Fatemi–Socie (1988)**: combines the shear strain amplitude
//!   with the peak normal stress on the maximum-shear plane,
//!   `F = (Δγ_max / 2) (1 + k ⟨σ_n_max⟩ / σ_y)`.
//! - **Smith–Watson–Topper (1970)**:
//!   `F = σ_max · (Δε / 2)`.
//! - **Crystallographic slip** (Tanaka–Mura 1981): critical
//!   accumulated shear `γ_c`; initiation when
//!   `Σ_α γ_acc^α ≥ γ_c`.
//!
//! # Pipeline
//!
//! 1. Run a cycle-by-cycle CP-FEM simulation on an RVE (one grain per
//!    integration point, using [`tpt_mat_crystal_plasticity::CpFemSolver`]).
//! 2. At each load step extract
//!    [`tpt_mat_crystal_plasticity::CpFemResult`] (stresses +
//!    accumulated shear per slip system).
//! 3. Compute the FIP field with [`fatigue_indicator_parameter`].
//! 4. Drive [`predict_crack_initiation`] to identify the critical
//!    grain and estimate cycles-to-initiation.
//!
//! # References
//!
//! - Findley, W. N. (1957).  "Fatigue of metals under combined
//!   bending and torsion."  Proc. ASTM 57, 880–886.
//! - Fatemi, A., & Socie, D. F. (1988).  "A critical plane approach
//!   to multiaxial fatigue damage including out-of-phase loading."
//!   Fatigue Fract. Eng. Mater. Struct. 11(3), 149–165.
//! - Smith, K. N., Watson, P., & Topper, T. H. (1970).  "A
//!   stress-strain function for the fatigue of metals."  J. Mater.
//!   5(4), 767–778.
//! - Tanaka, K., & Mura, T. (1981).  "A dislocation model for
//!   fatigue crack initiation."  J. Appl. Mech. 48(1), 97–103.

#![warn(missing_docs)]

mod crystallographic;
mod fatemi_socie;
mod findley;
mod fip_field;
mod prediction;
mod swt;

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_mat_rve::Rve;

/// Fatigue-indicator parameter (FIP) criterion.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum FatigueCriterion {
    /// Findley (1957): critical-plane search over stress states.
    Findley {
        /// Sensitivity coefficient `k` (typical `0.3` for steels).
        k: f64,
    },
    /// Fatemi–Socie (1988): strain-based, normal-stress amplified.
    FatemiSocie {
        /// Yield stress `σ_y` (MPa).
        sigma_y: f64,
        /// Sensitivity coefficient `k` (typical `0.5–1.0`).
        k: f64,
    },
    /// Smith–Watson–Topper (1970): energy-based, `σ_max · (Δε / 2)`.
    SmithWatsonTopper,
    /// Crystallographic-slip (Tanaka–Mura 1981): critical accumulated
    /// shear on a slip system.
    CrystallographicSlip {
        /// Critical accumulated shear `γ_c` for crack initiation.
        critical_accumulated_shear: f64,
    },
}

impl Default for FatigueCriterion {
    fn default() -> Self {
        Self::Findley { k: 0.3 }
    }
}

/// Microstructural-fatigue problem: an RVE + a chosen FIP criterion.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MicrostructuralFatigue {
    /// The representative volume element (per-grain orientation + stiffness).
    pub rve: Rve,
    /// FIP criterion.
    pub criterion: FatigueCriterion,
}

impl MicrostructuralFatigue {
    /// Construct a microstructural-fatigue problem.
    pub fn new(rve: Rve, criterion: FatigueCriterion) -> Self {
        Self { rve, criterion }
    }
}

/// Result of a crack-initiation prediction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CrackInitiationResult {
    /// Estimated number of cycles to crack initiation at the critical
    /// grain.  `∞` if no initiation is predicted within the supplied
    /// load history.
    pub cycles_to_initiation: f64,
    /// Index of the critical grain (highest FIP).
    pub critical_grain: usize,
    /// Maximum FIP value reached (at the critical grain).
    pub critical_fip: f64,
    /// FIP field over all grains.
    pub fip_field: Vec<f64>,
}

/// Errors raised by the microstructural-fatigue pipeline.
#[derive(Debug, Error)]
pub enum FatigueMicroError {
    /// The RVE has zero grains.
    #[error("empty RVE")]
    EmptyRve,
}

/// Re-export the per-criterion FIP functions for convenience.
pub use crystallographic::{critical_shear_initiation, cycles_to_initiation_coffin_manson};
pub use fatemi_socie::{fatemi_socie_fip, fatemi_socie_search};
pub use findley::{findley_fip, findley_from_multiaxial_stress, findley_search};
pub use fip_field::fatigue_indicator_parameter;
pub use prediction::predict_crack_initiation;
pub use swt::swt_fip;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_criterion_is_findley() {
        let c = FatigueCriterion::default();
        assert!(matches!(c, FatigueCriterion::Findley { .. }));
    }
}
