//! Battery active-material modelling.
//!
//! This crate implements the classical models for
//! lithium-ion-battery degradation:
//!
//! - **Li diffusion in a spherical particle** (Newman 1963):
//!   `∂c/∂t = D/r² ∂/∂r (r² ∂c/∂r)` solved on a 1-D radial grid.
//! - **SEI growth** (single-particle, Peled 1979):
//!   `dδ/dt = k · exp(−Q / RT)` (Arrhenius) on the particle surface;
//!   capacity loss scales linearly with `δ`.
//! - **Particle cracking** (descriptor): cracking initiates when
//!   the diffusion-induced hoop stress exceeds a critical value.
//! - **Lithium plating** at low temperature / high C-rate: triggered
//!   when the surface overpotential falls below `0 V` vs Li/Li⁺.
//! - **Transition-metal dissolution**: a first-order decay of the
//!   active material.
//!
//! These are composed to a [`DegradationCurve`] that the `tpt-energy`
//! crate consumes (spec §6).
//!
//! # References
//!
//! - Newman, J. (1963).  "The Electrochemical Reaction."  (PhD
//!   thesis, UC Berkeley).  Spherical-particle diffusion model.
//! - Peled, E. (1979).  "The Electrochemical Behavior of Alkali and
//!   Alkaline Earth Metals in Nonaqueous Battery Systems."  J.
//!   Electrochem. Soc. 126(12), 2047–2051.
//! - Christensen, J., & Newman, J. (2004).  "Stress generation and
//!   fracture in lithium insertion materials."  J. Solid State
//!   Electrochem. 10, 293–319.

#![warn(missing_docs)]

mod active_material;
mod chemistry;
mod degradation;
mod diffusion;
mod fade;

use serde::{Deserialize, Serialize};

pub use active_material::ActiveMaterial;
pub use chemistry::{default_active_material, BatteryChemistry};
pub use degradation::{
    capacity_loss_from_sei, lithium_plating_risk, particle_cracking_risk,
    transition_metal_dissolution, DegradationMechanism, ParticleState,
};
pub use diffusion::{simulate_diffusion, DiffusionParams, DiffusionResult};
pub use fade::{capacity_fade_curve, simulate_diffusion_stress, DegradationCurve};

/// Library-wide error type.
#[derive(Debug, thiserror::Error)]
pub enum BatteryError {
    /// A required parameter is non-positive.
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),
}
