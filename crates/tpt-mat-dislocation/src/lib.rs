//! Dislocation-density-based hardening.
//!
//! - [`DislocationDensityState`] — per-slip-system `ρ_SSD`,
//!   `ρ_GND`, and forest density.
//! - [`KocksMeckingEvolution`] — `dρ/dγ = k₁ √ρ − k₂ ρ` with
//!   Taylor stress `τ = α μ b √ρ`.
//! - [`gnd_from_curvature`] — Nye-tensor trace `||α||` to scalar
//!   GND density.
//! - [`back_stress`] — Armstrong–Frederick kinematic term from
//!   the GND density.
//!
//! These are intended to feed into [`tpt_mat_hardening`] as an
//! additional `HardeningLaw::DislocationDensity` variant.

#![warn(missing_docs)]

mod back_stress;
mod kocks_mecking;
mod nye;
mod state;

pub use back_stress::{armstrong_frederick_step, BackStressParams};
pub use kocks_mecking::{
    kocks_mecking_step, taylor_stress, KocksMeckingParams,
};
pub use nye::gnd_from_curvature;
pub use state::DislocationDensityState;