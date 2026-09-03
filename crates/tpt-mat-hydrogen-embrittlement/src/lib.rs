//! Hydrogen embrittlement: transport with trapping, embrittlement
//! indicators, susceptibility.
//!
//! - [`hydrogen_diffusivity_effective`] — Oriani local-equilibrium
//!   `D_eff = D_L / (1 + ∂C_T/∂C_L)`.
//! - [`mcnabb_foster_kinetics`] — kinetic trapping rate
//!   `∂C_T/∂t = k (C_L (1 − θ_T) − θ_T / K)`, `K = exp(−E_B/RT)`.
//! - [`stress_driven_flux`] — hydrostatic-stress-driven
//!   `J = −D C V_H / (RT) ∇σ_h`.
//! - [`hede_threshold`] and [`help_threshold`] — HEDE / HELP
//!   embrittlement criteria.
//! - [`susceptibility_index`] — combined H + triaxiality indicator.

#![warn(missing_docs)]

mod embrittlement;
mod flux;
mod trapping;

pub use embrittlement::{
    hede_threshold, help_threshold, susceptibility_index, HedeParams, HelpParams,
};
pub use flux::stress_driven_flux;
pub use trapping::{
    hydrogen_diffusivity_effective, mcnabb_foster_trapping_rate, McNabbFosterParams,
    OrianiParams,
};