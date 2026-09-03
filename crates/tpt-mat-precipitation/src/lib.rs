//! Precipitation kinetics: classical nucleation, KWN size-class
//! solver, and LSW coarsening.
//!
//! - [`classical_nucleation_rate`] — Turnbull–Fisher
//!   `I = I_0 exp(−ΔG*/kT) exp(−Q/kT)` with the spherical-cap
//!   barrier `ΔG* = 16π γ³ / (3 Δg_v²)` (for a stress-free,
//!   coherent nucleus).
//! - [`KwnState`] / [`kwn_step`] — discretised size-class
//!   Kampmann–Wagner solver (1-D mass balance).
//! - [`lsw_coarsening_rate`] — `dR³/dt = K`, `K` from
//!   `γ D c_eq / V_m`.
//! - [`strengthening_increment`] — Orowan bypass + shearing
//!   hand-off to [`tpt_mat_hardening`].

#![warn(missing_docs)]

mod cnt;
mod kwn;
mod lsw;
mod strengthening;

pub use cnt::{classical_nucleation_rate, critical_radius, nucleation_barrier};
pub use kwn::{kwn_step, KwnParams, KwnState};
pub use lsw::{lsw_coarsening_rate, lsw_radius_cubed_growth};
pub use strengthening::{
    orowan_bypass_strengthening, shearing_strengthening,
    StrengtheningIncrement, StrengtheningMechanism,
};