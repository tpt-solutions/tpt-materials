//! Hydrogel swelling equilibrium and kinetics.
//!
//! - [`flory_rehner_swelling_ratio`] — solves the Flory–Rehner
//!   equilibrium `π_mix + π_elastic = 0`.
//! - [`equilibrium_swelling_ratio`] — convenience wrapper that
//!   accepts `χ` and crosslink density.
//! - [`poro_uptake_curve`] — Fickian uptake on a slab.

#![warn(missing_docs)]

mod flory_rehner;
mod kinetics;

pub use flory_rehner::{
    equilibrium_swelling_ratio, flory_rehner_swelling_ratio, FloryRehnerParams,
};
pub use kinetics::{poro_uptake_curve, slab_fickian_uptake};

use serde::{Deserialize, Serialize};

/// Swelling ratio `Q = V_swollen / V_dry` (≥ 1).
pub type SwellingRatio = f64;