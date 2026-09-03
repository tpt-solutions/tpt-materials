//! Continuum damage mechanics.
//!
//! This crate implements the classical analytical damage models
//! required by Phase 6:
//!
//! - **Effective-stress concept (Kachanov 1958)**: the scalar
//!   damage variable `D ∈ [0, 1]` reduces the load-bearing area so
//!   that `σ̃ = σ / (1 - D)`.
//! - **Lemaitre (1971) coupled damage-viscoelasticity** for ductile
//!   damage: `Ḋ = ((σ̃_eq / σ_u)^s / B) (1 - D)^(-α)` where
//!   `σ̃_eq` is the equivalent stress in the *undamaged* configuration.
//! - **Kachanov (1958) creep-damage**: `Ḋ = A σ^n / (1 - D)^k`
//!   with rupture at `D = 1`.
//! - **Miner (1945) linear damage accumulation** for high-cycle
//!   fatigue: `D = Σ n_i / N_i` with failure at `D = 1`.
//! - **Chaboche (1977) nonlinear kinematic hardening** companion
//!   placeholder.
//!
//! # References
//!
//! - Kachanov, L. M. (1958).  "Time of the rupture process under
//!   creep conditions."  Izv. Akad. Nauk SSSR, Otd. Tekh. Nauk 8,
//!   26–31.
//! - Lemaitre, J. (1971).  "Evaluation of dissipation and damage
//!   in metals submitted to dynamic loading."  Proc. I.C.M. 1.
//! - Miner, M. A. (1945).  "Cumulative damage in fatigue."
//!   J. Appl. Mech. 12, A159–A164.
//!
//! # Inner attributes
#![warn(missing_docs)]

mod effective_stress;
mod kachanov;
mod lemaitre;
mod miner;

pub use effective_stress::{
    effective_stress, effective_youngs_modulus, stress_triaxiality,
};
pub use kachanov::{
    kachanov_damage_rate, kachanov_rupture_time, KachanovParams,
};
pub use lemaitre::{
    lemaitre_damage_rate, lemaitre_damage_step, LemaitreParams,
};
pub use miner::{
    miner_damage_accumulation, miner_remaining_life, MinerCyclicInputs,
};