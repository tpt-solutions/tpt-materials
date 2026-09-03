//! Single-crystal plasticity constitutive model and per-integration-point
//! update.
//!
//! This crate implements the small-strain, rate-dependent crystal
//! plasticity flow rule of Peirce, Asaro & Needleman (1982):
//!
//! `γ̇^α = γ̇_0 |τ^α / τ_c^α|^n sign(τ^α)`
//!
//! where `τ^α = σ : P^α` is the resolved shear stress, `τ_c^α` is the
//! slip-system CRSS (evolved by a [`HardeningLaw`](tpt_mat_hardening::HardeningLaw)),
//! `γ̇_0` is a reference slip rate and `n` is the rate sensitivity.
//!
//! The plastic velocity gradient follows as
//!
//! `L^p = Σ_α γ̇^α s^α ⊗ n^α`
//!
//! and (when combined with a stress integration scheme and the elastic
//! stiffness `C`) drives an FEM Newton-Raphson solve.  FEM assembly
//! itself lives in the cross-repo `tpt-fem` substrate; this crate
//! provides everything needed to evaluate the constitutive model at a
//! single integration point.

#![warn(missing_docs)]

mod elastic;
mod fem;
mod flow;
mod model;
mod state;

pub use elastic::{ElasticStiffness, ElasticStiffnessError, SymmetricFourthOrder};
pub use fem::{
    solve_increment_single_point, BoundaryConditions, CpFemResult, CpFemSolver, FemError, LoadStep,
    ReactionForce, StressUpdate,
};
pub use flow::{power_law_slip_rate, viscoplastic_velocity_gradient, PlasticIncrement};
pub use model::{CrystalPlasticityModel, RateSensitivity};
pub use state::SingleCrystalState;
