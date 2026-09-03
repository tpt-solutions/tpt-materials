//! Heat-treatment simulation pipeline.
//!
//! This crate composes the analytical Phase 4 kinetic models
//! (`tpt-mat-phase-transform`) with simple empirical hardness
//! regressions to provide a first-pass heat-treatment simulator.
//!
//! ## Process taxonomy
//!
//! - **Annealing** — slow cool from austenitising temperature;
//!   equilibrium phase fractions from CALPHAD, grain growth via
//!   `tpt-mat-grain-growth` is *not* modelled here.
//! - **Quenching** — rapid cool to a quench medium temperature;
//!   austenite → martensite via Koistinen–Marburger plus any
//!   diffusion-controlled product formed during the slower part
//!   of the cool.
//! - **Tempering** — reheat below `A_1` after quench; martensite
//!   fraction monotonically decreases with hold time.
//! - **Aging** — low-temperature precipitation-controlled
//!   hardening; uses an Avrami precipitation-style kinetic that
//!   only depends on a "peak-age" parameter.
//! - **Solution treatment** — hold above solvus; reports final
//!   composition-averaged phase fractions.
//!
//! ## Hardness model
//!
//! Following the Maynier-style regression tradition:
//!
//! ```text
//! HV = Σ_i f_i HV_i + ΔHV_age
//! ```
//!
//! where `f_i` are the final phase fractions and `HV_i` are
//! representative hardness values for each phase
//! (martensite ≈ 800 HV, bainite ≈ 450 HV, pearlite ≈ 250 HV,
//! ferrite ≈ 100 HV, austenite ≈ 200 HV).

#![warn(missing_docs)]

mod hardness;
mod process;

pub use hardness::{hardness_from_fractions, hardness_martensite, PhaseHardness};
pub use process::{
    simulate, AgingParams, AnnealingParams, HeatTreatmentProcess, HeatTreatmentResult,
    QuenchingParams, TemperingParams,
};
