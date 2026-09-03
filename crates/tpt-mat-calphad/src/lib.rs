//! CALPHAD-style thermodynamic modelling.
//!
//! Implements the analytical building blocks used in classical
//! CALPHAD assessments:
//!
//! - [`RedlichKister`]: Redlich–Kister polynomial for the binary
//!   excess Gibbs energy `^E G = Σ L_ν (x_A − x_B)^ν`.
//! - [`GibbsEnergyModel`]: end-member + ideal-mixing + excess
//!   contributions assembled per phase.
//! - [`SublatticeModel`]: Muggianu-style sublattice phase with two or
//!   three sublattices, site fractions, and ideal configurational
//!   entropy.
//! - [`PhaseDiagram`]: two-phase equilibrium construction via
//!   common-tangent construction on `G(x)` at a fixed temperature
//!   (binary `T-x` slice), plus a temperature sweep that emits a
//!   `T-x` phase boundary.
//!
//! These are deliberately minimal models — enough to plot a
//! `G-x` curve, find the miscibility gap at a fixed `T`, and emit a
//! `T-x` phase boundary.  Full TDB / PAC parsing is deferred; the
//! intended consumer for this crate is the diffusion / transformation
//! crate that drives CALPHAD-derived mobilities into
//! [`tpt_mat_diffusion`](https://docs.rs/tpt-mat-diffusion).

#![warn(missing_docs)]

mod gibbs;
mod phase_diagram;
mod redlich_kister;
mod sublattice;

pub use gibbs::{EndMember, GibbsEnergyModel, IdealMixing};
pub use phase_diagram::{two_phase_equilibrium, PhaseBoundary, PhaseDiagram, TwoPhaseEquilibrium};
pub use redlich_kister::{RedlichKister, RedlichKisterParams};
pub use sublattice::{Sublattice, SublatticeConfig, SublatticeModel, SublatticeSpecies};
