//! Solid-state phase-transformation kinetics.
//!
//! This crate implements the workhorse analytical models for
//! diffusion-controlled and diffusionless solid-state transformations:
//!
//! - [`AvramiModel`]: classical JMAK (Johnson–Mehl–Avrami–Kolmogorov)
//   isothermal transformation kinetics,
//!   `f(t) = 1 − exp(−k t^n)`, with the standard rate-temperature
//!   coupling `k(T)` (e.g. C-curve / TTT diagram).
//! - [`KoistinenMarburger`]: diffusionless (martensitic)
//!   transformation `f = 1 − exp(−α (M_s − T))` for `T ≤ M_s`.
//! - [`TransformationSolver`]: combined isothermal ↔ continuous-cooling
//!   simulation driver; integrates `df/dt` over a piecewise-linear
//!   thermal history using the Scheil additivity rule.
//!
//! These models are intentionally analytical (closed-form) rather than
//! numerical — that lets them serve as a *reference* solver against
//! which any future phase-field coupling (Phase 4+) can be
//! benchmarked.  The user-facing API is fully serialisable so that
//! TTT/CCT diagrams can be plotted directly from the JSON output.

#![warn(missing_docs)]

mod avrami;
mod koistinen;
mod solver;

pub use avrami::{AvramiModel, AvramiParams, TttPoint};
pub use koistinen::{KMParams, KoistinenMarburger};
pub use solver::{
    ThermalHistoryStep, TransformationError, TransformationSolver, TransformationState,
};
