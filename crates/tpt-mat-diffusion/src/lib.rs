//! Diffusion solvers and supporting thermodynamics.
//!
//! Implements Fick's second law on regular 1D / 2D / 3D Cartesian grids
//! for both single-component and multi-component substitutional
//! systems.  Diffusivities follow an Arrhenius law and grain-boundary
//! short-circuit diffusion is exposed via
//! [`GrainBoundaryDiffusion::effective_diffusivity`].
//!
//! The single-component solver is the textbook explicit forward-Euler
//! time-stepper:
//!
//! `∂c / ∂t = D ∇²c`
//!
//! with Neumann zero-flux boundaries (handled by `tpt-science`).  The
//! multi-component solver uses the same stencil but with a constant
//! diagonal diffusivity matrix (no cross-term; an
//! `OffDiagonalDiffusivity` pluggable point is reserved for Phase 4
//! expansion once CALPHAD mobility data is wired in).
//!
//! Stability limit (explicit Euler): `Δt ≤ Δx² / (2 d D)` where `d` is
//! the spatial dimension.

#![warn(missing_docs)]

mod arrhenius;
mod grain_boundary;
mod multicomponent;
mod scalar;

pub use arrhenius::{ArrheniusDiffusivity, ArrheniusParams};
pub use grain_boundary::{GBDiffusivity, GrainBoundaryDiffusion};
pub use multicomponent::{MultiComponentDiffusionSolver, MultiComponentError};
pub use scalar::{DiffusionSolver, DiffusionSolverError, DiffusionStepResult};
