//! Phase-field solvers on regular grids.
//!
//! Provides:
//! - [`PhaseFieldSolver`]: explicit first-order time integration of
//!   [`PhaseFieldModel::AllenCahn`], [`PhaseFieldModel::CahnHilliard`],
//!   [`PhaseFieldModel::Kobayashi`] (Allen-Cahn with latent-heat
//!   thermal coupling), and [`PhaseFieldModel::MultiPhase`].
//! - [`FreeEnergyFunctional`] + [`BulkEnergy`] (`DoubleWell`,
//!   `Polynomial`, `RegularSolution`).
//! - [`PhaseFieldResult`] snapshot after each step.
//!
//! Discretisation:
//! - Spatial: central finite differences on the [`Grid2D`] /
//!   [`Grid3D`] supplied by `tpt-science` (Neumann zero-flux
//!   boundaries).
//! - Temporal: forward Euler with explicit stability limits
//!   (`Δt ≤ Δx⁴ / (κ D)` for Cahn-Hilliard, `Δt ≤ Δx² / (2κ)` for
//!   Allen-Cahn).
//!
//! [`Grid2D`]: tpt_science::Grid2D
//! [`Grid3D`]: tpt_science::Grid3D

#![warn(missing_docs)]

mod energy;
mod model;
mod result;
mod solver;

pub use energy::{BulkEnergy, FreeEnergyFunctional, RegularSolutionParams};
pub use model::PhaseFieldModel;
pub use result::PhaseFieldResult;
pub use solver::PhaseFieldSolver;
