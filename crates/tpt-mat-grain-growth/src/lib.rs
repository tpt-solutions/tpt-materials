//! Multi-grain phase-field grain-growth simulation.
//!
//! Each grain is represented by a non-conserved order parameter `η_i`;
//! the model is the standard multi-phase Allen-Cahn
//! `∂η_i/∂t = −L_i (δF/δη_i)` with `Σ η_i = 1` enforced softly via a
//! pairwise coupling term.  Misorientation-dependent grain-boundary
//! mobility is supplied through [`GrainBoundaryMobility`].

#![warn(missing_docs)]

mod mobility;
mod solver;

pub use mobility::{GrainBoundaryMobility, MobilityParams};
pub use solver::{GrainGrowthSolver, GrainSizeDistribution, GrainSizeStats};
