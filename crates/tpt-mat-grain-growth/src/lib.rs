//! Multi-grain phase-field grain-growth simulation.
//!
//! Each grain is represented by a non-conserved order parameter `η_i`;
//! the model is the standard multi-phase Allen-Cahn
//! `∂η_i/∂t = −L_i (δF/δη_i)` with `Σ η_i = 1` enforced softly via a
//! pairwise coupling term.  Misorientation-dependent grain-boundary
//! mobility is supplied through [`GrainBoundaryMobility`].
//!
//! # Interfaces & grain boundaries
//!
//! The [`interfaces`] module provides consolidated models for
//! grain-boundary energy (`GrainBoundaryEnergy`), the grain-boundary
//! character distribution (`GrainBoundaryCharacterDistribution`),
//! triple-junction force balance (`TripleJunction`) and
//! Langmuir–McLean solute segregation (`LangmuirMcLean`).

#![warn(missing_docs)]

mod interfaces;
mod mobility;
mod recrystallization;
mod solver;
mod stereology;

pub use interfaces::{
    GBCDEntry, GrainBoundaryCharacterDistribution, GrainBoundaryEnergy, LangmuirMcLean,
    TripleJunction,
};
pub use mobility::{GrainBoundaryMobility, MobilityParams};
pub use recrystallization::{DrxKinetics, JmakRecrystallization, ZenerHollomon};
pub use solver::{GrainGrowthSolver, GrainSizeDistribution, GrainSizeStats};
pub use stereology::{
    saltykov_size_distribution, AreaFraction2D, LinearIntercept, NumberPerArea, VolumeFraction3D,
};
