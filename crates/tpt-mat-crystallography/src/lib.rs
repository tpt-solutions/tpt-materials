//! Crystal structures, slip systems, and the Schmid tensor.
//!
//! Phase 1 covers the structures the rest of the Phase 1 verification
//! tests need: FCC, BCC, HCP, Diamond, SimpleCubic, BCT, and a
//! `Custom` escape hatch.

#![warn(missing_docs)]

mod crystal_structure;
mod lattice;
mod miller;
mod slip_system;

pub use crystal_structure::CrystalStructure;
pub use lattice::LatticeParameters;
pub use miller::MillerIndex;
pub use slip_system::{SlipFamily, SlipSystem, SlipSystemError};
