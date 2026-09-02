//! Core material data model.
//!
//! Phase 1 provides:
//! - [`MaterialMicrostructure`]
//! - [`Phase`]
//! - [`Grain`]
//! - [`Composition`] / [`CompositionBasis`]
//! - [`CrystalOrientation`] / [`OrientationRepresentation`]
//!
//! These types are the data backbone for every solver crate in later
//! phases.

#![warn(missing_docs)]

mod composition;
mod grain;
mod microstructure;
mod orientation;
mod phase;

pub use composition::{Composition, CompositionBasis, CompositionError};
pub use grain::{Grain, GrainId};
pub use microstructure::MaterialMicrostructure;
pub use orientation::{CrystalOrientation, OrientationConversionError, OrientationRepresentation};
pub use phase::{Phase, PhaseId};
