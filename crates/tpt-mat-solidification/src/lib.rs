//! Phase-field dendritic solidification.
//!
//! Wraps `tpt-mat-phase-field::PhaseFieldSolver` with anisotropic
//! interface stiffness for dendrite growth (Kobayashi, 1993), plus
//! post-processing to extract tip velocity and secondary-arm spacing.
//!
//! The anisotropy is the classical 4-fold cubic form
//! `a(n) = α (1 + ε cos(4 θ))` where `n` is the interface normal
//! direction in 2D.

#![warn(missing_docs)]

mod anisotropy;
mod solver;

pub use anisotropy::{AnisotropyMode, AnisotropyModel};
pub use solver::{SolidificationResult, SolidificationSolver};
