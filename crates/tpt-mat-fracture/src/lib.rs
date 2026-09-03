//! Fracture mechanics.
//!
//! This crate bundles the classical analytical fracture models with
//! a regularised phase-field fracture functional.  The phase-field
//! solve itself reuses the [`tpt_science`] Laplacian helpers via
//! the integration layer; the local dissipation functional,
//! AT1/AT2 degradation functions, and the cohesive-zone traction
//! separation laws are provided here.
//!
//! ## Models
//!
//! - [`StressIntensityFactor`] — `K_I` for a centre crack in an
//!   infinite plate (`K = σ √(π a)`) and an edge crack
//!   (`K = 1.12 σ √(π a)`).
//! - [`EnergyReleaseRate`] — Irwin relation `G = K_I² / E'` (plane
//!   strain `E' = E / (1 − ν²)`, plane stress `E' = E`).
//! - [`CohesiveZoneModel`] — bilinear traction–separation law with
//!   mixed-mode Benzeggagh–Kenane decomposition.
//! - [`PhaseFieldFracture`] — AT1 and AT2 degradation functions
//!   `g(d)` and crack-surface dissipation functional
//!   `∫ (G_c / c_w) [w(d)/l_0 + l_0 |∇d|²] dV`.
//! - [`fracture_toughness_transition`] — ASTM E1921 master curve
//!   reference temperature `T_0`.
//!
//! # References
//!
//! - Griffith, A. A. (1921).  "The phenomena of rupture and flow
//!   in solids."  Phil. Trans. R. Soc. A 221, 163–198.
//! - Irwin, G. R. (1957).  "Analysis of stresses and strains near
//!   the end of a crack traversing a plate."  J. Appl. Mech. 24,
//!   361–364.
//! - Benzeggagh, M. L., & Kenane, M. (1996).  "Measurement of
//!   mixed-mode delamination fracture toughness of unidirectional
//!   glass/epoxy composites with mixed-mode bending apparatus."
//!   Compos. Sci. Technol. 56, 439–449.
//! - Bourdin, B., Francfort, G., & Marigo, J.-J. (2000).  "Numerical
//!   experiments in revisited brittle fracture."  JMPS 48, 797–826.
//! - Miehe, C., Welschinger, F., & Hofacker, M. (2010).  "Thermally
//!   consistent phase-field models of fracture: variational
//!   principle and multi-field FE implementations."  IJNME 83,
//!   1273-1311.

#![warn(missing_docs)]

mod cohesive;
mod energy_release_rate;
mod phase_field;
mod sif;
mod transition;

pub use cohesive::{CohesiveZoneModel, MixedModeDecomposition, TractionSeparation};
pub use energy_release_rate::{
    energy_release_rate_irwin, energy_release_rate_j_integral, IrwinModulus,
};
pub use phase_field::{at1_degradation, at2_degradation, dissipation_density, PhaseFieldFracture};
pub use sif::{
    k_i_centre_crack, k_i_edge_crack, k_ii_centre_crack, k_iii_centre_crack, StressIntensityFactor,
};
pub use transition::{fracture_toughness_transition, master_curve_k_jc, MasterCurveParams};
