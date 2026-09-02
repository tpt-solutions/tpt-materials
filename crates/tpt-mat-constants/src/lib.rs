//! Physical constants (CODATA 2018) and a small periodic table.
//!
//! Values follow CODATA 2018 recommended values and the IUPAC 2021
//! standard atomic weights.  All quantities use SI units (kg, m, s, K,
//! mol, J).

#![warn(missing_docs)]

mod constants;
mod periodic_table;

pub use constants::PhysicalConstants;
pub use periodic_table::{AtomicData, PeriodicTable};
