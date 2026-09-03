//! Hydrogen storage materials.
//!
//! - [`HydrideType`] — metal hydride, chemical hydride, porous
//!   adsorbent.
//! - [`HydrogenStorageMaterial`] — material record.
//! - [`pct_isotherm`] — pressure–composition–temperature
//!   curve with a plateau and van 't Hoff temperature
//!   dependence.

#![warn(missing_docs)]

mod material;
mod pct;

pub use material::{HydrideType, HydrogenStorageMaterial, PorousAdsorbent};
pub use pct::{pct_isotherm, van_t_hoff_pressure, PctParams, PctPoint};