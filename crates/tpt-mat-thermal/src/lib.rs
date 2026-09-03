//! Effective thermal & transport-property homogenization.
//!
//! Symmetric-positive-definite transport tensors (thermal
//! conductivity `k`, electrical conductivity `σ`, ionic diffusivity
//! `D`) obey the same Voigt / Reuss / Hashin–Shtrikman /
//! self-consistent / Maxwell–Garnett bounds as the elastic
//! stiffness tensor.  This crate reuses those bounds for the
//! scalar case and adds:
//!
//! - [`effective_conductivity`] — series/parallel/HS/self-consistent
//!   bounds for scalar conductivities.
//! - [`effective_cte`] — Turner / Kerner / Rosen–Hashin bounds for
//!   the linear coefficient of thermal expansion.
//! - [`effective_specific_heat`] — mass-weighted rule of mixtures.
//! - [`effective_diffusivity`] — tortuosity-corrected Bruggeman
//!   relation for porous media.
//! - [`interface_thermal_resistance`] — Kapitza-resistance correction
//!   for laminate composites.
//!
//! # References
//!
//! - Hashin, Z., & Shtrikman, S. (1962).  "A variational approach
//!   to the theory of the effective magnetic permeability of
//!   multiphase materials."  J. Appl. Phys. 33, 3125–3131.
//! - Turner, P. S. (1946).  "Thermal-expansion stresses in
//!   reinforced plastics."  J. Res. NBS 37, 239–247.
//! - Rosen, B. W., & Hashin, Z. (1970).  "Effective thermal
//!   expansion coefficients and specific heats of particulate
//!   composites."  Int. J. Eng. Sci. 8, 157–173.
//! - Bruggeman, D. A. G. (1935).  "Berechnung verschiedener
//!   physikalischer Konstanten von heterogenen Substanzen."  Ann.
//! Phys. 24, 636-664.

#![warn(missing_docs)]

mod conductivity;
mod cte;
mod diffusivity;
mod resistance;
mod specific_heat;

pub use conductivity::{effective_conductivity, hashin_shtrikman_k, ConductivityBound};
pub use cte::{effective_cte, CteBound};
pub use diffusivity::{bruggeman_diffusivity, effective_diffusivity};
pub use resistance::kapitza_correction;
pub use specific_heat::effective_specific_heat;
