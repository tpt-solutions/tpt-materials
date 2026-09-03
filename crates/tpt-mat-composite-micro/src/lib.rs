//! Micromechanics of composite microstructures.
//!
//! This crate provides mean-field homogenization schemes for
//! multi-phase composites that go beyond the simple Voigt / Reuss
//! bounds:
//!
//! - [`mori_tanaka`]: the **Mori–Tanaka (1973)** scheme, in which the
//!   inclusion phase is embedded in the *average matrix* seen by an
//!   inclusion.  Gives an estimate for the effective stiffness
//!   `C_eff` that matches Hashin–Shtrikman for spherical inclusions
//!   and is widely used for particle- and short-fiber-reinforced
//!   composites.
//! - [`dilute_estimate`]: the **dilute (Maxwell, 1873)** scheme, in
//!   which inclusions are treated as non-interacting (valid only for
//!   `f ≪ 1`).
//! - [`rule_of_mixtures`]: the engineering rule-of-mixtures
//!   weighted by volume fraction.
//!
//! All three accept `(C_r, f_r)` slices and return a 6×6 effective
//! stiffness `C_eff` in Voigt form.
//!
//! # References
//!
//! - Mori, T. & Tanaka, K. (1973).  "Average stress in matrix and
//!   average elastic energy of materials with misfitting
//!   inclusions."  Acta Metall. 21, 571–574.
//! - Benveniste, Y. (1987).  "A new approach to the application of
//!   Mori–Tanaka's theory in composite materials."  Mech. Mater. 6,
//!   147–157.

#![warn(missing_docs)]

mod dilute;
mod mori_tanaka;
mod rule_of_mixtures;

pub use dilute::{dilute_estimate, eshelby_for_isotropic_sphere};
pub use mori_tanaka::{mori_tanaka, mori_tanaka_iterative, MoriTanakaResult};
pub use rule_of_mixtures::rule_of_mixtures;
