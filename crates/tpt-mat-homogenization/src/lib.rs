//! Analytical micromechanical homogenization.
//!
//! This crate delivers the classical closed-form estimates that bound
//! and approximate the effective elastic stiffness `C_eff` of a
//! heterogeneous medium given the stiffnesses `C_r` and volume
//! fractions `f_r` of its `N` phases.
//!
//! # Bounds
//!
//! - [`voigt`] / [`reuss`]: first-order arithmetic and harmonic
//!   bounds on the effective stiffness, in Voigt (stiffness) and
//!   Reuss (compliance) averages respectively.
//! - [`voigt_reuss_average`]: the Hill average `M_VRH` = `½ (M_V + M_R)`
//!   commonly used for polycrystalline elastic moduli.
//! - [`hashin_shtrikman`]: tighter second-order bounds for an
//!   isotropic two-phase composite with spherical inclusions
//!   (Hashin & Shtrikman, 1962).  Generalises to `N` phases by
//!   sequential pairwise homogenisation.
//!
//! # Eshelby inclusion
//!
//! - [`eshelby_tensor`]: Eshelby's tensor `S` for a spherical
//!   inclusion in an isotropic matrix.  Relates the eigenstrain
//!   `ε*` inside the inclusion to the constrained strain `ε^c`
//!   via `ε^c = S : ε*`.
//! - [`dilute_strain_concentration`]: the dilute strain-concentration
//!   tensor `A = [I + S C_0^{-1} (C_1 - C_0)]^{-1}` for a single
//!   ellipsoidal inclusion (Mori–Tanaka, 1973 special case).
//!
//! # Conventions
//!
//! All 6x6 stiffness and compliance matrices use engineering-shear
//! Voigt order, identical to [`tpt_math_linalg_fixed::Vec6`].  The
//! implementations use the [`tpt_mat_crystal_plasticity::SymmetricFourthOrder`]
//! stiffness representation.

#![warn(missing_docs)]

mod bounds;
mod eshelby;
mod fft;
mod hashin_shtrikman;
mod voigt_reuss;

pub use bounds::{voigt_reuss_bounds, VoigtReussBounds};
pub use eshelby::{
    dilute_strain_concentration, eshelby_spherical, EshelbySpherical, StrainConcentrationTensor,
};
pub use fft::{fft2d, moulinec_suquet_2d};
pub use hashin_shtrikman::{
    hashin_shtrikman_k_g, hashin_shtrikman_spherical_bulk, hashin_shtrikman_spherical_shear,
    hashin_shtrikman_two_phase, HashinShtrikmanResult,
};
pub use voigt_reuss::{g_from_e_nu, k_from_e_nu, reuss, voigt, voigt_reuss_average};
