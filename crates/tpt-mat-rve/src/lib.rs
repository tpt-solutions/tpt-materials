//! Representative Volume Element (RVE) framework.
//!
//! An RVE in micromechanics is a microstructural sample large enough
//! to be statistically representative of the macroscopic behaviour of
//! the heterogeneous medium.  This crate provides a *minimal* data
//! model — sufficient to drive the analytical homogenizers in
//! [`tpt_mat_homogenization`] and the Bishop–Hill Taylor-factor solver
//! in [`tpt_mat_rve::bishop_hill`] — without committing to a specific
//! FEM solver.  A `Rve` carries:
//!
//! - `n_grains` grains, each with a crystallographic orientation
//!   (axis-angle rotation matrix in the sample frame),
//! - a per-grain volume fraction (`w_i`), and
//! - a per-grain stiffness tensor `C_i`.
//!
//! The [`HomogenizationScheme`] enum lets the caller choose between
//! `Voigt` (iso-strain), `Reuss` (iso-stress), and the proper
//! [`SelfConsistent`](SelfConsistent) one-site mean-field scheme
//! (Kröner, 1958; Budiansky & Wu, 1962).
//!
//! ## Bishop–Hill Taylor factor
//!
//! [`bishop_hill_taylor_factor`] solves the small linear
//! complementarity problem
//!
//! `min_γ  σ : ε`   s.t.   `σ = C : ε^p`,   `|τ^α| ≤ τ_c^α`,
//!
//! where the plasticity is fully accommodated by exactly 5 active slip
//! systems in each grain; this is the Taylor/Bishop–Hill assumption.
//! The Taylor factor `M = σ_eq / τ_c` is the ratio of the macroscopic
//! stress intensity to the critical resolved shear stress.  For a
//! random FCC polycrystal `M ≈ 3.06` (Taylor, 1938).

#![warn(missing_docs)]

mod bishop_hill;
mod lemke;
mod rve;

pub use bishop_hill::{
    bishop_hill_lemke, bishop_hill_lemke_with_slips, bishop_hill_taylor_factor,
    bishop_hill_taylor_factor_axis, BishopHillResult,
};
pub use lemke::{lemke_solve, LemkeResult};
pub use rve::{HomogenizationScheme, Rve, RveGrain, RveStats, SimpleHomogenizer};

pub use tpt_mat_homogenization as homogenization;
