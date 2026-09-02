//! Fixed-size linear algebra types for `tpt-materials`.
//!
//! This crate deliberately ships **no** BLAS / LAPACK / `nalgebra`
//! dependency. All types are `Copy`, `f64`-backed, `serde`-friendly, and
//! optimised for the 3x3 / 6x6 problems that dominate crystal plasticity,
//! phase-field, and homogenization.
//!
//! Larger dense or sparse linear systems belong in the downstream
//! `tpt-fem` and `tpt-science` substrate crates.

#![deny(missing_docs)]

mod mat3;
mod sym_mat3;
mod vec3;
mod vec6;

pub use mat3::Mat3;
pub use sym_mat3::SymMat3;
pub use vec3::Vec3;
pub use vec6::Vec6;
