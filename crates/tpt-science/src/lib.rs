//! `tpt-science` — substrate crate for `tpt-materials`.
//!
//! Houses small, dependency-light numerical helpers shared by every
//! solver crate in the workspace:
//   - [`grid`]: regular Cartesian grids + finite-difference / spectral
//!     Laplacian operators for Allen-Cahn, Cahn-Hilliard, Fick's
//!     second law, and thermal diffusion.
//!
//! Larger dense / sparse linear systems live in the cross-repo
//! `tpt-fem` substrate (mesh, assembly, Newton-Raphson); this crate
//! intentionally ships only the small, self-contained helpers needed
//! by Phase 3 (phase-field) and Phase 4 (diffusion / CALPHAD).

#![warn(missing_docs)]

pub mod grid;

pub use grid::{Grid1D, Grid2D, Grid3D};
