//! Crystallographic texture analysis.
//!
//! Phase 2 delivers:
//!
//! - [`TextureAnalyzer`]: a collection of grain orientations + weights
//!   (volume fractions).
//! - [`PoleFigure`]: equal-area or stereographic projection of a
//!   crystallographic pole (e.g. `{111}` for FCC, `{0001}` for HCP).
//! - [`OrientationDistributionFunction`]: kernel-density-estimate ODF
//!   `f(g)` on `SO(3)` (axis-angle parameterisation for now).
//! - [`taylor_factor`]: the Sachs/Taylor Taylor factor `M` averaged
//!   over a random texture (`M ≈ 3.06` for FCC, see Taylor, 1938).
//!
//! The module is `serde`-friendly so ODF / pole-figure data can be
//! persisted to JSON for downstream plotting.

#![warn(missing_docs)]

mod odf;
mod pole_figure;
mod taylor;
mod texture;

pub use odf::{OrientationDistributionFunction, OdfKernel};
pub use pole_figure::{PoleFigure, PoleFigureGrid, PoleKind};
pub use taylor::taylor_factor;
pub use texture::TextureAnalyzer;