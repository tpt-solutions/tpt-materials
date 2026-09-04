# tpt-science

[![Crates.io](https://img.shields.io/crates/v/tpt-science.svg)](https://crates.io/crates/tpt-science)
[![Documentation](https://docs.rs/tpt-science/badge.svg)](https://docs.rs/tpt-science)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Substrate crate for tpt-materials: regular-grid finite-difference / spectral
> PDE helpers (Allen-Cahn, Cahn-Hilliard, Fick's laws, electrochemistry).

`tpt-science` provides the regular-grid infrastructure used by every
phase-field, diffusion, and transport solver in the workspace:

- `Grid1D`, `Grid2D`, `Grid3D` — owned regular grids with `dx`, `shape`,
  `lini`, range-checked indexing.
- `BoundaryCondition` — `NeumannZeroFlux`, `Periodic`, `Dirichlet` applied
  uniformly or per-axis.
- `compute_laplacian` — second-order central differences with the chosen
  boundary treatment, on 1D / 2D / 3D grids.
- `compute_biharmonic` — `∇⁴` operator for Cahn–Hilliard and
  phase-field fracture regularised functionals.
- `compute_gradient` — central-difference vector field on a regular grid.

This crate has no materials-science semantics — it is purely a numerical
  substrate on top of `tpt-math-linalg-fixed`. Keeping it separate lets the
  domain crates (`tpt-mat-phase-field`, `tpt-mat-diffusion`,
  `tpt-mat-hydrogen-embrittlement`, `tpt-mat-fracture`) share one
  well-tested Laplacian implementation.

## Example

```rust
use tpt_science::{Grid2D, compute_laplacian, BoundaryCondition};

let grid = Grid2D::new(64, 64, 1.0);
let mut field = vec![0.0; grid.shape().0 * grid.shape().1];
field[0] = 1.0; // point source
let lap = compute_laplacian(&field, &grid, &BoundaryCondition::NeumannZeroFlux);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).