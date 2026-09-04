# tpt-mat-diffusion

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-diffusion.svg)](https://crates.io/crates/tpt-mat-diffusion)
[![Documentation](https://docs.rs/tpt-mat-diffusion/badge.svg)](https://docs.rs/tpt-mat-diffusion)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Single- and multi-component diffusion on regular grids (Fick's laws,
> Arrhenius diffusivity, grain-boundary diffusion).

`tpt-mat-diffusion` solves Fickian diffusion on `tpt-science` grids:

- `ArrheniusDiffusivity { D_0, Q, T }` — `D = D_0 exp(−Q / (R T))`.
- `DiffusionSolver` — single-component forward-Euler Fick's second law
  with Neumann zero-flux boundaries; mass conserved to ~1 × 10⁻⁶.
- `MultiComponentDiffusionSolver` — independent species, per-species CFL
  stability limit.
- `GrainBoundaryDiffusion` — Fisher Regime-A effective diffusivity
  `D_eff = D_l + (w / L) (δ D_gb / D_l)^{1/2}` form (Fisher 1951).

The driving example is `examples/diffusion-carbon-steel`: 1-D
carburisation of a steel slab; `D(T)` is Arrhenius-driven and case
depth grows from 0.01 mm at t = 0 to 0.62 mm at t = 60 s.

## Example

```rust
use tpt_mat_diffusion::{ArrheniusDiffusivity, DiffusionSolver};

let d = ArrheniusDiffusivity { d_0: 2.3e-5, q: 142_000.0 }.at(1200.0);
let solver = DiffusionSolver::new(grid, d);
let profile = solver.solve(initial_carbon_profile, 1.0, 60);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).