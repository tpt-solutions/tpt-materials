# tpt-mat-crystal-plasticity

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-crystal-plasticity.svg)](https://crates.io/crates/tpt-mat-crystal-plasticity)
[![Documentation](https://docs.rs/tpt-mat-crystal-plasticity/badge.svg)](https://docs.rs/tpt-mat-crystal-plasticity)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Single-crystal and polycrystal plasticity solver: rate-sensitive
> flow rule, radial-return integration, FIP export.

`tpt-mat-crystal-plasticity` is the constitutive heart of the workspace.
It provides everything needed for a single-crystal material point driver
and (optionally) a CP-FEM coupling:

- `CrystalPlasticityModel { crystal, slip_systems, hardening, rate,
  elastic }` — the complete material description.
- `RateSensitivity { gamma_dot_0, n }` — power-law viscoplastic flow
  `γ̇ = γ̇_0 (τ / τ_c)^{1/n} sign(τ)`.
- `power_law_slip_rate(tau, tau_c, rate) -> f64`.
- `resolved_shear_stresses(stress, slip_systems) -> Vec<f64>`.
- `viscoplastic_velocity_gradient(gamma_dots, slip_systems) -> Mat3` —
  `L^p = Σ_α γ̇^α s^α ⊗ n^α`.
- `solve_increment_single_point(model, state, dt, delta_L) -> PlasticIncrement`
  — radial-return with hardening update.
- `CpFemSolver { BoundaryConditions, LoadStep, CpFemResult }` — single-
  point integration glue.  Full Newton–Raphson FEM assembly is wired
  through optional `tpt-fem` handles (off by default; the integration
  test requires `tpt-fem` from another repo and is deferred).

Verified: under uniaxial tension along `x` an FCC single crystal
activates 5 / 12 slip systems with the correct Schmid ranking (see
`examples/fcc-single-crystal-tension`).

## Example

```rust
use tpt_mat_crystal_plasticity::{CrystalPlasticityModel, solve_increment_single_point};

let mp = CrystalPlasticityModel::fcc_aluminium_default();
let dp = solve_increment_single_point(&mp, &state, 1e-4, &delta_L);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).