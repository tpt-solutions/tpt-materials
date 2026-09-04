# tpt-mat-phase-field

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-phase-field.svg)](https://crates.io/crates/tpt-mat-phase-field)
[![Documentation](https://docs.rs/tpt-mat-phase-field/badge.svg)](https://docs.rs/tpt-mat-phase-field)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Phase-field method: Allen-Cahn, Cahn-Hilliard, Kobayashi, multi-phase on
> regular grids.

`tpt-mat-phase-field` is the phase-field framework the spinodal,
dendritic-solidification, and phase-field fracture solvers all build on:

- `PhaseFieldSolver` — explicit forward-Euler driver with CFL stability
  limits.
- `PhaseFieldModel` enum — `AllenCahn`, `CahnHilliard`, `Kobayashi`,
  `MultiPhase`.
- `BulkEnergy` — `DoubleWell`, `Polynomial`, `RegularSolution`.
- `FreeEnergyFunctional` — central-difference evaluation in 2D.
- `step_allen_cahn`, `step_cahn_hilliard`, `step_kobayashi`,
  `step_multi_phase`.
- `PhaseFieldResult` — order parameter, concentration, temperature,
  free energy, interface area.

Verified: free energy monotonically decreases under Allen-Cahn step;
mean concentration conserved under Cahn–Hilliard step (to ~1 × 10⁻⁶).

## Example

```rust
use tpt_mat_phase_field::{PhaseFieldSolver, PhaseFieldModel, step_cahn_hilliard};

let solver = PhaseFieldSolver::new(Grid2D::new(128, 128, 1.0), PhaseFieldModel::CahnHilliard);
let mut c = initial_concentration();
for _ in 0..300 {
    c = step_cahn_hilliard(&solver, &c, 1e-3);
}
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).