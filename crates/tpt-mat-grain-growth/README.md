# tpt-mat-grain-growth

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-grain-growth.svg)](https://crates.io/crates/tpt-mat-grain-growth)
[![Documentation](https://docs.rs/tpt-mat-grain-growth/badge.svg)](https://docs.rs/tpt-mat-grain-growth)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Grain growth, interfaces, recrystallisation, stereology.

`tpt-mat-grain-growth` is the microstructural-evolution companion to
`tpt-mat-phase-field` with a particular focus on grain-boundary physics:

- `GrainBoundaryMobility` — Read–Shockley low-angle mobility + HAGB
  multiplier.
- `GrainGrowthSolver` — curvature-driven front-tracking driver.
- `grain_size_distribution` → histogram + `GrainSizeStats`.
- `interfaces` (module): GB energy (Read–Shockley), GBCD / CSL
  fractions, triple-junction Herring balance, Langmuir–McLean
  segregation.
- `recrystallization` (module): static JMAK RX, Zener–Hollomon +
  Sellars–Tegart, dynamic RX (DrxKinetics, Cahn–Hagel form).
- `stereology` (module): ASTM E112 linear intercept, area / volume
  fraction + counting uncertainty, Saltykov 3-D reconstruction.

## Example

```rust
use tpt_mat_grain_growth::{GrainGrowthSolver, GrainBoundaryMobility, ReadShockley};

let mob = GrainBoundaryMobility::ReadShockley { hagb: 1.0e-6, theta_c: 15.0_f64.to_radians() };
let solver = GrainGrowthSolver::new(mob);
let final_field = solver.run(initial_grain_ids, 1000);
let stats = grain_size_distribution(&final_field);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).