# tpt-mat-corrosion

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-corrosion.svg)](https://crates.io/crates/tpt-mat-corrosion)
[![Documentation](https://docs.rs/tpt-mat-corrosion/badge.svg)](https://docs.rs/tpt-mat-corrosion)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Environmental degradation: Butler–Volmer mixed-potential corrosion
> rate, Tafel polarization curves, Wagner parabolic oxidation scale
> growth.

`tpt-mat-corrosion` covers electrochemical degradation (consumed by
`tpt-medical` for implant selection and `tpt-transport` for
alloy/environment screens):

- `ElectrodeKinetics { exchange_current_density, tafel_slope, equilibrium_potential }`.
- `CorrosionModel { anode, cathode, electrolyte }`.
- `corrosion_rate(T, pH) -> CorrosionRate` — Butler–Volmer
  mixed-potential solve returning `current_density`, `penetration_rate`
  (mm / yr), `mass_loss_rate` (g / m² · day).
- `polarization_curve(potential_range) -> PolarizationCurve` — anodic /
  cathodic Tafel branches.
- **Oxidation module** (`oxidation` submodule):
  - Wagner parabolic growth (`ParabolicRateConstant`, `ScaleGrowth`,
    `DopingEffect`, `BreakawayCriterion`, `LinearBreakawayRate`,
    `OxidationModel`).

## Example

```rust
use tpt_mat_corrosion::{CorrosionModel, ElectrodeKinetics, corrosion_rate};

let fe = ElectrodeKinetics { i_0: 1e-7, beta: 0.05, e_eq: -0.44 };
let rate = corrosion_rate(&fe, &other, 298.0, 0.0); // mm/yr
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).