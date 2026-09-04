# tpt-mat-additive

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-additive.svg)](https://crates.io/crates/tpt-mat-additive)
[![Documentation](https://docs.rs/tpt-mat-additive/badge.svg)](https://docs.rs/tpt-mat-additive)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Additive-manufacturing process models: LPBF, DED, EBM thermal
> histories (Rosenthal / Gaussian), Hunt CET map, residual-stress field.

`tpt-mat-additive` packages the analytic AM thermal-history and
microstructure models used by the additive-manufacturing demo:

- `AmProcess` enum — `LaserPowderBedFusion { laser_power, scan_speed, hatch_spacing, layer_thickness }`,
  `DirectedEnergyDeposition`, `ElectronBeamMelting`.
- `thermal_history(location) -> ThermalHistory` — Rosenthal
  moving-point-source or Gaussian-beam analytic solution.
- `predict_microstructure(&ThermalHistory) -> PredictedMicrostructure` —
  grain size, phase fractions, texture, porosity; columnar / equiaxed
  from the Hunt G–R solidification map.
- `residual_stress(&ThermalHistory) -> ResidualStressField` —
  thermal-contraction eigenstrain form.

Verified: higher cooling rate → finer predicted grain size (Hall–Petch
trend).

## Example

```rust
use tpt_mat_additive::{AmProcess, predict_micro};

let process = AmProcess::LaserPowderBedFusion {
    laser_power: 200.0, scan_speed: 1.2, hatch_spacing: 0.10, layer_thickness: 0.03,
};
let th = process.thermal_history([0.0, 0.0, 0.0]);
let ms = predict_microstructure(&th);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).