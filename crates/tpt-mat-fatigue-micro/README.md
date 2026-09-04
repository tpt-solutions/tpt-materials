# tpt-mat-fatigue-micro

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-fatigue-micro.svg)](https://crates.io/crates/tpt-mat-fatigue-micro)
[![Documentation](https://docs.rs/tpt-mat-fatigue-micro/badge.svg)](https://docs.rs/tpt-mat-fatigue-micro)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Microstructural fatigue indicator parameters (FIPs) and
> crack-initiation prediction: Findley, Fatemi–Socie,
> Smith–Watson–Topper, Tanaka–Mura / crystallographic slip.

`tpt-mat-fatigue-micro` adds the microstructure-aware fatigue pathway to
the classical `tpt-mat-fatigue` crate:

- `MicrostructuralFatigue { rve, criterion }`.
- `FatigueCriterion` enum — `Findley`, `FatemiSocie`,
  `SmithWatsonTopper`, `CrystallographicSlip { critical_accumulated_shear }`.
- `fatigue_indicator_parameter(&CpFemResult) -> Vec<f64>` — FIP field
  per grain.
- `predict_crack_initiation(&[LoadStep]) -> CrackInitiationResult` —
  `cycles_to_initiation`, `critical_grain`, `critical_location`, FIP
  field.
- Cycle-by-cycle CP-FEM driver that accumulates plastic slip at grain
  boundaries.

Driving example: `examples/fatigue-crack-initiation` runs a 5-grain
polycrystal Findley scan, identifies grain 2 as the critical grain,
and returns `N_i ≈ 12` cycles to initiation.

## Example

```rust
use tpt_mat_fatigue_micro::{MicrostructuralFatigue, FatigueCriterion, predict_crack_initiation};

let criterion = FatigueCriterion::Findley { k: 1.0 };
let fip = fatigue_indicator_parameter(&cp_result);
let init = predict_crack_initiation(&load_steps);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).