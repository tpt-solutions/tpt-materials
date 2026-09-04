# tpt-mat-welding

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-welding.svg)](https://crates.io/crates/tpt-mat-welding)
[![Documentation](https://docs.rs/tpt-mat-welding/badge.svg)](https://docs.rs/tpt-mat-welding)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Welding heat-affected zone (HAZ) modelling: Rosenthal peak-temperature
> profile, CG-HAZ cooling rate, grain coarsening, Avrami + KM phase
> fractions across HAZ sub-zones.

`tpt-mat-welding` is the welding companion to
`tpt-mat-heat-treatment`; it provides the spatial HAZ discretisation that
the heat-treatment crate does not:

- `WeldModel { base_metal, filler_metal, process }`.
- `heat_affected_zone(heat_input) -> HazResult` — HAZ width + peak-
  temperature profile.
- `predict_haz_microstructure(cooling_rate) -> HazMicrostructure` —
  grain coarsening + transformation using `tpt-mat-phase-transform`
  (Avrami + KM).

## Example

```rust
use tpt_mat_welding::{WeldModel, WeldProcess, predict_haz_microstructure};

let weld = WeldModel::gmaw_default();
let haz = weld.heat_affected_zone(1500.0); // J / mm
let ms = predict_haz_microstructure(haz.cooling_rate(0.005));
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).