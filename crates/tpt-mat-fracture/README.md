# tpt-mat-fracture

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-fracture.svg)](https://crates.io/crates/tpt-mat-fracture)
[![Documentation](https://docs.rs/tpt-mat-fracture/badge.svg)](https://docs.rs/tpt-mat-fracture)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Fracture mechanics — LEFM stress intensity factors, Irwin energy
> release rate, J-integral, cohesive-zone models, phase-field
> fracture (AT1 / AT2), DBTT master curve.

`tpt-mat-fracture` covers the three classical fracture modelling
schools used by every structural-life consumer:

- `StressIntensityFactor` — `K_I`, `K_II`, `K_III` with geometry
  factors (centre / penny-shaped, edge crack, Mode II / Mode III),
  `k_from_load(geometry, stress, crack_length)`.
- `EnergyReleaseRate` — `G`, Irwin `G = K² / E'`, J-integral
  (domain-integral form on a `tpt-science` grid).
- `CohesiveZoneModel` — bilinear / exponential traction–separation
  (`t_0`, `δ_c`, `G_c`); mixed-mode Benzeggagh–Kenane.
- `PhaseFieldFracture` — Griffith / AT1 / AT2 regularised
  functionals (`κ`, `G_c`, `l_0`); staggered solve reusing the
  `tpt-mat-phase-field` infrastructure.
- `fracture_toughness_transition` — DBTT / master-curve
  (ASTM E1921) helper.

Verified: phase-field fracture recovers the Griffith load for a 1-D
bar; Irwin `K → G` consistency.

## Example

```rust
use tpt_mat_fracture::{k_from_load, Geometry};

let k = k_from_load(Geometry::CentreCrack { a: 0.01 }, 200.0e6, 0.01);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).