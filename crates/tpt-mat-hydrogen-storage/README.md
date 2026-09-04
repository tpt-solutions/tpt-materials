# tpt-mat-hydrogen-storage

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-hydrogen-storage.svg)](https://crates.io/crates/tpt-mat-hydrogen-storage)
[![Documentation](https://docs.rs/tpt-mat-hydrogen-storage/badge.svg)](https://docs.rs/tpt-mat-hydrogen-storage)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Hydrogen storage materials: metal hydride, chemical hydride, and
> porous-adsorbent PCT isotherms with van 't Hoff temperature
> dependence.

`tpt-mat-hydrogen-storage` is the small soft-matter / energy crate
covering sorption and absorption of hydrogen for storage system
sizing:

- `HydrideType` enum — `MetalHydride { alloy }`,
  `ChemicalHydride { compound }`, `PorousMaterial { surface_area }`.
- `HydrogenStorageMaterial { hydride_type, storage_capacity_wt_pct, absorption_kinetics }`.
- `pct_isotherm(temperature) -> Vec<(f64, f64)>` — pressure–composition–
  temperature curve: plateau + van 't Hoff temperature dependence.

Driving example: `examples/metal-hydride-pct` produces the LaNi₅ PCT
isotherm family at 298 / 348 / 398 K.

## Example

```rust
use tpt_mat_hydrogen_storage::{HydrogenStorageMaterial, HydrideType};

let lani5 = HydrogenStorageMaterial::lani5_default();
let isotherm = lani5.pct_isotherm(298.0);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).