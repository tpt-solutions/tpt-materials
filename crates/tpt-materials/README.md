# tpt-materials

[![Crates.io](https://img.shields.io/crates/v/tpt-materials.svg)](https://crates.io/crates/tpt-materials)
[![Documentation](https://docs.rs/tpt-materials/badge.svg)](https://docs.rs/tpt-materials)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Umbrella / facade crate re-exporting every `tpt-mat-*` domain crate
> behind optional per-domain cargo features (`crystal-plasticity`,
> `phase-field`, …, `full`).

`tpt-materials` is the single entry-point crate that the spec §6 and
§13 import snippets resolve against:

```rust
use tpt_materials::crystallography::{CrystalStructure, slip_systems};
use tpt_materials::crystal_plasticity::{solve_increment_single_point, CrystalPlasticityModel};
use tpt_materials::energy_materials::{ActiveMaterial, BatteryChemistry};
```

## Features

| Feature | Pulls in |
|---|---|
| `crystallography` | `tpt-mat-crystallography` |
| `crystal-plasticity` | `tpt-mat-crystal-plasticity`, `tpt-mat-hardening` |
| `phase-field` | `tpt-mat-phase-field` |
| `phase-transform` | `tpt-mat-phase-transform` |
| `calphad` | `tpt-mat-calphad` |
| `homogenization` | `tpt-mat-homogenization` |
| `rve` | `tpt-mat-rve` |
| `damage` / `fatigue` / `fatigue-micro` / `creep` / `corrosion` | the matching crates |
| `battery` / `additive` / `welding` / `heat-treatment` | the matching crates |
| `database` / `machine-learning` | the matching crates |
| `fracture` / `thermal` / `dislocation` / `precipitation` | the matching crates |
| `hydrogen-embrittlement` / `polymer` / `hydrogel` / `hydrogen-storage` | the matching crates |
| `inverse` | `tpt-mat-inverse` |
| `wasm` | `tpt-mat-wasm` |
| `adapters` | `battery` |
| `full` | every feature above |

## Cargo features quick-start

```toml
[dependencies]
tpt-materials = { version = "0.1", features = ["crystal-plasticity", "phase-field"] }
```

```bash
cargo build -p tpt-materials --features full
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).