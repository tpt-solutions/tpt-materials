# tpt-mat-thermal

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-thermal.svg)](https://crates.io/crates/tpt-mat-thermal)
[![Documentation](https://docs.rs/tpt-mat-thermal/badge.svg)](https://docs.rs/tpt-mat-thermal)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Effective thermal & transport-property homogenization — conductivity,
> CTE, specific heat, diffusivity, Kapitza interface resistance.

`tpt-mat-thermal` is the transport-property companion to
`tpt-mat-homogenization` (which is purely elastic). The same
symmetric-positive tensor machinery is generalised to any 2nd-order
transport property:

- `effective_conductivity` — series / parallel / VRH / Hashin–
  Shtrikman / self-consistent / Maxwell–Garnett (2-phase and N-phase).
- `effective_cte` — Turner, Kerner, Rosen–Hashin bounds for composite
  thermal expansion.
- `effective_specific_heat` — mass-weighted rule of mixtures.
- `effective_diffusivity` — tortuosity / Bruggeman for porous &
  multiphase media (shared math with electrical conductivity).
- `interface_thermal_resistance` — Kapitza-resistance correction.

Verified: conductivity HS bounds enclose the self-consistent estimate
and collapse at zero contrast.

## Example

```rust
use tpt_mat_thermal::{effective_conductivity, Scheme};

let k = effective_conductivity(&[k_matrix, k_inclusion], &[1.0 - f, f], Scheme::SelfConsistent);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).