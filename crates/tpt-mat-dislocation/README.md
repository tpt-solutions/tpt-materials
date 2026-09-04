# tpt-mat-dislocation

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-dislocation.svg)](https://crates.io/crates/tpt-mat-dislocation)
[![Documentation](https://docs.rs/tpt-mat-dislocation/badge.svg)](https://docs.rs/tpt-mat-dislocation)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Dislocation-density-based hardening — per-system SSD / GND / forest
> state, Kocks–Mecking evolution, Taylor stress, Armstrong–Frederick
> back-stress, Nye-tensor GND reconstruction.

`tpt-mat-dislocation` provides the physically-based alternative to
phenomenological Voce hardening:

- `DislocationDensityState` — per-slip-system `ρ_SSD`, `ρ_GND`,
  forest density.
- `KocksMeckingEvolution` — `dρ/dγ = k_1 √ρ − k_2 ρ` (storage vs
  dynamic recovery); Taylor stress `τ = α μ b √ρ`.
- `back_stress` — Armstrong–Frederick kinematic term from GND
  gradients.
- `gnd_from_curvature` — Nye tensor → `ρ_GND` from a
  lattice-curvature field.

Verified: single-slip response reproduces Voce-like saturation;
`ρ` stays non-negative under the discrete Kocks–Mecking update.

## Example

```rust
use tpt_mat_dislocation::{DislocationDensityState, KocksMeckingEvolution};

let mut rho = DislocationDensityState::initial(1e10, 12);
let params = KocksMeckingEvolution { k_1: 1e8, k_2: 5.0 };
rho.advance(0.01, &params);
let tau = rho.taylor_stress(70e9, 2.5e-10, 0.3);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).