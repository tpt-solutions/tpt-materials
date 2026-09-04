# tpt-mat-battery

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-battery.svg)](https://crates.io/crates/tpt-mat-battery)
[![Documentation](https://docs.rs/tpt-mat-battery/badge.svg)](https://docs.rs/tpt-mat-battery)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Li-ion battery active-material degradation: 1-D radial Li diffusion +
> diffusion-induced stress, SEI growth, particle cracking, lithium
> plating, transition-metal dissolution.

`tpt-mat-battery` is the source of the `DegradationCurve` artifact
consumed by `tpt-energy`:

- `ActiveMaterial { chemistry, particle_radius, diffusion_coefficient, partial_molar_volume }`.
- `BatteryChemistry` — NMC811, NMC622, LFP, NCA, graphite, silicon, …
- `DegradationMechanism` enum — `SeiGrowth { rate_constant, activation_energy }`,
  `ParticleCracking { critical_stress }`,
  `LithiumPlating { plating_potential }`,
  `TransitionMetalDissolution { dissolution_rate }`.
- `simulate_diffusion_stress(c_rate, num_cycles) -> DiffusionStressResult`
  — radial Li diffusion coupled with diffusion-induced stress in a
  spherical particle (reuses `tpt-mat-diffusion`).
- `capacity_fade_curve(cycles, temperature) -> DegradationCurve` —
  returns `cycles`, `capacity_retention`, `resistance_growth`.

Verified: capacity retention monotonically decreasing; √t SEI-limited
fade at low C-rate.

## Example

```rust
use tpt_mat_battery::{ActiveMaterial, BatteryChemistry};

let nmc = ActiveMaterial::default(BatteryChemistry::NMC811);
let curve = nmc.capacity_fade_curve(10_000, 298.0);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).