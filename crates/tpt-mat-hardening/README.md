# tpt-mat-hardening

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-hardening.svg)](https://crates.io/crates/tpt-mat-hardening)
[![Documentation](https://docs.rs/tpt-mat-hardening/badge.svg)](https://docs.rs/tpt-mat-hardening)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Hardening laws for crystal plasticity: Voce, Power-Law, Kocks–Mecking,
> dislocation-density, combined / latent hardening.

`tpt-mat-hardening` provides the isotropic / kinematic / latent hardening
laws consumed by `tpt-mat-crystal-plasticity`:

- `Voce` — `τ_c = τ_0 + (τ_s − τ_0)(1 − e^{−γ/γ_c}) + θ_0 γ`.
- `PowerLaw` — `Δτ_c = h_0 (τ_s − τ_c) |Δγ| / τ_s`.
- `KocksMecking` — implicit exponential saturation.
- `DislocationDensity` — Kocks–Mecking evolution `dρ/dγ = k₁ √ρ − k₂ ρ`
  with Taylor stress `τ = α μ b √ρ`.
- `CombinedHardening` — additive sum of any two hardening laws.
- `LatentHardeningMatrix` — `h_{αβ} = h_0 (q for α ≠ β, 1 for α = β)`
  with the convention `q = 1.0` for coplanar slip, `q = 1.4` for
  non-coplanar (default).

All laws are trait-free `pub fn` calls plus a `Hardening` enum variant for
convenience.

## Example

```rust
use tpt_mat_hardening::{Voce, voce_update};

let law = Voce { tau_0: 30.0, tau_s: 80.0, gamma_c: 0.05, theta_0: 200.0 };
let tau_c = voce_update(&law, 0.0);     // initial CRSS
let tau_c = voce_update(&law, 0.10);    // after 0.10 shear
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).