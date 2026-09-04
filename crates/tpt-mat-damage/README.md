# tpt-mat-damage

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-damage.svg)](https://crates.io/crates/tpt-mat-damage)
[![Documentation](https://docs.rs/tpt-mat-damage/badge.svg)](https://docs.rs/tpt-mat-damage)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Continuum damage mechanics (Kachanov / Lemaitre / Chaboche / Miner)
> plus Gurson–Tvergaard–Needleman porous plasticity.

`tpt-mat-damage` groups two distinct schools of continuum damage:

- **Classical CDM**
  - Kachanov effective-stress `σ̃ = σ / (1 − D)`.
  - Lemaitre ductile damage evolution.
  - Kachanov creep-damage coupling.
  - Miner linear damage accumulation.
  - Chaboche non-linear kinematic / isotropic hardening placeholder.
- **GTN porous plasticity** (added Phase 6)
  - `GursonTvergaardNeedleman { f_0, f_c, f_f, q_1, q_2, q_3, … }`.
  - `yield_function(σ, f, σ_y)` — `Φ = (σ_eq / σ_y)² + 2 q₁ f cosh(3 q₂ σ_m / 2 σ_y) − (1 + q₃ f²)`.
  - `update_porosity` — void growth `(1 − f) dε^p_kk` + strain-controlled
    nucleation (Chu–Needleman).
  - Verification: GTN yield → von Mises as `f → 0`.

## Example

```rust
use tpt_mat_damage::{Lemaitre, kachanov_effective_stress};
let sigma_tilde = kachanov_effective_stress(100.0, 0.3);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).