# tpt-mat-precipitation

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-precipitation.svg)](https://crates.io/crates/tpt-mat-precipitation)
[![Documentation](https://docs.rs/tpt-mat-precipitation/badge.svg)](https://docs.rs/tpt-mat-precipitation)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Precipitation hardening kinetics — classical nucleation theory
> (Turnbull–Fisher), Kampmann–Wagner numerical (KWN) size-class
> solver, LSW coarsening, Orowan / shearing strengthening.

`tpt-mat-precipitation` covers the precipitation pathway that
turns an alloy CALPHAD description into a yield-strength increment:

- `ClassicalNucleation` — `I = I_0 exp(−ΔG*/kT) exp(−Q/kT)`; `ΔG*`
  from interfacial energy `γ` and driving force `Δg_v` (driving
  force from `tpt-mat-calphad`).
- `KwnModel` — Kampmann–Wagner–Numerical: discretised size classes,
  coupled nucleation + growth (`dR/dt`) + capillarity.
- `LswCoarsing` — `R̄³ − R̄_0³ = K t`, `K` from `γ`, `D`, `c_eq`,
  `V_m`.
- `PrecipitateState` — number density, mean radius, volume
  fraction, matrix supersaturation vs time.
- `strengthening_increment` — Orowan bypass + shearing → `Δσ_y`,
  hand-off to `tpt-mat-hardening`.

Verified: KWN conserves solute mass; late-stage slope → LSW `t^{1/3}`.

## Example

```rust
use tpt_mat_precipitation::{KwnModel, ClassicalNucleation, strengthening_increment};

let state = KwnModel::default().run(temperatures, times);
let dsigma = strengthening_increment(&state, &alloy);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).