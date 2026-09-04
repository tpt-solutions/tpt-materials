# tpt-mat-creep

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-creep.svg)](https://crates.io/crates/tpt-mat-creep)
[![Documentation](https://docs.rs/tpt-mat-creep/badge.svg)](https://docs.rs/tpt-mat-creep)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Classical analytical creep models: Norton–Bailey power law,
> θ-projection (Wilshire–Burt), Monkman–Grant, Larson–Miller,
> Sherby–Dorn.

`tpt-mat-creep` packages the analytical steady-state / primary-creep
models used by `tpt-energy` and `tpt-transport` for alloy selection:

- `NortonBailey { a, n }` — `ε̇_ss = A σ^n exp(−Q / RT)`.
- `ThetaProjection { theta_1, theta_2, beta_1, beta_2 }` — Wilshire–Burt
  primary + secondary projection.
- `MonkmanGrant { log_c, m }` — `t_r = C / ε̇_ss^m`.
- `LarsonMiller { c }` — `LM = T (C + log t_r)`.
- `SherbyDorn { a, n }` — `ln t_r − Q/(RT) = f(σ)`.

## Example

```rust
use tpt_mat_creep::{NortonBailey, norton_strain_rate};

let law = NortonBailey { a: 1.5e-5, n: 5.0, q: 200_000.0 };
let edot = norton_strain_rate(&law, 100.0, 1073.0); // s⁻¹
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).