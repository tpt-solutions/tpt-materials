# tpt-mat-composite-micro

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-composite-micro.svg)](https://crates.io/crates/tpt-mat-composite-micro)
[![Documentation](https://docs.rs/tpt-mat-composite-micro/badge.svg)](https://docs.rs/tpt-mat-composite-micro)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Composite micromechanics: rule of mixtures, dilute estimate,
> Mori–Tanaka two-phase closed-form and N-phase iterative driver.

`tpt-mat-composite-micro` is the lightweight companion to
`tpt-mat-homogenization` and `tpt-mat-rve` focused on inclusion-based
composites:

- `rule_of_mixtures(c_matrix, c_inclusion, f) -> f64` — Voigt average.
- `dilute_estimate(c_matrix, c_inclusion, f, s) -> SymmetricFourthOrder` —
  non-interacting inclusions using the Eshelby strain-concentration tensor.
- `mori_tanaka(c_matrix, c_inclusion, f) -> SymmetricFourthOrder` —
  Benveniste (1987) closed-form two-phase MT; recovers matrix at `f = 0`
  and inclusion at `f = 1`.
- `mori_tanaka_iterative(phases) -> SymmetricFourthOrder` — N-phase
  iterative driver.

## Example

```rust
use tpt_mat_composite_micro::{mori_tanaka, SymmetricFourthOrder};

let c_al = SymmetricFourthOrder::isotropic(70.0, 0.33);
let c_sic = SymmetricFourthOrder::isotropic(450.0, 0.17);
let c_eff = mori_tanaka(&c_al, &c_sic, 0.20);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).