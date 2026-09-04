# tpt-mat-inverse

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-inverse.svg)](https://crates.io/crates/tpt-mat-inverse)
[![Documentation](https://docs.rs/tpt-mat-inverse/badge.svg)](https://docs.rs/tpt-mat-inverse)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Inverse-problem calibration helpers — closed-form fits for the
> standard constitutive laws used across the workspace, and a
> generic Levenberg–Marquardt driver for everything else.

`tpt-mat-inverse` is the small crate that lets every other `tpt-mat-*`
crate expose its model as a calibration target. It bundles:

- **Closed-form fits** (preferred — numerically robust, no iteration
  needed):
  - `fit_arrhenius(t, D) -> { D_0, Q }` — Arrhenius diffusivity fit
    (used by `tpt-mat-diffusion`).
  - `fit_norton_bailey(sigma, edot) -> { A, n, Q }`.
  - `fit_voce(gamma, tau_c) -> { τ_0, τ_s, γ_c }`.
  - `fit_basquin(cycles, stress) -> { σ'_f, b }`.
  - `fit_avrami(t, f) -> { k, n }`.
  - `fit_coffin_manson(cycles, strain) -> { ε'_f, c }`.
- **Generic Levenberg–Marquardt** for everything else — gradient and
  residual-vector API.

Driving example: `examples/backlog-modules-demo` recovers ground-truth
parameters to ~1 part in 10⁴.

## Example

```rust
use tpt_mat_inverse::fit_voce;

let params = fit_voce(&gamma, &tau_c);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).