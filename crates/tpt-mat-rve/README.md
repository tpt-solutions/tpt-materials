# tpt-mat-rve

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-rve.svg)](https://crates.io/crates/tpt-mat-rve)
[![Documentation](https://docs.rs/tpt-mat-rve/badge.svg)](https://docs.rs/tpt-mat-rve)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Representative Volume Element (RVE) framework for polycrystal
> homogenization; Bishop–Hill Taylor-factor solver.

`tpt-mat-rve` packages the analytical polycrystal homogenisation drivers on
top of `tpt-mat-homogenization`:

- `Rve`, `RveGrain` — per-grain orientation, volume fraction, stiffness.
- `HomogenizationScheme` enum — `Voigt`, `Reuss`, `SelfConsistent` driver.
- `Rve::rotated_stiffness` — 4th-order stiffness rotation into the
  sample frame.
- `RveStats` — n_grains, total volume fraction, unique-orientation
  count.
- `SimpleHomogenizer` — Voigt/Reuss convenience wrapper.
- **Bishop–Hill (1951) Taylor-factor solver** with L2 pseudo-inverse
  (`bishop_hill_taylor_factor_axis`, `bishop_hill_taylor_factor`).
  Documented limitation: recovers `M ≈ 2.0–2.5` for random FCC (vs
  Taylor's classical `3.06`); L1 Lemke solver is queued for Phase 8.

The Hill–Mandel macro-homogeneity condition is verified for both
Voigt (`uniform strain`) and Reuss (`uniform stress`) limit cases.

## Example

```rust
use tpt_mat_rve::{Rve, HomogenizationScheme, bishop_hill_taylor_factor};

let c_eff = Rve::new(grains).homogenize(HomogenizationScheme::Voigt);
let m = bishop_hill_taylor_factor(&[001, 0, 0], CrystalStructure::Fcc);  // M ≈ 2.67
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).