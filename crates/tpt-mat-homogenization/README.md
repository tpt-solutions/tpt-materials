# tpt-mat-homogenization

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-homogenization.svg)](https://crates.io/crates/tpt-mat-homogenization)
[![Documentation](https://docs.rs/tpt-mat-homogenization/badge.svg)](https://docs.rs/tpt-mat-homogenization)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Analytical micromechanics: Voigt / Reuss / VRH bounds, Hashin–Shtrikman
> variational bounds, Eshelby spherical inclusion, dilute strain
> concentration.

`tpt-mat-homogenization` is the analytical micromechanics core the rest of
the homogenisation / RVE / composite stack consumes:

- `voigt`, `reuss`, `voigt_reuss_average` — N-phase first-order bounds on
  stiffness.
- `voigt_reuss_bounds` — Hill `M_VRH` for isotropic `K, G`.
- `hashin_shtrikman_two_phase`, `hashin_shtrikman_k_g` — variational
  bounds; enclose Voigt for low-contrast mixtures.
- `EshelbySpherical`, `eshelby_spherical(matrix_E, matrix_ν)` —
  closed-form spherical-inclusion Eshelby tensor.
- `dilute_strain_concentration` — `A = [I + S C_0⁻¹ ΔC]⁻¹`, 6×6
  strain-concentration tensor.
- `k_from_e_nu`, `g_from_e_nu` — isotropic `K, G` from `E, ν`.

The tensor machinery is generalised in `tpt-mat-thermal` to any
symmetric-positive transport property (conductivity, diffusivity,
permittivity).

## Example

```rust
use tpt_mat_homogenization::{voigt_reuss_bounds, hashin_shtrikman_two_phase};

let (k_vrh_lo, k_vrh_hi) = voigt_reuss_bounds(70.0, 0.33, 200.0, 0.30, 0.5);
let (k_hs_lo, k_hs_hi) = hashin_shtrikman_two_phase(70.0, 0.33, 200.0, 0.30, 0.5);
assert!(k_hs_lo <= k_vrh_lo && k_vrh_hi <= k_hs_hi);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).