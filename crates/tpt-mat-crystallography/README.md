# tpt-mat-crystallography

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-crystallography.svg)](https://crates.io/crates/tpt-mat-crystallography)
[![Documentation](https://docs.rs/tpt-mat-crystallography/badge.svg)](https://docs.rs/tpt-mat-crystallography)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Crystal structures, slip systems, Schmid tensors, resolved shear stress.

`tpt-mat-crystallography` is the home of the classical single-crystal
kinematics and elasticity primitives the entire crystal-plasticity stack
depends on:

- `CrystalStructure` — `Fcc`, `Bcc`, `Hcp` with the canonical slip /
  twinning system enumerations (12 FCC, 12 BCC {110}, 24 BCC {112},
  3 basal + 3 prismatic HCP, …).
- `slip_systems(crystal) -> Vec<SlipSystem>` — `(s, n)` direction and
  normal pairs.
- `schmid_tensor(slip) -> Mat3` — symmetric rank-1 tensor `s ⊗ n`,
  verified to be symmetric to 1 part in 10⁻¹².
- `resolved_shear_stresses(stress, slip_systems) -> Vec<f64>` — per-system
  τ = σ_ij s_i n_j.
- `SymmetricFourthOrder` — 6×6 stiffness with cubic, hexagonal, isotropic
  constructors used by `tpt-mat-crystal-plasticity` for the elasticity step.

Verification: 12 FCC slip systems, 24 BCC, 12 HCP, Schmid-tensor symmetry
checked under `cargo test`.

## Example

```rust
use tpt_mat_crystallography::{CrystalStructure, slip_systems, schmid_tensor};

let systems = slip_systems(CrystalStructure::Fcc);
assert_eq!(systems.len(), 12);
let schmid = schmid_tensor(systems[0]);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).