# tpt-mat-polymer

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-polymer.svg)](https://crates.io/crates/tpt-mat-polymer)
[![Documentation](https://docs.rs/tpt-mat-polymer/badge.svg)](https://docs.rs/tpt-mat-polymer)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Polymer-chain constitutive equations — Arruda–Boyce 8-chain rubber
> elasticity (inverse Langevin), Worm-Like-Chain Marko–Siggia
> force–extension, Freely-Jointed-Chain inverse-Langevin.

`tpt-mat-polymer` provides the standard single-chain constitutive
models used for elastomers and biopolymers:

- `ChainModel` enum — `FreelyJointedChain { num_segments, segment_length }`,
  `WormLikeChain { persistence_length, contour_length }`,
  `ArrudaBoyce { n_segments, shear_modulus }`.
- `PolymerModel { chain_model, crosslink_density }`.
- `stress_strain(stretch) -> f64` — Arruda–Boyce 8-chain via inverse
  Langevin; WLC force–extension.

Verified: Arruda–Boyce → neo-Hookean at small stretch (Taylor expansion
of the inverse Langevin gives the standard μ = nₖᵦT shear modulus).

## Example

```rust
use tpt_mat_polymer::{PolymerModel, ChainModel, stress_strain};

let rubber = PolymerModel::new(ChainModel::ArrudaBoyce { n_segments: 8.0, shear_modulus: 0.3e6 }, 1.0);
let sigma = stress_strain(&rubber, 1.5);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).