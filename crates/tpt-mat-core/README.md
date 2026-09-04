# tpt-mat-core

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-core.svg)](https://crates.io/crates/tpt-mat-core)
[![Documentation](https://docs.rs/tpt-mat-core/badge.svg)](https://docs.rs/tpt-mat-core)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Core material data model: microstructure, phase, grain, composition,
> orientation.

`tpt-mat-core` defines the foundational data structures that every other
`tpt-mat-*` crate builds on:

- `Composition` — mole / mass fractions with normalisation and helpers.
- `Phase` — phase identifier with `name`, `formula`, density, color tag.
- `Grain` — single crystallite with orientation, phase id, and shape data.
- `Microstructure` — collection of grains + phases + global metadata.
- `Orientation` — passive rotation matrix in `SO(3)` (column-major
  convention).
- `Lattice` — Bravais lattice type (`Cubic`, `Hexagonal`, …) with
  `a`, `c` lattice parameters and angle conventions used by
  `tpt-mat-crystallography`.

This crate has no PDE or constitutive semantics. It is the canonical
"shared vocabulary" crate — every other `tpt-mat-*` domain crate accepts
and returns types defined in `tpt-mat-core`.

## Example

```rust
use tpt_mat_core::{Composition, Phase, Grain, Microstructure, Orientation};

let comp = Composition::from_mole_fractions(&[("Fe", 0.95), ("C", 0.05)])?;
let phase = Phase::new("Ferrite", "Fe-C", 7.87);
let grain = Grain::new(0, Orientation::identity(), 1.0e-6);
let ms = Microstructure::new(vec![phase], vec![grain], comp);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).