# tpt-mat-texture

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-texture.svg)](https://crates.io/crates/tpt-mat-texture)
[![Documentation](https://docs.rs/tpt-mat-texture/badge.svg)](https://docs.rs/tpt-mat-texture)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Crystallographic texture analysis: EBSD import, pole figures,
> orientation distribution function, Taylor factor.

`tpt-mat-texture` is the analysis companion to `tpt-mat-crystal-plasticity`.
It consumes orientation data and produces the standard texture
quantities:

- `TextureAnalyzer` — orientation + weight container with
  `from_ebsd(orientations, weights)` and `from_random(n)`.
- `PoleFigure` — equal-area / stereographic projection with
  pre-binned crystal-direction sets (`Fcc111`, `Fcc200`, `Bcc110`,
  `Hcp0001`, `Hcp10T10`) and a custom-direction builder.
- `OrientationDistributionFunction` — geodesic-Gaussian KDE on SO(3).
- `taylor_factor` — single-Schmid proxy used for screening; the
  Bishop–Hill L2 solver lives in `tpt-mat-rve`.

The Taylor-factor proxy over-counts because it assumes single-system
activation — the true Taylor factor requires finding the 5 slip
systems that accommodate the macroscopic strain. The full Bishop–Hill
solver (L1 Lemke) is queued for the Phase 5 / Phase 8 workstream.

## Example

```rust
use tpt_mat_texture::{TextureAnalyzer, PoleFigure, CrystalDirection};

let ta = TextureAnalyzer::from_random(1024);
let pf = PoleFigure::equal_area(&ta, CrystalDirection::Fcc111, 256);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).