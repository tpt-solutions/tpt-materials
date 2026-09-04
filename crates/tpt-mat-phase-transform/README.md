# tpt-mat-phase-transform

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-phase-transform.svg)](https://crates.io/crates/tpt-mat-phase-transform)
[![Documentation](https://docs.rs/tpt-mat-phase-transform/badge.svg)](https://docs.rs/tpt-mat-phase-transform)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Solid-state phase transformation kinetics: Avrami (JMAK),
> Koistinen–Marburger (martensite), TTT / CCT diagrams, additive
> thermal-history integration.

`tpt-mat-phase-transform` is the steel / alloy heat‑treatment kinetics
crate:

- `AvramiModel { k, n }` — `f = 1 − exp(−k t^n)` with Arrhenius
  temperature-dependent `k(T)`.
- `KoistinenMarburger { alpha, m_s }` — diffusionless martensite
  `f = 1 − exp(−α (M_s − T)⁺)` with retained-at-M_s cap.
- `TransformationSolver` — Scheil-additive integration over piecewise-
  linear thermal histories (combined Avrami + KM).
- `time_to_fraction(f, T) -> f64` — TTT / CCT diagram helper.

Driving example: `examples/phase-transform-jmak` shows Avrami fraction
reaching 1.0 in 600 s at 900 K and a martensite fraction of 0.983 after
2000 s of cooling from 1100 K → 200 K.

## Example

```rust
use tpt_mat_phase_transform::{AvramiModel, KoistinenMarburger, TransformationSolver};

let avrami = AvramiModel { k: 1e-3, n: 2.0 };
let km = KoistinenMarburger { alpha: 0.011, m_s: 600.0 };
let solver = TransformationSolver::new(avrami, km);
let f = solver.transform_fraction(&thermal_history, 1200.0);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).