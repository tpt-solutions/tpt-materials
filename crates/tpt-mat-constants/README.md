# tpt-mat-constants

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-constants.svg)](https://crates.io/crates/tpt-mat-constants)
[![Documentation](https://docs.rs/tpt-mat-constants/badge.svg)](https://docs.rs/tpt-mat-constants)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Physical constants and periodic-table data.

`tpt-mat-constants` centralises every numerical constant used across the
`tpt-mat-*` domain crates:

- CODATA physical constants (`R_GAS`, `K_BOLTZMANN`, `N_AVOGADRO`,
  `STEFAN_BOLTZMANN`, `PLANCK`, …).
- Atomic masses and ground-state properties for the first 92 elements
  (sufficient for every alloy and ceramic system the workspace targets).
- Crystal-structure reference data (Bravais lattice name, point-group
  Schoenflies / Hermann–Mauguin symbols, prototype material).
- Burgers-vector magnitudes for FCC / BCC / HCP metals and a few
  intermetallics.
- A handful of unit conversions (eV ↔ J, GPa ↔ Pa, Å ↔ m).

All constants are `pub const` items so they can be used in `const fn`
contexts and inlined into hot loops.

## Example

```rust
use tpt_mat_constants::{R_GAS, BURGERS_VECTOR_FCC};

let d = 1.0e-3 * (-150_000.0 / (R_GAS * 1173.0)).exp();
let burgers = BURGERS_VECTOR_FCC; // m, scaled by `a / sqrt(2)`
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).