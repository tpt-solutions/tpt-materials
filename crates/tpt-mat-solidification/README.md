# tpt-mat-solidification

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-solidification.svg)](https://crates.io/crates/tpt-mat-solidification)
[![Documentation](https://docs.rs/tpt-mat-solidification/badge.svg)](https://docs.rs/tpt-mat-solidification)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Dendritic solidification: Kobayashi wrapper with anisotropy models
> (cubic / hexagonal / isotropic) + tip velocity and arm-spacing metrics.

`tpt-mat-solidification` packages the dendritic-growth use case around
`tpt-mat-phase-field`'s Kobayashi driver:

- `AnisotropyModel` enum — `Cubic4Fold`, `Hexagonal6Fold`, `Isotropic`.
- `SolidificationSolver` — wraps Kobayashi with the chosen anisotropy and
  the undercooling boundary condition.
- `simulate_dendrite(seed, num_steps) -> DendriteResult` — runs the
  solver and returns the solid fraction, tip velocity, and primary /
  secondary arm spacing.
- `secondary_arm_spacing` post-processor.

Verified: solid fraction grows monotonically from a single seed under
constant undercooling (see `examples/dendritic-solidification`).

## Example

```rust
use tpt_mat_solidification::{AnisotropyModel, SolidificationSolver};

let solver = SolidificationSolver::new(anisotropy::Cubic4Fold { strength: 0.05 }, 0.45);
let result = solver.simulate_dendrite(0, 500);
println!("v_tip = {}, λ₂ = {}", result.tip_velocity, result.secondary_arm_spacing);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).