# tpt-mat-calphad

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-calphad.svg)](https://crates.io/crates/tpt-mat-calphad)
[![Documentation](https://docs.rs/tpt-mat-calphad/badge.svg)](https://docs.rs/tpt-mat-calphad)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> CALPHAD-style thermodynamic models: Redlich–Kister polynomials,
> end-member + ideal mixing + excess Gibbs energy, Muggianu sublattice
> model, two-phase common-tangent equilibrium.

`tpt-mat-calphad` provides the analytical thermodynamic descriptors used
by every phase-transformation, precipitation, and homogenisation crate:

- `RedlichKister { coefficients }` — `Σ L_ν (x_A − x_B)^ν` polynomial
  with analytic derivative.
- `GibbsEnergyModel { end_members, ideal, excess }` — full Gibbs-energy
  expression and `derivative(x)`.
- `SublatticeModel` — Muggianu multi-sublattice model with ideal
  configurational Gibbs energy.
- `PhaseDiagram` + `two_phase_equilibrium` — grid-search common-tangent
  construction, returning `PhaseBoundary` per temperature.

The end-members + interaction parameters come from a user-supplied
database; `tpt-mat-calphad` does not bundle a TDB file. The
`examples/calphad-binary-phase-diagram` example feeds a hand-coded
binary A–B Redlich–Kister parameter set and prints the equilibrium
  phase-fraction curve.

## Example

```rust
use tpt_mat_calphad::{RedlichKister, GibbsEnergyModel, two_phase_equilibrium};

let l0 = RedlichKister::new(vec![(-150_000.0, 0.0)]); // regular solution
let g_liquid = GibbsEnergyModel::ideal_2(8_000.0, 12_000.0);
let g_bcc = GibbsEnergyModel::excess(g_liquid, l0);
let boundary = two_phase_equilibrium(&g_liquid, &g_bcc, 800.0);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).