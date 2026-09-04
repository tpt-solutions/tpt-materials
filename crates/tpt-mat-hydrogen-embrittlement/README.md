# tpt-mat-hydrogen-embrittlement

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-hydrogen-embrittlement.svg)](https://crates.io/crates/tpt-mat-hydrogen-embrittlement)
[![Documentation](https://docs.rs/tpt-mat-hydrogen-embrittlement/badge.svg)](https://docs.rs/tpt-mat-hydrogen-embrittlement)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Hydrogen transport with trapping + embrittlement criteria — Oriani
> local-equilibrium `D_eff`, McNabb–Foster kinetic trapping,
> stress-driven uphill hydrogen flux, HEDE / HELP thresholds,
> susceptibility index.

`tpt-mat-hydrogen-embrittlement` bridges the diffusion and fracture
crates to give a complete pipeline from hydrogen ingress to
embrittlement indicator:

- `HydrogenTransport` — Fick + trapping (Oriani local equilibrium,
  McNabb–Foster kinetic trapping with `N_T`, `E_B`, occupancy `θ_T`).
- `stress_driven_diffusion` —
  `∂C/∂t = ∇·(D∇C − D C V_H ∇σ_h / RT)` (hydrostatic-stress uphill
  flux) on a `tpt-science` grid.
- `EmbrittlementCriterion` — HEDE critical-lattice-decohesion and
  HELP local-plasticity indicators; threshold `C_crit(σ_h)`.
- `susceptibility_index` — local H concentration + triaxiality field
  (feeds `tpt-mat-fracture` / `tpt-mat-fatigue-micro`).

Verified: trapping retards effective diffusivity by
`D_eff = D_L / (1 + ∂C_T/∂C_L)`.

## Example

```rust
use tpt_mat_hydrogen_embrittlement::{HydrogenTransport, OrianiTrapping, EmbrittlementCriterion};

let h = HydrogenTransport::default();
let c_eff = h.oriani_effective_diffusivity(&state);
let ind = h.susceptibility_index(&stress_field, &c_eff);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).