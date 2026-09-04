# tpt-mat-fatigue

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-fatigue.svg)](https://crates.io/crates/tpt-mat-fatigue)
[![Documentation](https://docs.rs/tpt-mat-fatigue/badge.svg)](https://docs.rs/tpt-mat-fatigue)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Classical analytical fatigue: Basquin S–N, Coffin–Manson LCF,
> Manson–Coffin–Basquin strain-life, Walker mean-stress correction,
> Paris / Forman / Walker / NASGRO crack growth, ASTM E1049 rainflow
> counting.

`tpt-mat-fatigue` is the analytical high-cycle / low-cycle fatigue stack
that complements the microstructural FIP-based model in
`tpt-mat-fatigue-micro`:

- `Basquin { sigma_f_prime, b }` — high-cycle S–N `σ_a = σ'_f (2N_f)^b`.
- `CoffinManson { epsilon_f_prime, c }` — LCF `Δε_p / 2 = ε'_f (2N_f)^c`.
- `StrainLife` — Manson–Coffin–Basquin total strain-life with elastic
  and plastic components.
- `Walker` — mean-stress correction `σ_ar = σ_a (1 − R)^(1 − m_walker)`.
- Crack growth: `Paris`, `Forman`, `Walker`, `NASGRO`.
- `rainflow` — ASTM E1049 rainflow cycle-counting on a load history.

## Example

```rust
use tpt_mat_fatigue::{Basquin, paris_da_dn, rainflow};

let n_f = basquin_cycles_to_failure(300.0, 700.0, -0.1);  // 1e6+ cycles
let cycles = rainflow(&load_history);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).