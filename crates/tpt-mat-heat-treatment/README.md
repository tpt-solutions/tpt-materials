# tpt-mat-heat-treatment

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-heat-treatment.svg)](https://crates.io/crates/tpt-mat-heat-treatment)
[![Documentation](https://docs.rs/tpt-mat-heat-treatment/badge.svg)](https://docs.rs/tpt-mat-heat-treatment)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Heat-treatment schedule simulation: anneal / quench / temper / age /
> solution-treat, CALPHAD + TTT-driven final microstructure, Maynier
> hardness regression.

`tpt-mat-heat-treatment` composes `tpt-mat-phase-transform`,
`tpt-mat-calphad`, and `tpt-mat-grain-growth` into the user-facing heat-
treatment workflow:

- `HeatTreatmentProcess` enum — `Annealing`, `Quenching { medium }`,
  `Tempering`, `Aging`, `SolutionTreatment`.
- `simulate(&HeatTreatmentProcess) -> HeatTreatmentResult` — final
  phase fractions + grain size, driven by CALPHAD equilibrium and
  TTT / CCT kinetics.
- `predict_hardness(&PredictedMicrostructure) -> f64` — rule-of-mixtures
  / Maynier-type regression.

Driving case: `examples/heat-treatment-steel` runs a 4-step
anneal / quench / temper / age schedule on a low-alloy steel and
reports final phase fractions + Vickers hardness.

## Example

```rust
use tpt_mat_heat_treatment::{HeatTreatmentProcess, simulate};

let schedule = vec![
    HeatTreatmentProcess::Annealing { temperature: 1173.0, duration: 3600.0 },
    HeatTreatmentProcess::Quenching { medium: "water".into() },
    HeatTreatmentProcess::Tempering { temperature: 873.0, duration: 3600.0 },
    HeatTreatmentProcess::Aging { temperature: 423.0, duration: 14_400.0 },
];
let result = simulate(&schedule, &MaterialId::default());
let hv = predict_hardness(&result);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).