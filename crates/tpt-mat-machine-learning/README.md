# tpt-mat-machine-learning

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-machine-learning.svg)](https://crates.io/crates/tpt-mat-machine-learning)
[![Documentation](https://docs.rs/tpt-mat-machine-learning/badge.svg)](https://docs.rs/tpt-mat-machine-learning)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Lightweight ML surrogates for materials-property prediction:
> OLS / ridge / kernel-ridge regression, polynomial-feature
> augmentation, multi-class softmax phase predictor.

`tpt-mat-machine-learning` is the on-workspace replacement for the
`linfa` ecosystem when the only thing required is small, transparent,
scikit-learn-style surrogates:

- `MlModel` enum — `PropertyPredictor { features }`,
  `PhasePredictor`, `SurrogateModel`.
- `train(&[(Vec<f64>, f64)]) -> TrainingResult` — closed-form normal
  equations with optional L2 ridge penalty.
- `predict(&Composition) -> f64`.
- Polynomial-feature augmentation (degree ≤ 5).
- Multi-class softmax phase predictor.
- RBF kernel ridge regression.

Verified: surrogate reproduces a CP-FEM / homogenisation sweep
within tolerance (see `examples/ml-surrogate`).

## Example

```rust
use tpt_mat_machine_learning::{MlModel, train, predict};

let model = MlModel::PropertyPredictor { features: 2 };
let trained = train(&model, &dataset);
let y = trained.predict(&[0.95, 0.05]);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).