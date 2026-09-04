# Changelog — tpt-mat-machine-learning

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `MlModel` enum — `PropertyPredictor { features }`,
  `PhasePredictor`, `SurrogateModel`.
- OLS / ridge / kernel-ridge regression with closed-form normal
  equations.
- Polynomial-feature augmentation (degree ≤ 5).
- Multi-class softmax phase predictor.

### Verified
- Surrogate reproduces a CP-FEM / homogenisation sweep within
  tolerance.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-machine-learning-v0.1.0