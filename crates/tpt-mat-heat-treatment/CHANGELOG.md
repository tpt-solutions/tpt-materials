# Changelog — tpt-mat-heat-treatment

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `HeatTreatmentProcess` enum — `Annealing`, `Quenching { medium }`,
  `Tempering`, `Aging`, `SolutionTreatment`.
- `simulate(&[HeatTreatmentProcess]) -> HeatTreatmentResult` — final
  phase fractions + grain size, CALPHAD + TTT driven.
- `predict_hardness(&PredictedMicrostructure) -> f64` — rule-of-mixtures
  / Maynier-type regression.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-heat-treatment-v0.1.0