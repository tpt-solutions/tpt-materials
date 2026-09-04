# Changelog — tpt-mat-corrosion

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `ElectrodeKinetics`, `CorrosionModel`, `CorrosionRate`.
- `corrosion_rate(T, pH) -> CorrosionRate` — Butler–Volmer
  mixed-potential solve.
- `polarization_curve(potential_range) -> PolarizationCurve` — Tafel
  anodic / cathodic branches.

## [0.1.1] — 2026-09-04

### Added
- `oxidation` submodule: Wagner parabolic scale growth
  (`ParabolicRateConstant`, `ScaleGrowth`, `DopingEffect`,
  `BreakawayCriterion`, `LinearBreakawayRate`, `OxidationModel`).
- Crate re-themed as "environmental degradation".

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-corrosion-v0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-corrosion-v0.1.0