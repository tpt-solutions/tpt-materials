# Changelog — tpt-mat-fracture

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `StressIntensityFactor` with `K_I`, `K_II`, `K_III` and centre /
  edge / Mode II / Mode III geometry factors.
- `EnergyReleaseRate` — `G`, Irwin `K → G`, J-integral domain form.
- `CohesiveZoneModel` — bilinear / exponential traction–separation;
  Benzeggagh–Kenane mixed mode.
- `PhaseFieldFracture` — Griffith / AT1 / AT2 regularised functionals.
- `fracture_toughness_transition` — ASTM E1921 master curve.

### Verified
- Phase-field fracture recovers the Griffith load for a 1-D bar.
- Irwin `K → G` consistency.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-fracture-v0.1.0