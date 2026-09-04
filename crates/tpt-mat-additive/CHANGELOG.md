# Changelog — tpt-mat-additive

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `AmProcess` enum — `LaserPowderBedFusion`, `DirectedEnergyDeposition`,
  `ElectronBeamMelting` with linear + volumetric energy densities.
- `thermal_history(location) -> ThermalHistory` — Rosenthal
  moving-point-source or Gaussian-beam.
- `predict_microstructure` — Hunt CET map, phase fractions, texture,
  porosity.
- `residual_stress` — thermal-contraction eigenstrain.

### Verified
- Higher cooling rate → finer predicted grain size (Hall–Petch trend).

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-additive-v0.1.0