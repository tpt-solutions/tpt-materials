# Changelog — tpt-mat-grain-growth

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `GrainBoundaryMobility` — Read–Shockley low-angle mobility + HAGB
  multiplier.
- `GrainGrowthSolver` — curvature-driven front-tracking driver.
- `grain_size_distribution` returning `GrainSizeStats`.

## [0.1.1] — 2026-09-04

### Added
- `interfaces` module: GB energy (Read–Shockley), GBCD / CSL fractions,
  triple-junction Herring balance, Langmuir–McLean segregation.
- `recrystallization` module: static JMAK RX, Zener–Hollomon +
  Sellars–Tegart, dynamic RX (DrxKinetics, Cahn–Hagel form).
- `stereology` module: ASTM E112 linear intercept, area / volume
  fraction + counting uncertainty, Saltykov 3-D reconstruction.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-grain-growth-v0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-grain-growth-v0.1.0