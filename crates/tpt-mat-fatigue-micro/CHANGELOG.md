# Changelog — tpt-mat-fatigue-micro

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `MicrostructuralFatigue { rve, criterion }`.
- `FatigueCriterion` enum — `Findley`, `FatemiSocie`,
  `SmithWatsonTopper`, `CrystallographicSlip`.
- `fatigue_indicator_parameter` — per-grain FIP field from a CP-FEM
  result.
- `predict_crack_initiation` — cycles-to-initiation, critical grain,
  critical_location.
- Cycle-by-cycle CP-FEM driver accumulating plastic slip at grain
  boundaries.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-fatigue-micro-v0.1.0