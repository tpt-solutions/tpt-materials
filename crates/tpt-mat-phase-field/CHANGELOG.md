# Changelog — tpt-mat-phase-field

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `PhaseFieldSolver` explicit forward-Euler driver.
- `PhaseFieldModel` enum — `AllenCahn`, `CahnHilliard`, `Kobayashi`,
  `MultiPhase`.
- `BulkEnergy` — `DoubleWell`, `Polynomial`, `RegularSolution`.
- `FreeEnergyFunctional` 2D central-difference evaluator.
- `step_allen_cahn`, `step_cahn_hilliard`, `step_kobayashi`,
  `step_multi_phase`.
- `PhaseFieldResult` — order parameter, concentration, temperature,
  free energy, interface area.

### Verified
- Free energy monotonically decreases under Allen-Cahn step.
- Mean concentration conserved under Cahn–Hilliard step to ~1 × 10⁻⁶.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-phase-field-v0.1.0