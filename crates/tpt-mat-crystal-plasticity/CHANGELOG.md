# Changelog — tpt-mat-crystal-plasticity

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `CrystalPlasticityModel` aggregating crystal structure, slip systems,
  hardening law, rate sensitivity, and elastic tensor.
- `RateSensitivity`, `power_law_slip_rate`.
- `resolved_shear_stresses` and `viscoplastic_velocity_gradient`
  (L^p = Σ_α γ̇^α s^α ⊗ n^α).
- `solve_increment_single_point` — radial-return kernel with hardening
  update.
- `CpFemSolver`, `BoundaryConditions`, `LoadStep`, `CpFemResult`,
  `PlasticIncrement` glue for single-point integration.

### Deferred
- Full FEM assembly + Newton–Raphson is wired through optional
  `tpt-fem-solve` / `tpt-fem-assembly` dependencies (off by default);
  no integration test yet because the `tpt-fem` mesh handle types are
  not in this workspace.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-crystal-plasticity-v0.1.0