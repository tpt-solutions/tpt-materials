# Changelog — tpt-mat-diffusion

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `ArrheniusDiffusivity { d_0, q }` — `D = D_0 exp(−Q / (R T))`.
- `DiffusionSolver` — single-component forward-Euler Fick's second law
  on regular 2D grid with Neumann zero-flux boundaries; mass
  conserved to ~1 × 10⁻⁶.
- `MultiComponentDiffusionSolver` — independent species, per-species
  CFL stability limit.
- `GrainBoundaryDiffusion` — Fisher Regime-A effective diffusivity.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-diffusion-v0.1.0