# Changelog — tpt-mat-calphad

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `RedlichKister` — `Σ L_ν (x_A − x_B)^ν` polynomial with analytic
  derivative.
- `GibbsEnergyModel` — end-member + ideal mixing + Redlich–Kister
  excess; `total(x)` and `derivative(x)`.
- `SublatticeModel` — Muggianu multi-sublattice model with ideal
  configurational Gibbs energy.
- `PhaseDiagram` + `two_phase_equilibrium` — grid-search common-
  tangent construction, returning `PhaseBoundary` per temperature.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-calphad-v0.1.0