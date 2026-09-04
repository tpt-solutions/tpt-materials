# Changelog — tpt-mat-phase-transform

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `AvramiModel` — JMAK `f = 1 − exp(−k t^n)` with Arrhenius `k(T)`.
- `KoistinenMarburger` — diffusionless `f = 1 − exp(−α (M_s − T)⁺)`
  with retained-at-M_s cap.
- `TransformationSolver` — Scheil-additive integration over
  piecewise-linear thermal histories; combined Avrami + KM.
- `time_to_fraction` TTT / CCT diagram helper.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-phase-transform-v0.1.0