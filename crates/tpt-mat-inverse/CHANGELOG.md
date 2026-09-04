# Changelog — tpt-mat-inverse

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-04

### Added
- Closed-form fits: Arrhenius, Norton–Bailey, Voce, Basquin S–N,
  Avrami, Coffin–Manson.
- Generic `LevenbergMarquardt` driver for arbitrary residual-vector
  problems.
- `examples/backlog-modules-demo` — calibration fits recover
  ground-truth parameters to 1 part in 10⁴.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-inverse-v0.1.0