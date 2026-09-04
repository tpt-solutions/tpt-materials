# Changelog — tpt-mat-rve

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `Rve`, `RveGrain` data model with per-grain orientation, volume
  fraction, stiffness.
- `HomogenizationScheme::{Voigt, Reuss, SelfConsistent}` driver.
- `Rve::rotated_stiffness` 4th-order stiffness rotation into the sample
  frame.
- `RveStats` (n_grains, total volume fraction, unique-orientation
  count).
- `SimpleHomogenizer` Voigt/Reuss convenience wrapper.
- **Bishop–Hill (1951) Taylor-factor solver** with L2 pseudo-inverse.

## [0.1.1] — 2026-09-04

### Added
- `hill_mandel_voigt_uniform_strain_energy_consistency` and
  `hill_mandel_reuss_uniform_stress_energy_consistency` unit tests
  verifying the Hill–Mandel macro-homogeneity condition in both
  limit cases.

### Deferred
- Full LCP (Lemke) Bishop–Hill solver for `M = 3.06` recovery.
- FFT (Moulinec–Suquet 1998) homogenisation.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-rve-v0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-rve-v0.1.0