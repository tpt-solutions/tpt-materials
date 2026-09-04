# Changelog — tpt-mat-solidification

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `AnisotropyModel` enum — `Cubic4Fold`, `Hexagonal6Fold`, `Isotropic`.
- `SolidificationSolver` — Kobayashi wrapper with anisotropy + undercooling.
- `simulate_dendrite` + `snapshot` returning solid fraction, tip velocity,
  primary / secondary arm spacing.
- `secondary_arm_spacing` post-processor.

### Verified
- Solid fraction grows monotonically from a single seed under undercooling.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-solidification-v0.1.0