# Changelog — tpt-mat-dislocation

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `DislocationDensityState` — per-system SSD / GND / forest.
- `KocksMeckingEvolution` — `dρ/dγ = k_1 √ρ − k_2 ρ`.
- Taylor stress `τ = α μ b √ρ`.
- Armstrong–Frederick back-stress from GND gradients.
- `gnd_from_curvature` from Nye tensor.

### Verified
- Single-slip response reproduces Voce-like saturation.
- `ρ` stays non-negative under the discrete update.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-dislocation-v0.1.0