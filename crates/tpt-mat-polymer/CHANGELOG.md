# Changelog — tpt-mat-polymer

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `ChainModel` enum — `FreelyJointedChain`, `WormLikeChain`,
  `ArrudaBoyce`.
- `PolymerModel { chain_model, crosslink_density }`.
- `stress_strain(stretch) -> f64` — Arruda–Boyce 8-chain via inverse
  Langevin; WLC Marko–Siggia.

### Verified
- Arruda–Boyce → neo-Hookean at small stretch.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-polymer-v0.1.0