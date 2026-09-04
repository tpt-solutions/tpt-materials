# Changelog — tpt-mat-precipitation

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `ClassicalNucleation` — Turnbull–Fisher `I = I_0 exp(−ΔG*/kT) exp(−Q/kT)`.
- `KwnModel` — Kampmann–Wagner numerical size-class driver.
- `LswCoarsing` — `R̄³ − R̄_0³ = K t`.
- `PrecipitateState` — number density, mean radius, volume fraction,
  matrix supersaturation.
- `strengthening_increment` — Orowan bypass + shearing.

### Verified
- KWN conserves solute mass.
- Late-stage slope → LSW `t^{1/3}`.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-precipitation-v0.1.0