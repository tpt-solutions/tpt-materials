# Changelog — tpt-mat-hydrogel

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- Flory–Rehner swelling equilibrium (mixing + elastic osmotic balance)
  with root-find over polymer volume fraction.
- Poroelastic swelling kinetics — Fickian solvent uptake on a
  `tpt-science` 1-D grid.
- `equilibrium_swelling_ratio(chi, crosslink_density) -> f64`.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-hydrogel-v0.1.0