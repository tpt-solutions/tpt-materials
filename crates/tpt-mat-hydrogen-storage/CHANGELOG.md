# Changelog — tpt-mat-hydrogen-storage

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `HydrideType` enum — `MetalHydride`, `ChemicalHydride`,
  `PorousMaterial`.
- `HydrogenStorageMaterial { hydride_type, storage_capacity_wt_pct, absorption_kinetics }`.
- `pct_isotherm(temperature)` returning the PCT curve with plateau and
  van 't Hoff temperature dependence.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-hydrogen-storage-v0.1.0