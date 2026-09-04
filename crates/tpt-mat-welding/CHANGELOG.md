# Changelog — tpt-mat-welding

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `WeldModel { base_metal, filler_metal, process }`.
- `heat_affected_zone(heat_input) -> HazResult` — HAZ width +
  peak-temperature profile.
- `predict_haz_microstructure(cooling_rate) -> HazMicrostructure` —
  grain coarsening + Avrami / KM phase fractions across 8 HAZ
  sub-zones.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-welding-v0.1.0