# Changelog — tpt-mat-thermal

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `effective_conductivity` — series / parallel / VRH / Hashin–Shtrikman
  / self-consistent / Maxwell–Garnett.
- `effective_cte` — Turner / Kerner / Rosen–Hashin bounds.
- `effective_specific_heat` — mass-weighted rule of mixtures.
- `effective_diffusivity` — Bruggeman / tortuosity for porous media.
- `interface_thermal_resistance` — Kapitza correction.

### Verified
- Conductivity HS bounds enclose the self-consistent estimate.
- Bounds collapse at zero contrast.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-thermal-v0.1.0