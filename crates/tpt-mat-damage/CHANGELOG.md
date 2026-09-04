# Changelog — tpt-mat-damage

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- Kachanov effective-stress concept.
- Lemaitre ductile damage evolution.
- Kachanov creep-damage coupling.
- Miner linear damage accumulation.
- Chaboche placeholder.

## [0.1.1] — 2026-09-03

### Added
- `GursonTvergaardNeedleman` porous plasticity (`f_0`, `f_c`, `f_f`,
  `q_1`, `q_2`, `q_3`, `f_n`, `s_n`, `e_n`).
- `yield_function(stress, porosity, matrix_yield)`.
- `update_porosity` — void growth + Chu–Needleman nucleation.

### Verified
- GTN yield surface → von Mises as `f → 0`.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-damage-v0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-damage-v0.1.0