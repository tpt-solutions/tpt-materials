# Changelog — tpt-mat-homogenization

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `voigt`, `reuss`, `voigt_reuss_average` — N-phase first-order bounds.
- `voigt_reuss_bounds` — Hill `M_VRH` for isotropic `K, G`.
- `hashin_shtrikman_two_phase` + `hashin_shtrikman_k_g` — variational
  bounds.
- `EshelbySpherical` + `eshelby_spherical(matrix_E, matrix_ν)` —
  closed-form spherical-inclusion tensor.
- `dilute_strain_concentration` — 6×6 strain-concentration tensor.
- `k_from_e_nu`, `g_from_e_nu` isotropic helpers.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-homogenization-v0.1.0