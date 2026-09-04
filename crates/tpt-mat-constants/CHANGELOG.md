# Changelog — tpt-mat-constants

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- CODATA physical constants (`R_GAS`, `K_BOLTZMANN`, `N_AVOGADRO`,
  `STEFAN_BOLTZMANN`, `PLANCK`).
- Atomic masses and ground-state properties for the first 92 elements.
- Bravais-lattice and Burgers-vector reference data for FCC / BCC / HCP.
- Unit-conversion helpers (eV ↔ J, GPa ↔ Pa, Å ↔ m).

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-constants-v0.1.0