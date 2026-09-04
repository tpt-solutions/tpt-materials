# Changelog — tpt-mat-core

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `Composition` (mole / mass fractions, normalisation).
- `Phase`, `Grain`, `Microstructure` data model.
- `Orientation` (passive rotation matrix in `SO(3)`).
- `Lattice` Bravais lattice descriptor (`Cubic`, `Hexagonal`, `Tetragonal`, …).

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-core-v0.1.0