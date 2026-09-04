# Changelog — tpt-science

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `Grid1D`, `Grid2D`, `Grid3D` regular-grid types with `dx`, shape, indexing.
- `compute_laplacian` on 1D / 2D / 3D grids with Neumann zero-flux,
  periodic, and Dirichlet boundary conditions.
- `compute_biharmonic` (`∇⁴`) for Cahn–Hilliard and phase-field fracture.
- `compute_gradient` central-difference vector field.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-science-v0.1.0