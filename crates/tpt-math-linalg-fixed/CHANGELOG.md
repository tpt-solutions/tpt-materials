# Changelog — tpt-math-linalg-fixed

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `Vec2`, `Vec3` fixed-size Euclidean vectors.
- `Mat3` column-major 3×3 matrix with transpose, determinant, inverse,
  `sym` / `skew` decomposition.
- `SymmetricFourthOrder` 6×6 Voigt-notation elastic stiffness with
  cubic / hexagonal / isotropic constructors.

### Fixed
- `Mat3::inverse`, `Mat3::sym`, and `Mat3::skew` had wrong index mappings
  (column-major storage was treated as row-major). Corrected in
  `src/mat3.rs`; existing `transpose_inverse` and `sym_skew_decompose`
  tests now pass.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-math-linalg-fixed-v0.1.0