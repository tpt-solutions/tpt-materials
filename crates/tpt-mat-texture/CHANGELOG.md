# Changelog — tpt-mat-texture

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `TextureAnalyzer` with `from_ebsd` and `from_random` constructors.
- `PoleFigure` equal-area / stereographic projection with pre-binned
  crystal directions (`Fcc111`, `Fcc200`, `Bcc110`, `Hcp0001`,
  `Hcp10T10`) and custom-direction builder.
- `OrientationDistributionFunction` — geodesic-Gaussian KDE on SO(3).
- `taylor_factor` single-Schmid proxy (documented limitation; the
  Bishop–Hill L1 solver is queued for the Phase 5/8 workstream).

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-texture-v0.1.0