# Changelog — tpt-mat-crystallography

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `CrystalStructure` enum (`Fcc`, `Bcc`, `Hcp`) with slip / twinning
  system enumerations.
- `slip_systems(crystal)` returning 12 FCC, 24 BCC, 12 HCP systems.
- `schmid_tensor(slip) -> Mat3` symmetric rank-1 `s ⊗ n`.
- `resolved_shear_stresses(stress, slip_systems) -> Vec<f64>`.
- `SymmetricFourthOrder` 6×6 stiffness with cubic, hexagonal, isotropic
  constructors.

### Verified
- FCC slip system count = 12.
- Schmid-tensor symmetry to machine precision.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-crystallography-v0.1.0