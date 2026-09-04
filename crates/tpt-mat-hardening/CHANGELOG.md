# Changelog — tpt-mat-hardening

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `Voce` — `τ_c = τ_0 + (τ_s − τ_0)(1 − e^{−γ/γ_c}) + θ_0 γ`.
- `PowerLaw` — `Δτ_c = h_0 (τ_s − τ_c) |Δγ| / τ_s`.
- `KocksMecking` — implicit exponential saturation.
- `CombinedHardening` — additive combination of two laws.
- `LatentHardeningMatrix` — `h_{αβ}` latent-hardening matrix with
  `q = 1.0` coplanar / `q = 1.4` non-coplanar defaults.

## [0.1.1] — 2026-09-04

### Added
- `DislocationDensityHardening` wired into `Hardening::DislocationDensity`
  variant, with `HardeningState.extra` carrying the per-system density
  triple `(ρ_SSD, ρ_GND, ρ_forest)`.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-hardening-v0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-hardening-v0.1.0