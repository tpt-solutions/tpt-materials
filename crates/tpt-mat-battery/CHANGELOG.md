# Changelog — tpt-mat-battery

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `ActiveMaterial` with `BatteryChemistry` enum (NMC811 / NMC622 / LFP /
  NCA / graphite / silicon, …).
- `DegradationMechanism` enum: `SeiGrowth`, `ParticleCracking`,
  `LithiumPlating`, `TransitionMetalDissolution`.
- `simulate_diffusion_stress(c_rate, num_cycles)`.
- `capacity_fade_curve(cycles, temperature) -> DegradationCurve` —
  cycles, capacity_retention, resistance_growth.

### Verified
- Capacity retention monotonically decreasing.
- √t SEI-limited fade at low C-rate.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-battery-v0.1.0