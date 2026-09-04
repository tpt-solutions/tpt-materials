# Changelog — tpt-mat-hydrogen-embrittlement

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `HydrogenTransport` — Fick + trapping (Oriani local equilibrium,
  McNabb–Foster kinetic trapping with `N_T`, `E_B`, `θ_T`).
- `stress_driven_diffusion` — `∂C/∂t = ∇·(D∇C − D C V_H ∇σ_h / RT)`.
- `EmbrittlementCriterion` — HEDE / HELP thresholds `C_crit(σ_h)`.
- `susceptibility_index` from local H + triaxiality.

### Verified
- Trapping retards `D_eff = D_L / (1 + ∂C_T/∂C_L)`.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-hydrogen-embrittlement-v0.1.0