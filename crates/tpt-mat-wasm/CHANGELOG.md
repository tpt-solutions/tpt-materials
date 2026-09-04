# Changelog — tpt-mat-wasm

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- Crate scaffold with `wasm-bindgen` dependency wiring for
  `tpt-mat-core`, `tpt-mat-constants`, `tpt-mat-crystallography`,
  `tpt-math-linalg-fixed`.
- `WasmRveSolver` — Voigt / Reuss / self-consistent RVE homogenisation.
- `WasmPhaseField` — `step` and `get_order_parameter` accessors.

### Deferred
- Full WASM coverage of crystal plasticity, microstructural fatigue,
  hydrogen embrittlement (Phase 8 follow-up).

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-wasm-v0.1.0