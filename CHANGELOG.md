# Changelog

All notable changes to `tpt-materials` are documented in this file. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Workspace scaffolding (Phase 0 of [`todo.md`](./todo.md)).
- Phase 1 crate `tpt-math-linalg-fixed` providing `Vec3`, `Mat3`, `SymMat3`,
  `Vec6`, and helpers for fourth-order elasticity/Schmid tensor packing.
- Phase 1 crate `tpt-mat-core` providing `MaterialMicrostructure`, `Phase`,
  `Grain`, `Composition`, `CompositionBasis`, `CrystalOrientation`, and
  `OrientationRepresentation` (`EulerBunge`, `Quaternion`, `RotationMatrix`,
  `Rodrigues`, `AxisAngle`).
- Phase 1 crate `tpt-mat-constants` providing `PhysicalConstants` (CODATA
  2018 values) and a minimal `PeriodicTable` with `atomic_mass` and
  `atomic_radius` lookup.
- Phase 1 crate `tpt-mat-crystallography` providing `CrystalStructure`,
  `MillerIndex`, `LatticeParameters`, `SlipSystem`, and
  `CrystalStructure::slip_systems()` for **FCC** (12: {111}⟨110⟩),
  **BCC** (12: {110}⟨111⟩), and **HCP** (basal, prismatic, pyramidal ⟨a⟩,
  pyramidal ⟨c+a⟩).
- Phase 1 crate `tpt-mat-wasm` scaffolding only; full bindings arrive in
  Phase 8.
- Verification tests:
  - FCC has exactly 12 slip systems.
  - BCC has exactly 12 slip systems.
  - HCP has 3 basal + 3 prismatic + 6 pyramidal ⟨a⟩ + 12 pyramidal ⟨c+a⟩.
  - Schmid tensor `P = sym(s ⊗ n)` is symmetric.
  - `resolved_shear_stress(σ, s, n) == σ · P · σ̂` (component form).
- CI workflows: `ci.yml`, `license.yml`, `benchmark.yml`, `docs.yml`,
  `release.yml`.
- Issue and pull request templates.

### Changed

- _none._

### Deprecated

- _none._

### Removed

- _none._

### Fixed

- _none._

### Security

- _none._

## [0.0.0] — _not yet released_

Initial scaffold commit. No crates yet.

---

**Release cadence:** SemVer, every 6 weeks. The
[`.github/workflows/release.yml`](./.github/workflows/release.yml) workflow
automates cutting releases and publishing to crates.io.

**Maintainer:** TPT Solutions &lt;[email protected]&gt;