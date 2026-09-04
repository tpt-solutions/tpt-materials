# Changelog — tpt-materials

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- Umbrella / facade crate re-exporting every `tpt-mat-*` domain crate.
- Per-domain cargo features (`crystallography`, `crystal-plasticity`,
  `phase-field`, `homogenization`, `rve`, `battery`, …) and the `full`
  meta-feature.
- Spec §6 / §13 doc-tests for `use tpt_materials::…` import snippets.

## [0.1.1] — 2026-09-04

### Added
- `inverse` feature re-exporting `tpt-mat-inverse`.
- Public spec §6 / §13 doc-test wrappers.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-materials-v0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-materials-v0.1.0