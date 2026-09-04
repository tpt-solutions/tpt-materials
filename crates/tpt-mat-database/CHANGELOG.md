# Changelog — tpt-mat-database

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-03

### Added
- `MaterialRecord` (composition / mechanical / thermal / electrical).
- `DataSource` provenance enum — `Textbook`, `Standard`, `Nist`,
  `Lab`, `Datasheet`.
- `PropertyQuery` — property-range queries.
- `MaterialsDatabase::load_builtin()` — bundled MIT-clean property set.
- `search_by_property` returning the filtered record list.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-mat-database-v0.1.0