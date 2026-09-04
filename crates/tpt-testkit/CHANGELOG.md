# Changelog — tpt-testkit

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Updated README.md to match the workspace-wide documentation standard.

## [0.1.0] — 2026-09-02

### Added
- `assert_close`, `assert_close_loose`, and `within` floating-point
  comparison helpers used by every `tpt-mat-*` crate's unit tests.
- MIT-clean replacement for `approx` so the workspace satisfies
  `cargo deny check licenses` without exceptions.

[Unreleased]: https://github.com/tpt-solutions/tpt-materials/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-materials/releases/tag/tpt-testkit-v0.1.0