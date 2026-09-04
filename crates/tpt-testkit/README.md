# tpt-testkit

[![Crates.io](https://img.shields.io/crates/v/tpt-testkit.svg)](https://crates.io/crates/tpt-testkit)
[![Documentation](https://docs.rs/tpt-testkit/badge.svg)](https://docs.rs/tpt-testkit)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Minimal from-scratch floating-point test assertions for the tpt-materials
> workspace (MIT-clean approx replacement).

`tpt-testkit` provides the assertion helpers that every other `tpt-mat-*` crate
depends on in its `dev-dependencies`. It is intentionally tiny: a handful of
`assert_close`, `assert_close_loose`, and similar floating-point comparison
functions plus a `within` extension trait. No third-party testing libraries are
used, which keeps the workspace's MIT/Apache-2.0 dependency chain short and
satisfies `cargo deny check licenses` without exceptions.

## When to use it

This crate exists only to be referenced from other workspace crates as
`dev-dependencies = { tpt-testkit = { workspace = true } }`. It is not a
general-purpose testing library; for new code prefer
[`approx`](https://crates.io/crates/approx) or the standard
[`float_eq`](https://crates.io/crates/float_eq) crates.

## API overview

| Function | Purpose |
|---|---|
| `assert_close(a, b)` | Strict absolute / relative comparison |
| `assert_close_loose(a, b)` | Looser tolerance for iterative solvers |
| `within(a, b, eps)` | Boolean comparison returning `bool` |

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).