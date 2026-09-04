# tpt-mat-wasm

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-wasm.svg)](https://crates.io/crates/tpt-mat-wasm)
[![Documentation](https://docs.rs/tpt-mat-wasm/badge.svg)](https://docs.rs/tpt-mat-wasm)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> WebAssembly bindings for `tpt-materials` (scaffold; full bindings
> queued for the Phase 8 ecosystem / informatics work).

`tpt-mat-wasm` exposes the Rust micro-mechanics solvers to a browser
runtime via `wasm-bindgen`.  Currently published (for `tpt-materials`
`full` feature builds):

- `WasmRveSolver` — Voigt / Reuss / self-consistent homogenisation of a
  user-supplied RVE.
- `WasmPhaseField` — `step` and `get_order_parameter` accessors.

Full WASM coverage of `tpt-mat-crystal-plasticity`,
`tpt-mat-fatigue-micro`, and `tpt-mat-hydrogen-embrittlement` is queued
for the Phase 8 follow-up.

## Building

```bash
# requires the wasm32-unknown-unknown target
cargo build -p tpt-mat-wasm --target wasm32-unknown-unknown --release
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).