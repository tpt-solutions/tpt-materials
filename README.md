# tpt-materials

> **Micro-Scale Physics. Macro-Scale Impact. MIT-licensed.**

`tpt-materials` is a fully open-source, MIT OR Apache-2.0 dual-licensed
computational engine for **micro-scale material physics**: crystal plasticity,
phase-field modeling, diffusion, micro-mechanics, homogenization, and
degradation. It is the "micro-to-macro bridge" that produces the effective
material properties, degradation curves, and fatigue limits consumed by the
rest of the TPT Solutions stack (`tpt-energy`, `tpt-transport`,
`tpt-electronics`, `tpt-medical`).

| | |
|---|---|
| **Organization** | TPT Solutions |
| **Crate prefix** | `tpt-mat-*` |
| **License** | MIT OR Apache-2.0 (dual) |
| **MSRV** | Rust 1.75 |
| **Release cadence** | SemVer, every 6 weeks |
| **Status** | Phase 1 / Foundation (active development) |

## Why this exists

Materials simulation is trapped between academic GPL code
(DAMASK, MOOSE, FEniCS, LAMMPS) and commercial black boxes
(Abaqus UMATs, ANSYS Material Designer, Thermo-Calc). Neither is suitable
for the commercial R&D teams at steel mills, battery makers, and aerospace
suppliers who need **device-grade output** (degradation curves, fatigue
limits, homogenized stiffness) **from micro-scale physics**.

`tpt-materials` is that bridge:

- **Pure Rust**, compiles to native and to **WebAssembly** (run RVE
  homogenization or phase-field directly in the browser).
- **MIT OR Apache-2.0**. No GPL traps, no academic lock-in.
- **100 % open source**. Every crate, every constitutive model, every
  solver is in this repo.

## Crates (current)

| Crate | Phase | Description |
|---|---|---|
| `tpt-math-linalg-fixed` | 1 | Fixed-size vector/matrix/tensor types (Vec3, Mat3, SymMat3, Vec6, stiffness/Schmid tensors) |
| `tpt-mat-core` | 1 | `MaterialMicrostructure`, `Phase`, `Grain`, `Composition`, `CrystalOrientation` |
| `tpt-mat-constants` | 1 | `PhysicalConstants`, `PeriodicTable` |
| `tpt-mat-crystallography` | 1 | `CrystalStructure`, `MillerIndex`, `SlipSystem`, Schmid tensor, RSS |
| `tpt-mat-wasm` | 1 | WASM bindings (scaffold; full bindings in Phase 8) |

See [`todo.md`](./todo.md) for the full 8-phase roadmap.

## Quick start

```toml
# Cargo.toml
[dependencies]
tpt-mat-core = "0.1"
tpt-mat-crystallography = "0.1"
tpt-mat-constants = "0.1"
```

```rust
use tpt_mAT_crystallography::{CrystalStructure, SlipSystem};
use tpt_math_linalg_fixed::Vec3;

let fcc = CrystalStructure::FCC;
let slips = fcc.slip_systems();
assert_eq!(slips.len(), 12);

// Schmid tensor P = sym(s ⊗ n)
let s = Vec3::new(1.0, 0.0, 0.0);
let n = Vec3::new(0.0, 0.0, 1.0);
let p = SlipSystem::schmid_tensor(s, n);
```

## Building from source

```bash
git clone https://github.com/tpt-solutions/tpt-materials
cd tpt-materials
cargo build --workspace
cargo test  --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
cargo deny check licenses advisories bans
```

## Examples

```bash
cargo run --example fcc_single_crystal_tension
```

(`examples/fcc-single-crystal-tension/` contains Phase 1 data; the full
single-crystal-tension solver lands in Phase 2.)

## Documentation

- `docs/book/` — mdBook user guide (in progress)
- `docs/api/` — rustdoc-generated API reference
- `spec.txt` — full design specification
- `todo.md` — phase-by-phase roadmap with checkboxes
- `rfcs/` — Requests for Comments for new constitutive models and APIs

## Contributing

We welcome contributions. Read [`CONTRIBUTING.md`](./CONTRIBUTING.md) — the
short version is:

1. Fork and branch.
2. Write code **plus** tests.
3. `cargo fmt`, `cargo clippy -- -D warnings`, `cargo test` all clean.
4. `cargo deny check licenses` clean.
5. Sign off your commits (`git commit -s`) per the [DCO](https://developercertificate.org/).
6. Open a PR. New constitutive models require an RFC discussion and **two
   maintainer approvals** before merge.

## Governance

- **License:** MIT OR Apache-2.0 (dual). See [`LICENSE-MIT`](./LICENSE-MIT) and
  [`LICENSE-APACHE`](./LICENSE-APACHE).
- **Contributions:** MIT OR Apache-2.0, accepted under DCO (CLA-free).
- **Trademark:** "TPT Materials" is reserved by TPT Solutions. See
  [`TRADEMARK.md`](./TRADEMARK.md).
- **Process:** Benevolent dictator + RFC process (see `rfcs/`).
- **Roadmap:** Public GitHub Projects board (TBD).
- **Security:** Private disclosure via [`SECURITY.md`](./SECURITY.md).

## Trademark

"TPT Materials" and "TPT Solutions" are trademarks of TPT Solutions and are
not licensed under MIT/Apache-2.0. See [`TRADEMARK.md`](./TRADEMARK.md).

## Code of Conduct

By participating, you agree to abide by our
[`CODE_OF_CONDUCT.md`](./CODE_OF_CONDUCT.md).

## License

```
MIT OR Apache-2.0
```

Copyright (c) 2026 TPT Solutions. Either LICENSE-MIT or LICENSE-APACHE may
apply at your option.