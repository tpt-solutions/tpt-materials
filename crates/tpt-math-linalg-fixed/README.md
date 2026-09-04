# tpt-math-linalg-fixed

[![Crates.io](https://img.shields.io/crates/v/tpt-math-linalg-fixed.svg)](https://crates.io/crates/tpt-math-linalg-fixed)
[![Documentation](https://docs.rs/tpt-math-linalg-fixed/badge.svg)](https://docs.rs/tpt-math-linalg-fixed)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Fixed-size vector / matrix / tensor types used across tpt-materials.

`tpt-math-linalg-fixed` is the substrate numerical crate of the workspace. It
provides:

- `Vec2`, `Vec3` — small Euclidean vectors with `+`, `-`, `*`, `dot`, `cross`,
  `norm`, `normalize`.
- `Mat3` — 3×3 matrix with explicit column-major storage, transpose, inverse,
  trace, determinant, `sym` / `skew` decomposition, and product with `Vec3`.
- `SymmetricFourthOrder` — 6×6 Voigt-notation elastic stiffness with
  cubic / hexagonal / isotropic constructors and Mandel / minor symmetry.
- Helper utilities for 4th-order tensor → 6×6 conversion used by
  `tpt-mat-crystallography` and `tpt-mat-crystal-plasticity`.

The crate is `no_std`-friendly (depends only on `num-traits`), has zero
unsafe code (`unsafe_code = "forbid"` workspace-wide), and is intended for
inlining into hot loops without monomorphisation overhead.

## API overview

```rust
use tpt_math_linalg_fixed::Mat3;

let a = Mat3::identity();
let b = Mat3::from_columns([Vec3::new(1.0, 0.0, 0.0),
                            Vec3::new(0.0, 2.0, 0.0),
                            Vec3::new(0.0, 0.0, 3.0)]);
let c = a * b;          // matrix product
let v = Vec3::new(1.0, 1.0, 1.0);
let w = c * v;          // matrix-vector product
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).