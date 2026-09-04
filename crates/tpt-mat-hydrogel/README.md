# tpt-mat-hydrogel

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-hydrogel.svg)](https://crates.io/crates/tpt-mat-hydrogel)
[![Documentation](https://docs.rs/tpt-mat-hydrogel/badge.svg)](https://docs.rs/tpt-mat-hydrogel)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Hydrogel swelling equilibrium and kinetics — Flory–Rehner
> mixing + elastic osmotic balance, Fickian solvent uptake on a
> `tpt-science` grid.

`tpt-mat-hydrogel` is the soft-matter hydrogel pair to
`tpt-mat-polymer`:

- Flory–Rehner swelling equilibrium (mixing + elastic osmotic
  pressure balance, root-find over polymer volume fraction).
- Poroelastic swelling kinetics — Fickian solvent uptake on a
  `tpt-science` grid (1-D series solution for slab geometry).
- `equilibrium_swelling_ratio(chi, crosslink_density) -> f64`.

## Example

```rust
use tpt_mat_hydrogel::equilibrium_swelling_ratio;

let q = equilibrium_swelling_ratio(0.4, 1e-3);
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).