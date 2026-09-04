# tpt-mat-database

[![Crates.io](https://img.shields.io/crates/v/tpt-mat-database.svg)](https://crates.io/crates/tpt-mat-database)
[![Documentation](https://docs.rs/tpt-mat-database/badge.svg)](https://docs.rs/tpt-mat-database)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](../LICENSE-MIT)

> Bundled materials database with property-range queries and
> ASTM / ISO / NIST data-source provenance.

`tpt-mat-database` is the on-disk material property registry used by
the consumer crates and the `examples/database-material-search`
example:

- `MaterialRecord { name, composition, mechanical, thermal, electrical, sources }`.
- `MaterialsDatabase { materials }` + `load_builtin()` — bundled
  MIT-clean property set.
- `search_by_property(PropertyQuery) -> Vec<MaterialRecord>` —
  property-range queries (e.g. yield strength ≥ 200 MPa, CTE ≤ 12e⁻⁶
  K⁻¹).
- `DataSource` enum — `Textbook`, `Standard`, `Nist`, `Lab`,
  `Datasheet` provenance for spec §9 traceability.

## Example

```rust
use tpt_mat_database::{MaterialsDatabase, PropertyQuery};

let db = MaterialsDatabase::load_builtin();
let hits = db.search_by_property(PropertyQuery::yield_strength_mpa_in(200.0, 1000.0));
```

## License

Licensed under either of [MIT](../LICENSE-MIT) or
[Apache-2.0](../LICENSE-APACHE) at your option.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).