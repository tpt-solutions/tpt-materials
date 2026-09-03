# Licensing policy

`tpt-materials` is released under **`MIT OR Apache-2.0`** (dual), at the
downstream user's option. We additionally guarantee that a downstream
user who chooses **MIT alone** is never forced onto any other license by
our dependency tree.

## Rules

1. **First-party code is dual `MIT OR Apache-2.0`.** Every crate sets
   `license.workspace = true`.

2. **Every dependency must offer MIT (or something strictly more
   permissive).** Acceptable dependency licenses:

   | Category | Examples |
   |---|---|
   | MIT and MIT-family | `MIT`, `MIT-0` |
   | Public-domain-equivalent | `CC0-1.0`, `Unlicense`, `0BSD` |
   | Permissive, MIT-compatible | `BSD-2-Clause`, `BSD-3-Clause`, `ISC`, `Zlib`, `BSL-1.0`, `Unicode-3.0`, `Unicode-DFS-2016` |

   A dependency licensed `MIT OR Apache-2.0` is fine — we consume it under
   MIT.

3. **Disallowed, no exceptions:**
   - `Apache-2.0` **as the only option** (forces Apache patent/notice
     terms onto a would-be MIT consumer)
   - `MPL-2.0` and any other file- or module-level copyleft
   - `GPL-*`, `LGPL-*`, `AGPL-*`, and any reciprocal license
   - Unlicensed / license-unknown crates

4. **Enforcement:** `deny.toml`'s `[licenses].allow` list omits
   `Apache-2.0`, so `cargo deny check licenses` fails on any Apache-only
   or copyleft crate anywhere in the tree. Run it in CI and locally
   before every release.

5. **Clean-room reimplementation.** Where the only mature crate for a
   capability is GPL/Apache-only (e.g. `approx` → `tpt-testkit`), we
   write a minimal MIT implementation from scratch rather than take the
   dependency. Such crates carry a header note pointing back to this
   file.

## Current first-party substitutions

| Instead of | We use | Why |
|---|---|---|
| `approx` (Apache-2.0) | `tpt-testkit` | test-only float asserts; trivial to reimplement |
| `nalgebra` (Apache-2.0) | `tpt-math-linalg-fixed` | fixed-size vector/matrix/tensor math |

## Third-party runtime dependencies (allowed)

`serde`, `serde_json`, `thiserror`, `num-traits`, `wasm-bindgen`,
`js-sys` — all `MIT OR Apache-2.0`, consumed under MIT.
