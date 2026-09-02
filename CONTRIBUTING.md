# Contributing to tpt-materials

Thank you for your interest in contributing to `tpt-materials`. This project
follows a lightweight, RFC-driven open-source process.

## Code of Conduct

By participating you agree to abide by the
[`CODE_OF_CONDUCT.md`](./CODE_OF_CONDUCT.md). Maintainers may remove
contributions or block contributors who violate it.

## Licensing & DCO

`tpt-materials` is dual-licensed **MIT OR Apache-2.0**. We do **not** use a
CLA; instead every commit must be signed off under the
[Developer Certificate of Origin (DCO)](https://developercertificate.org/).

Add a `Signed-off-by:` line to every commit:

```bash
git commit -s -m "Add BCC slip systems"
# produces:
# Add BCC slip systems
#
# Signed-off-by: Your Name <[email protected]>
```

By signing off you affirm that you wrote the contribution yourself (or have
permission to pass it on under the project license) and that you agree to
the DCO terms.

## Workflow

1. **Fork** the repository and create a topic branch:
   ```bash
   git checkout -b feat/my-constitutive-model
   ```
2. **Code + tests.** All new functionality comes with tests. New constitutive
   models are expected to include either an analytical reference case or a
   golden-file comparison in `test-data/golden/`.
3. **Local gates** — all of these must pass before opening a PR:
   ```bash
   cargo fmt --all -- --check
   cargo clippy --workspace --all-targets -- -D warnings
   cargo test  --workspace
   cargo deny check licenses advisories bans
   ```
4. **Open a Pull Request** using the
   [PR template](./.github/PULL_REQUEST_TEMPLATE.md). The CI workflow
   re-runs the same gates plus cross-platform builds.
5. **DCO sign-off** — every commit must carry `Signed-off-by:`. The CI bot
   will comment and block the merge if sign-off is missing.
6. **Review** — two maintainer approvals are required for merge. New
   constitutive models additionally require an accepted **RFC** in `rfcs/`.

## When do you need an RFC?

Open an RFC in `rfcs/` (using `rfcs/0000-template.md` as a starting point)
**before** opening a PR that:

- Adds a new constitutive model (Voce, power-law, GTN, Arruda-Boyce, etc.).
- Changes a public API in a non-backwards-compatible way.
- Adds a new top-level crate to the workspace.
- Changes the on-disk format of golden test data.

RFCs are discussed openly in the PR that introduces them. After **two
maintainer approvals** the RFC is merged and the implementing PR can move
forward.

Smaller changes — bug fixes, performance improvements, refactors that
preserve public API, documentation — do not need an RFC.

## Style

- `cargo fmt` with the workspace `rustfmt.toml`.
- `cargo clippy -- -D warnings`. No new warnings on existing code.
- MSRV is **Rust 1.75**. Do not require newer features without discussion.
- Prefer `f64` for material quantities; `f32` is acceptable only in
  performance-critical inner loops with explicit rationale in a comment.
- Prefer strong types — `Vec3`, `Mat3`, `SymMat3` from
  `tpt-math-linalg-fixed` — over bare `[f64; 3]` / `[f64; 9]`.
- Use `#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]` on
  configuration and result types unless there is a good reason not to.
- Errors via `thiserror` enums; no `unwrap`/`expect` in library code (tests
  and `examples/` are fine).
- **No comments unless they add something the code doesn't say.** If the
  code is self-explanatory, leave it.

## Tests

- Unit tests live next to the code they test (`#[cfg(test)] mod tests`).
- Golden tests live in `tests/golden.rs` and compare against JSON files in
  `test-data/golden/`. Golden outputs are committed; CI fails on drift.
- Property-based / fuzz tests are encouraged for tensor math.
- Verification tests (e.g. *"FCC has exactly 12 slip systems"*) are
  first-class tests, not documentation.

## Commit messages

We follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat(crystallography): add HCP pyramidal <c+a> slip systems
fix(core): clamp Euler angle range in CrystalOrientation::euler_bunge
docs(rfcs): publish 0001-crystal-plasticity-fem.md
test(core): add property test for composition sum invariant
chore: bump nalgebra to 0.33
```

Include the `Signed-off-by:` trailer on every commit (use `git commit -s`).

## Reporting bugs & requesting features

Open an issue using the
[bug report](./.github/ISSUE_TEMPLATE/bug_report.md) or
[feature request](./.github/ISSUE_TEMPLATE/feature_request.md) template.
For security issues, **do not** open a public issue — see
[`SECURITY.md`](./SECURITY.md).

## Release process

Maintainers cut releases on a **6-week SemVer cadence**. The
[`.github/workflows/release.yml`](./.github/workflows/release.yml) workflow
publishes to crates.io and creates a GitHub release with auto-generated
notes. See [`CHANGELOG.md`](./CHANGELOG.md) for prior releases.

## Getting help

- Discussions: GitHub Discussions (TBD).
- Discord / Matrix: TBD.
- Maintainer email: [email protected].