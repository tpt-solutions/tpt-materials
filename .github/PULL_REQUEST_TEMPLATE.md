## Summary

<!-- One-paragraph description of the change. -->

## Related issue / RFC

<!-- Link to the issue or `rfcs/NNNN-*.md` this PR addresses. -->

## Type of change

- [ ] Bug fix (non-breaking)
- [ ] New feature (non-breaking)
- [ ] Breaking change (requires RFC + CHANGELOG entry)
- [ ] Documentation
- [ ] Refactor (no functional change)
- [ ] Performance

## Checklist

- [ ] I have signed off every commit (`git commit -s`)
- [ ] I have read [`CONTRIBUTING.md`](../CONTRIBUTING.md)
- [ ] `cargo fmt --all -- --check` passes
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` passes
- [ ] `cargo test --workspace` passes
- [ ] `cargo deny check` passes
- [ ] I have added tests for new functionality
- [ ] I have updated [`CHANGELOG.md`](../CHANGELOG.md) under "Unreleased"
- [ ] I have updated rustdoc comments on public API
- [ ] If this is a new constitutive model, the linked RFC has **two maintainer approvals**

## Test plan

<!-- How was this verified? Include command output, before/after numbers, etc. -->