# Security Policy

## Supported versions

`tpt-materials` is in active early development (currently Phase 1). Security
fixes are issued for the latest released version only. There are no
Long-Term Support branches at this stage of the project.

| Version | Supported |
|---|---|
| latest released (`0.x.y`) | ✅ |
| older releases | ❌ |

## Reporting a vulnerability

**Please do not open a public GitHub issue for security vulnerabilities.**

Report privately via one of the following channels:

- **Email:** [email protected] (PGP key on request)
- **GitHub private disclosure:** use
  [the "Report a vulnerability" button](https://github.com/tpt-solutions/tpt-materials/security/advisories/new)
  on the Security tab of the repository.

Please include:

1. A clear description of the vulnerability and its impact.
2. Steps to reproduce, or a minimal proof-of-concept.
3. The affected version(s) and commit hash(es).
4. Your name / handle for the public advisory (or "anonymous" if you prefer).
5. Whether you would like to be credited in the advisory.

You should hear from a maintainer within **3 business days**. We aim to:

- Acknowledge receipt within 3 business days.
- Triage and assign a CVE within 10 business days.
- Coordinate disclosure and a fix timeline with you.

## Coordinated disclosure timeline

We follow a **90-day coordinated disclosure** window. After 90 days from
the initial report (or sooner if a fix is ready) we will publish a GitHub
Security Advisory and a CVE. If a fix is in flight near day 90 we will
negotiate an extension with you.

## Scope

In scope:

- Crates under `crates/*` published from this workspace.
- The WASM bindings under `crates/tpt-mat-wasm`.
- Example programs under `examples/`.
- GitHub Actions workflows under `.github/workflows/`.

Out of scope:

- Vulnerabilities in upstream dependencies that are already disclosed
  upstream — please open an issue and we will bump.
- Theoretical concerns without a demonstrable impact on this codebase.

## Recognized contributors

Security researchers who report valid issues will be credited in the
advisory (unless they prefer anonymity) and listed in a future
`SECURITY-ACKNOWLEDGEMENTS.md`.

## License of this policy

This document is dedicated to the public domain under
[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) so that it
can be freely reused by other projects.