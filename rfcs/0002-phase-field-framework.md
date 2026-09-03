# RFC 0002 — Phase-Field Framework

| | |
|---|---|
| **Status** | Accepted |
| **Author** | TPT Solutions |
| **Created** | 2026-09-02 |
| **Phase** | 3 |

## Summary

Adopt the standard Allen–Cahn / Cahn–Hilliard / Kobayashi
phase-field framework as the modelling backbone for
`tpt-mat-phase-field`, `tpt-mat-grain-growth`, and
`tpt-mat-solidification`, and wire the grid-based solvers through the
shared `tpt-science::grid` substrate.

## Scope

- **Phase-field models**: Allen–Cahn (non-conserved order parameter),
  Cahn–Hilliard (conserved concentration), Kobayashi
  (Allen–Cahn + latent-heat thermal coupling), multi-phase Allen–Cahn
  for grain growth.
- **Bulk free energies**: symmetric double-well, 6th-order polynomial,
  regular-solution `Ωc(1-c) + RT [c ln c + (1-c) ln(1-c)]`.
- **Discretisation**: 2nd-order central finite differences on regular
  grids (`Grid1D`, `Grid2D`, `Grid3D` from `tpt-science`), Neumann
  zero-flux boundaries, forward-Euler time integration.
- **Phase-field grains / dendrites**: dedicated crates
  (`tpt-mat-grain-growth`, `tpt-mat-solidification`) that wrap the
  phase-field solver with grain-boundary mobility (Read–Shockley
  low-angle, fixed HAGB mobility) and 4-fold / 6-fold / isotropic
  anisotropy respectively.

## Out of scope (deferred to later phases)

- Spectral / FFT-based solvers (planned for Phase 4 RVE).
- Adaptive mesh refinement.
- 3D dendrite simulations.

## Test coverage

- Allen–Cahn: free energy monotonically decreases.
- Cahn-Hilliard: mean concentration conserved.
- Multi-grain growth: number of grains preserved under flat-boundary
  evolution.
- Solidification: solid fraction grows from a seed under undercooling.

## References

- Allen & Cahn, 1979. *A microscopic theory for domain wall motion*.
- Cahn & Hilliard, 1958. *Free energy of a nonuniform system*.
- Kobayashi, 1993. *Modeling and numerical simulations of dendritic
  crystal growth*. Physica D 63, 410–423.
- Moelans et al., 2008. *Phase-field analysis of grain growth in a
  polycrystalline material*. Acta Mater. 56, 3471–3478.