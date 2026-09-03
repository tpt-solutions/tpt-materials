# RFC 0003 — Diffusion, Phase-Transform Kinetics, and CALPHAD Framework

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-diffusion`, `tpt-mat-phase-transform`, `tpt-mat-calphad`.

---

## 1. Summary

Phase 4 of the `tpt-materials` roadmap introduces three new crates
that capture the macroscale thermodynamics and kinetics of
solid-state phase transformations and substitutional diffusion:

- **`tpt-mat-diffusion`** — Fick's second-law solvers on regular
  grids (single- and multi-component), Arrhenius temperature
  dependence, and Fisher-style grain-boundary diffusion.
- **`tpt-mat-phase-transform`** — Analytical isothermal
  (Avrami/JMAK) and diffusionless (Koistinen–Marburger)
  transformation kinetics, with a piecewise-linear thermal-history
  driver that integrates the analytical models using the Scheil
  additivity rule.
- **`tpt-mat-calphad`** — Building blocks of CALPHAD-style
  thermodynamic modelling: Redlich–Kister excess polynomials,
  end-member + ideal-mixing + excess Gibbs energies, Muggianu
  sublattice models, and a grid-search common-tangent construction
  for binary two-phase equilibria.

This RFC ratifies the public APIs of these three crates as the
canonical Phase 4 surface and records the explicit non-goals
(deliberately deferred to Phase 7–8).

## 2. Motivation

Crystal-plasticity (Phase 2) and phase-field (Phase 3) simulations
need thermodynamic driving forces and diffusional mobilities that
come from an *assessed* thermodynamic database.  Without a CALPHAD
crate that can ingest (even minimally-parameterised) Redlich–Kister
or sublattice descriptions, every downstream crate ends up
re-implementing the same energy/mobility plumbing.

Phase 4 also enables the homogenisation stack (Phase 5) to compute
on-the-fly the diffusivities the RVE solver should evolve over an
    -scale temperature history, and gives the manufacturing-prediction
stack (Phase 7) a calibrated `TTT` / `CCT` diagram source for heat-
treatment simulation.

## 3. Scope

### In scope

- **Single-component diffusion:** Forward-Euler finite-difference
  solver on `tpt_science::Grid2D` with explicit CFL stability
  reporting (`Δt ≤ Δx² / (4 D)`).
- **Multi-component diffusion:** Independent species, each with its
  own scalar `D_i`; full diagonal `M` matrix (no cross-terms).
- **Arrhenius diffusivity:** `D(T) = D_0 · exp(−Q / (R T))` with the
  temperature derivative exposed for sensitivity studies.
- **Grain-boundary diffusion:** Fisher-form `D_eff = D_l + (π δ / L)
  D_gb` (Regime A only; B/C deferred).
- **Avrami / JMAK:** `f(t, T) = 1 − exp(−k(T) t^n)` with
  Arrhenius-temperature-dependent rate constant `k(T)`.  TTT diagram
  generation is a single forward call.
- **Koistinen–Marburger:** `f(T) = 1 − exp(−α (M_s − T))` for
  `T ≤ M_s`, with `retained_at_m_start` cap.
- **Transformation driver:** Isothermal holds + continuous-cooling
  segments (16 inner Scheil sub-steps per segment) with combined
  `f_total = 1 − (1 − f_av)(1 − f_km)`.
- **Redlich–Kister polynomial:** `^E G(x_B) = Σ L_ν (x_A − x_B)^ν`,
  with closed-form derivative for the common-tangent search.
- **Gibbs energy:** End-member + ideal mixing (configurational
  entropy `RT(x_B ln x_B + x_A ln x_A)`) + Redlich–Kister excess.
- **Sublattice:** Multi-sublattice Muggianu configurational Gibbs
  energy (`−T · S_cfg`).
- **Phase diagram:** Binary `T-x` boundary via grid-search
  common-tangent construction (1024 grid points per phase).  Returns
  `PhaseBoundary { t_k, x_alpha, x_beta }` per temperature.

### Out of scope (Phase 4 deferred)

- Cross-term (`L_{ij}`, `i ≠ j`) multi-component mobilities (will land
  in Phase 7 once CALPHAD mobility databases are wired in).
- TDB / PAC parser; CALPHAD inputs are constructed in code for now.
- Hillert–Staffansson diffusion-controlled moving-boundary
  transformations (Phase 7, alongside homogenisation).
- Bayesian / inverse CALPHAD assessment of unknown parameters (Phase 8).

## 4. API surface

### `tpt-mat-diffusion`

```rust
pub use arrhenius::{ArrheniusDiffusivity, ArrheniusParams};
pub use grain_boundary::{GrainBoundaryDiffusion, GBDiffusivity};
pub use multicomponent::{MultiComponentDiffusionSolver, MultiComponentError};
pub use scalar::{DiffusionSolver, DiffusionSolverError, DiffusionStepResult};
```

`DiffusionSolver::step()` is forward-Euler and reports
`total_mass` (Neumann-mass-conserved to ~1e-6 relative error per step
for fields with negligible boundary contribution) and
`max_abs_change`.

### `tpt-mat-phase-transform`

```rust
pub use avrami::{AvramiModel, AvramiParams, TttPoint};
pub use koistinen::{KoistinenMarburger, KMParams};
pub use solver::{
    ThermalHistoryStep, TransformationError, TransformationSolver,
    TransformationState,
};
```

The combined `TransformationSolver` integrates the Avrami fraction
using the additive form (`−ln(1 − f) = k t^n`) so consecutive
isothermal / cooling segments compose correctly.  Martensite is
applied as an instantaneous jump at each temperature boundary.

### `tpt-mat-calphad`

```rust
pub use gibbs::{GibbsEnergyModel, IdealMixing, EndMember};
pub use phase_diagram::{
    PhaseBoundary, PhaseDiagram, TwoPhaseEquilibrium, two_phase_equilibrium,
};
pub use redlich_kister::{RedlichKister, RedlichKisterParams};
pub use sublattice::{Sublattice, SublatticeConfig, SublatticeModel, SublatticeSpecies};
```

The common-tangent search in `two_phase_equilibrium` is a brute-
force grid method — adequate for phase-diagram sketching but not a
substitute for a full thermodynamic-equilibrium solver.  It is
called out as such in the doc-comment.

## 5. Verification

The crates ship with unit tests asserting:

| Property | Where |
|---|---|
| `∂c / ∂t = D ∇²c` is mass-conserving to ~1e-6 relative under Neumann | `tpt-mat-diffusion/src/scalar.rs` |
| Multi-component solver conserves per-species mass | `tpt-mat-diffusion/src/multicomponent.rs` |
| Arrhenius `D` increases with `T`; derivative is positive | `tpt-mat-diffusion/src/arrhenius.rs` |
| GB diffusion > lattice diffusion; recovers lattice for large grains | `tpt-mat-diffusion/src/grain_boundary.rs` |
| Avrami `f(t, T)` is monotone, in `[0, 1]`, with closed-form `t(f, T)` | `tpt-mat-phase-transform/src/avrami.rs` |
| Koistinen–Marburger is monotone on cooling, 0 above `M_s` | `tpt-mat-phase-transform/src/koistinen.rs` |
| Combined solver integrates an isothermal hold + cooling ramp | `tpt-mat-phase-transform/src/solver.rs` |
| Redlich–Kister derivative matches finite differences | `tpt-mat-calphad/src/redlich_kister.rs` |
| Ideal mixing is 0 at endpoints, negative inside | `tpt-mat-calphad/src/gibbs.rs` |
| Sublattice configurational entropy is non-positive | `tpt-mat-calphad/src/sublattice.rs` |
| `two_phase_equilibrium` returns `x_α < x_β` for non-convex `G(x)` | `tpt-mat-calphad/src/phase_diagram.rs` |

End-to-end examples:

- `examples/diffusion-carbon-steel`: 1-hour carburisation of a steel
  slab using Arrhenius `D(T)` and an explicit-Euler Fick's-second-law
  solver on a 100-cell 1-mm grid.
- `examples/phase-transform-jmak`: TTT diagram + isothermal hold +
  continuous-cooling transformation through `M_s`.
- `examples/calphad-binary-phase-diagram`: binary A-B `T-x` phase
  boundary from a Redlich–Kister-driven common-tangent construction.

## 6. Migration / breaking changes

None — this RFC introduces three new crates and three new examples.
No existing public API is modified.

## 7. Open questions

- Whether to expose a CALPHAD `TDB` parser (vs requiring all inputs
  in code).  Decision: defer until at least one published
  TDB-to-Rust workflow is demonstrated externally.
- Whether the multi-component diffusion crate should consume CALPHAD
  mobilities directly (current API: caller supplies per-species `D_i`
  explicitly).  Decision: expose a `pub fn effective_d_from_calphad`
  in a follow-up RFC once Phase 5 homogenisation lands.

## 8. Alternatives considered

- **Linking to PyCalphad or Thermochimica** for full CALPHAD
  support.  Rejected for now — adds a Python dependency that
  conflicts with the Rust-native design and prevents `no-std`
  deployment.
- **Implementing CALPHAD in C/C and** as bindings.  Rejected — the
  surface area is small enough that native Rust is the right
  long-term investment.
- **Phase-field diffusion coupling inside `tpt-mat-phase-field`**.
  Rejected as in-scope for `tpt-mat-diffusion` so the phase-field
  crate can stay focused on the Allen-Cahn / Cahn-Hilliard dynamics.