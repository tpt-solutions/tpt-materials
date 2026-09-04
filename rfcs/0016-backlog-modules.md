# RFC 0016 — Backlog modules (oxidation, GBs, recrystallization, stereology, inverse)

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-corrosion::oxidation`,
`tpt-mat-grain-growth::{interfaces, recrystallization, stereology}`,
`tpt-mat-inverse` (new crate).

---

## 1. Summary

This RFC documents the five modules/crate listed in the
`todo.md` "Backlog" section, all of which were deliberately
deferred after the spec's 10-phase rollout.  They are documented
together because they form a coherent cluster of
**classical / analytical** microstructure-quantification and
calibration utilities that complement (rather than extend) the
main spec-driven crates.

## 2. High-temperature oxidation (`tpt-mat-corrosion::oxidation`)

The corrosion crate is renamed thematically to
"environmental degradation" (electrochemical corrosion + high-
temperature oxidation).  The `oxidation` module provides Wagner
parabolic scale-growth kinetics:

- `ParabolicRateConstant` — Arrhenius `k_p(T) = k_p^0 exp(-Q/RT)`.
- `ScaleGrowth` — driver: thickness `x(t) = sqrt(2 k_p t + x_0^2)`,
  instantaneous rate `dx/dt = k_p / x`, mass gain via the oxide
  stoichiometry.
- `DopingEffect` — Wagner–Hauffe multiplicative correction
  (p-type oxide + acceptor cations / n-type oxide + donor
  cations accelerate growth).
- `BreakawayCriterion` — predicate over critical thickness
  and/or in-plane strain that triggers the parabolic-to-linear
  transition.
- `LinearBreakawayRate` — linear `dx/dt = k_l` once the scale
  fails.
- `OxidationModel` — combined Evans picture: parabolic until the
  breakaway criterion fires (bisection finds `t_b`), then linear
  from there.

References: Wagner (1933) Z. phys. Chem. B21 25; Payer et al.
(1980) "Oxidation of Metals"; Evans (1947) *An Introduction to
Metallic Corrosion*.

## 3. Interfaces / grain boundaries (`tpt-mat-grain-growth::interfaces`)

Consolidates GB-property models and adds:

- `GrainBoundaryEnergy` — Read–Shockley misorientation-dependent
  boundary energy (dual of GB mobility).
- `GBCDEntry` and `GrainBoundaryCharacterDistribution` — the
  `(θ, n_hat, w)` 5-parameter descriptor, with LAB / CSL
  fraction helpers.
- `TripleJunction` — dihedral-angle force balance (Herring 1951)
  via the law of cosines; returns `None` if the Neumann triangle
  is violated.
- `LangmuirMcLean` — equilibrium GB solute segregation
  `X_GB / (1 - X_GB) = (X_bulk / (1 - X_bulk)) exp(-ΔG_seg / RT)`.

The pre-existing `GrainBoundaryMobility` (Read–Shockley mobility
in `mobility.rs`) and `GrainBoundaryDiffusion` (Fisher Regime-A
effective diffusivity in `tpt-mat-diffusion`) remain where they
are; the new module is the home for the consolidated GB-property
*and chemistry* models.

## 4. Recrystallization (`tpt-mat-grain-growth::recrystallization`)

- `JmakRecrystallization` — static RX in JMAK form
  `X(t) = 1 - exp(-k(T) t^n)` with Arrhenius `k(T)`.
- `ZenerHollomon` — `Z = ε̇ exp(Q/RT)` and the Sellars–Tegart
  hyperbolic-sine steady-state stress `σ_s = (1/α) asinh((Z/A)^(1/n))`.
- `DrxKinetics` — Cahn–Hagel-style coupled nucleation + growth
  `X(t) = 1 - exp(-π Ṅ G³ t⁴ / 3)` with Arrhenius `Ṅ(T)` and `G(T)`.

Distinct from `tpt-mat-grain-growth`'s phase-field solver, which
captures curvature-driven *normal grain growth*; RX here is
driven by stored-energy recovery.

## 6. Stereology (`tpt-mat-grain-growth::stereology`)

ASTM-standard microstructure quantification from 2-D section
data:

- `LinearIntercept` — ASTM E112 mean-linear-intercept
  `ℓ̄ = L / N` and grain-size number
  `G = -6.6439 log10(ℓ̄_in_inches)`.
- `AreaFraction2D` — pixel-based `A_A` and binomial counting
  uncertainty `σ = sqrt(A_A (1 − A_A) / N)`.
- `VolumeFraction3D` — voxel-based `V_V`.
- `NumberPerArea` — `N_A` for point-counting analyses.
- `saltykov_size_distribution` — first-pass Saltykov 3-D size-
  class population from a histogram of 2-D section radii.

## 7. Inverse / calibration (`tpt-mat-inverse`, new crate)

Closed-form fits for the most common constitutive laws:

- `arrhenius_fit` — log-linear regression on `(1/T, ln D)` for
  Arrhenius `D_0, Q`.
- `norton_creep_fit` — log-linear regression on `(ln σ, ln ε̇)`
  for Norton `A, n`.
- `voce_fit` — Gauss–Newton fit of the 4-parameter Voce
  hardening law.
- `power_law_sn_fit` — log-linear regression on `(2N, σ_a)` for
  Basquin `σ_f', b`.
- `avrami_fit` — log-linear regression on `(t, f)` for JMAK
  `k, n`.
- `coffin_manson_fit` — log-linear regression on `(2N, Δε_p/2)`
  for Coffin–Manson `ε_f', c`.

For non-linear least-squares problems the crate also ships a
generic `LevenbergMarquardt` driver with forward-difference
Jacobians, used for the LM-recovered tests in this PR.

References: Press et al. (2007) *Numerical Recipes* 3rd ed.
§15.5; ASTM E112-13.

## 8. Example + tests

`examples/backlog-modules-demo` drives every module on synthetic
data:

- Oxidation thickness from 1 min to 1 yr at 900 K.
- Langmuir–McLean β factor across 600 / 800 / 1000 K.
- Triple-junction dihedral angles for asymmetric `γ`.
- LAB and CSL fractions of a representative GBCD.
- Read–Shockley GB energy across `2° – 45°`.
- JMAK RX curve at 900 K + Zener–Hollomon at 1000 K.
- ASTM E112 intercept + area fraction with counting uncertainty.
- Closed-form Arrhenius / Norton / Basquin fits that recover
  the synthetic ground-truth parameters to 1 part in 10⁴.

## 9. Workspace state

- **32 domain crates** (added `tpt-mat-inverse`).
- **26 example binaries** (added `backlog-modules-demo`).
- **16 RFCs** (this one).
- All examples compile and run; full workspace test suite green.