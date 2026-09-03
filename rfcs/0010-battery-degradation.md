# RFC 0010 — Battery Degradation, Diffusion-Induced Stress and Capacity-Fade Pipeline

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-battery`.

---

## 1. Summary

This RFC ships [`tpt-mat-battery`], the Phase-7 energy-materials
crate that produces the cross-repo `DegradationCurve` consumed by
`tpt-energy` (spec §6).

## 2. Active-material model

An [`ActiveMaterial`](tpt-mat-battery::ActiveMaterial) describes one
electrode particle (chemistry, particle radius, diffusion
coefficient, partial molar volume, elastic constants, fracture
stress, max concentration).  Built-in defaults are bundled for
NMC811 / NMC622 / LFP / NCA / graphite / silicon via
[`BatteryChemistry`](tpt-mat-battery::BatteryChemistry).

## 3. Solid-state Li diffusion

[`simulate_diffusion`](tpt-mat-battery::simulate_diffusion) solves
`∂c/∂t = D/r² ∂/∂r (r² ∂c/∂r)` on a 1-D radial grid (forward Euler
with sub-step CFL reduction).  The C-rate is converted to a
surface concentration Dirichlet value via the applied current and
the bulk solid-state diffusivity.

## 4. Degradation mechanisms

Four [`DegradationMechanism`](tpt-mat-battery::DegradationMechanism)
variants are wired through [`DegradationCurve`][`tpt-mat-battery::DegradationCurve`]
construction:

- `SeiGrowth` — Arrhenius `δ̇ = k · exp(−Q/RT)`;
- `ParticleCracking` — instant loss when hoop stress exceeds
  `σ_c`;
- `LithiumPlating` — flag when surface overpotential falls below
  `0 V` vs Li/Li⁺;
- `TransitionMetalDissolution` — first-order decay.

## 5. Capacity-fade curve

[`capacity_fade_curve`](tpt-mat-battery::capacity_fade_curve) returns
a [`DegradationCurve`](tpt-mat-battery::DegradationCurve) over
`num_cycles` cycles: cycles → capacity-retention (monotonically
decreasing) + resistance-growth curves.  At low C-rate the sqrt(t)
SEI limit dominates (`t^{1/2}` concavity).

## 6. Example

[`examples/battery-electrode-degradation`] produces an NMC811
`DegradationCurve` over 10 000 cycles with cycle-by-cycle
capacity-retention and resistance-growth milestones.

## 7. References

- Newman, J. (1963).  *The Electrochemical Reaction.*  PhD thesis,
  UC Berkeley.
- Peled, E. (1979).  *J. Electrochem. Soc.* 126(12), 2047–2051.
- Christensen, J., & Newman, J. (2004).  *J. Solid State
  Electrochem.* 10, 293–319.
- Doyle, M., Fuller, T. F., & Newman, J. (1993).  *J. Electrochem.
  Soc.* 140(6), 1526–1533.