# RFC 0015 — Hydrogen Transport with Trapping + Embrittlement Criteria

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-hydrogen-embrittlement`.

---

## 1. Summary

Hydrogen embrittlement is a degradation mechanism that bridges
`tpt-mat-diffusion` (transport) with `tpt-mat-fracture` /
`tpt-mat-fatigue-micro` (failure).  This RFC documents the
`tpt-mat-hydrogen-embrittlement` crate.

## 2. Models

- **Oriani local equilibrium** — effective diffusivity
  `D_eff = D_L / (1 + ∂C_T/∂C_L)` with the single-trap
  approximation `∂C_T/∂C_L = N_T K / (1 + K C_L)²`.
- **McNabb–Foster kinetics** — kinetic trapping rate
  `∂C_T/∂t = k C_L (N_T − C_T) − p C_T`.
- **Stress-driven uphill flux** — hydrostatic-stress term
  `J_stress = -D C V_H / (RT) ∇σ_h`.
- **HEDE / HELP thresholds** — critical-lattice-decohesion
  (HEDE) and HELP local-plasticity indicators.
- **Susceptibility index** — `C × T` normalised by the critical
  values.

## 3. Verification

- Zero-trap case recovers `D_eff = D_L`.
- Trapping retards the effective diffusivity by the predicted
  factor.
- McNabb–Foster rate is zero at local equilibrium.
- Stress-driven flux reverses sign with the hydrostatic gradient.

## 4. References

- Oriani, R. A. (1970).  *Acta Metall.* 18, 147–157.
- McNabb, A., & Foster, P. K. (1963).  *Trans. AIME* 227, 618.
- Troiano, A. R. (1960).  *Trans. ASM* 52, 54.
- Beachem, C. D. (1972).  *Metall. Trans.* 3, 441.
- Sofronis, P., & Birnbaum, H. K. (1995).  *J. Mech. Phys.
  Solids* 43, 49–90.