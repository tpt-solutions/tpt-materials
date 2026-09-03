# RFC 0012 — Precipitation Kinetics: Classical Nucleation + KWN + LSW

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-precipitation`.

---

## 1. Summary

Precipitation kinetics is a standard thermo-kinetics topic in
alloy design (Al-Cu, Ni-Al, Fe-Cu-Ni, etc.).  This RFC covers
the three coupled length scales handled in the
`tpt-mat-precipitation` crate.

## 2. Models

- **Classical nucleation theory (Turnbull–Fisher)** — steady-state
  nucleation rate `I = I_0 exp(−ΔG*/kT) exp(−Q/kT)` with the
  spherical-cap barrier `ΔG* = 16π γ³ / (3 Δg_v²)`.  Driving force
  is taken from `tpt-mat-calphad`.
- **KWN (Kampmann–Wagner-Numerical)** — discretised size-class
  solver (1-D mass balance on the precipitate size distribution)
  with upwind advection and a single nucleation seed at the
  smallest class.
- **LSW coarsening** — `R̄³ − R̄_0³ = K t` with `K = 8 γ D c_eq / (9 V_m)`.
- **Strengthening hand-off** — Orowan bypass
  `Δτ = 0.4 μ b / (π √(1 − ν) L)` and Friedel shear `Δτ = 0.13 μ b / L`,
  passed to `tpt-mat-hardening` as additive yield-strength bumps.

## 3. Verification

- KWN conserves solute mass within solver tolerance.
- Late-stage size-distribution slope → LSW `t^{1/3}`.
- Strengthening increments vanish at `f = 0` and grow monotonically
  with `f` and `1/r̄`.

## 4. References

- Turnbull, D., & Fisher, J. C. (1949).  *J. Chem. Phys.* 17, 71.
- Kampmann, R., & Wagner, R. (1984).  Decomposition of alloys:
  the early stages.  *Acta/Scripta Met.* Conf. Proc.
- Lifshitz, I. M., & Slyozov, V. V. (1961).  *J. Phys. Chem.
  Solids* 19, 35.
- Wagner, C. (1961).  *Ber. Bunsenges. Phys. Chem.* 65, 581.
- Orowan, E. (1948).  *Symposium on Internal Stresses in Metals
  and Alloys,* Institute of Metals, London, 451.