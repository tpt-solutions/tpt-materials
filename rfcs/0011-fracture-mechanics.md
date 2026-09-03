# RFC 0011 — Fracture Mechanics: LEFM + Cohesive Zones + Phase-Field Fracture

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-fracture`.

---

## 1. Summary

This RFC covers the classical analytical fracture-mechanics crate
`tpt-mat-fracture`, which is not in the spec's `spec.txt` crate
list but is a standard micro-scale topic required by any
mechanical-engineering application of the framework.

## 2. Models

- **Stress intensity factors** — `K_I`, `K_II`, `K_III` for
  centre / edge / penny-shaped crack geometries.
- **Energy release rate** — Irwin `G = K_I² / E'` (plane
  stress / plane strain) and J-integral pass-through.
- **Cohesive-zone models** — bilinear traction–separation with
  Benzeggagh–Kenane mixed-mode decomposition.
- **Phase-field fracture** — AT1 and AT2 degradation functions,
  crack-surface dissipation functional `∫ (G_c / c_w) [w(d) / l_0
  + l_0 |∇d|²] dV`.
- **Master curve** — ASTM E1921 reference temperature `T_0` for
  the ductile-to-brittle transition.

## 3. Verification

- Irwn `G → 0` as `K → 0`.
- AT1: `g(0) = 1`, `g(1) = 0`, monotone decreasing.
- Bilinear cohesive law area equals `G_c` (numerical integration
  test).
- ASTM E1921 master curve monotonicity.

## 4. References

- Griffith, A. A. (1921).  *Phil. Trans. R. Soc. A* 221, 163–198.
- Irwin, G. R. (1957).  *J. Appl. Mech.* 24, 361–364.
- Benzeggagh, M. L., & Kenane, M. (1996).  *Compos. Sci. Technol.*
  56, 439–449.
- Bourdin, B., Francfort, G., & Marigo, J.-J. (2000).  *JMPS* 48,
  797–826.
- Miehe, C., Welschinger, F., & Hofacker, M. (2010).  *IJNME* 83,
  1273–1311.
- ASTM E1921 (latest).  *Standard Test Method for Determination
  of Reference Temperature, T₀, for Ferritic Steels in the
  Transition Range.*