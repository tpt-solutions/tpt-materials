# RFC 0014 — Dislocation-Density-Based Hardening

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-dislocation`.

---

## 1. Summary

Dislocation-density-based hardening provides a physically-grounded
alternative to the phenomenological Voce / PowerLaw hardening
laws in `tpt-mat-hardening`.  This RFC documents the
`tpt-mat-dislocation` crate that exposes the standard Kocks–Mecking
evolution, Taylor stress, Nye-tensor → GND density, and
Armstrong–Frederick back-stress.

## 2. Models

- **DislocationDensityState** — per-slip-system `ρ_SSD`,
  `ρ_GND`, forest density.
- **KocksMeckingEvolution** — forward-Euler step
  `ρ_{n+1} = ρ_n + (k₁ √ρ_f − k₂ ρ_n) Δγ`; saturation density
  `ρ_sat = (k₁ / k₂)²`.
- **Taylor stress** — `τ = α μ b √ρ`.
- **Back stress** — Armstrong–Frederick kinematic term driven by
  the GND density.
- **gnd_from_curvature** — Nye-tensor trace `||α||` to scalar
  GND density from a lattice-curvature field.

## 3. Hardening integration

A new `HardeningLaw::DislocationDensity` variant is the planned
hand-off target into `tpt-mat-hardening` /
`tpt-mat-crystal-plasticity`.  That integration is queued for a
follow-up PR (the present crate ships the underlying evolution
laws and a verification test for single-slip Voce-like
saturation).

## 4. Verification

- `ρ` stays non-negative through any stress history.
- Late-stage saturation density matches `(k₁ / k₂)²` within
  10 %.
- Single-slip response reproduces Voce-like saturation.

## 5. References

- Kocks, U. F. (1976).  *J. Eng. Mater. Technol.* 98, 76–85.
- Mecking, H., & Kocks, U. F. (1981).  *Acta Metall.* 29,
  1865–1875.
- Taylor, G. I. (1934).  *Proc. R. Soc. A* 145, 362–387.
- Nye, J. F. (1953).  *Acta Metall.* 1, 153–162.
- Armstrong, P. J., & Frederick, C. O. (1966).  CEGB Report
  RD/B/N731.