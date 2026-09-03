# RFC 0007 — Welding & Heat-Treatment Transformation Pipeline

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-welding`, `tpt-mat-heat-treatment`,
`tpt-mat-phase-transform`.

---

## 1. Summary

This RFC ships the welding-process crate
[`tpt-mat-welding`](tpt-mat-welding) that composes:

- Rosenthal-line-source HAZ width & peak-temperature profile (Phase 7
  welding foundation),
- an empirical cooling-rate model at the coarse-grained HAZ
  (Kou 2003),
- a Hall–Petch-like CG-HAZ grain-size prediction,
- a sub-zone classification (base-metal, sub-critical,
  inter-critical, fine-grained, coarse-grained, martensite, bainite,
  pearlite) driven by the Avrami / Koistinen–Marburger kinetics from
  [`tpt-mat-phase-transform`](tpt-mat-phase-transform).

It sits between `tpt-mat-additive` (process → thermal history) and
`tpt-mat-heat-treatment` (microstructure → hardness).

## 2. HAZ width

[`heat_affected_zone`] bisects the perpendicular distance `y` at
which the Rosenthal peak temperature falls below the supplied
threshold (typically `A_1`).  [`peak_temperature_profile`] samples
the temperature-vs-distance curve on a log-spaced grid for plotting.

## 3. CG-HAZ microstructure

[`predict_haz_microstructure`] returns the eight-element phase-fraction
vector (matching the
[`MicrostructurePhase`](tpt-mat-welding::MicrostructurePhase) order)
plus the CG-HAZ grain size and cooling rate.  Phase fractions are
driven by:

- An Avrami isothermal hold at the inter-critical temperature
  (Bainite kinetics).
- A Koistinen–Marburger martensite kick-in when `T < M_s` and the
  cooling rate exceeds the bainite/pearlite bypass threshold.

## 4. Example

[`examples/welding-haz`] runs a GMAW on AISI 4140 steel: HAZ width,
peak-temperature profile, CG-HAZ grain-coarsening and phase-fraction
prediction.

## 5. References

- Kou, S. (2003).  *Welding Metallurgy.*  2nd ed., Wiley.
- Easterling, K. (1992).  *Introduction to the Physical Metallurgy
  of Welding.*  2nd ed., Butterworth-Heinemann.
- Avrami, M. (1939).  *J. Chem. Phys.* 7, 1103.
- Koistinen, D. P., & Marburger, R. E. (1959).  *Acta Metall.* 7, 59.