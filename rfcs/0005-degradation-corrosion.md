# RFC 0005 — Degradation & Failure: FIP-Based Microstructural Fatigue and Butler–Volmer Corrosion

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-fatigue-micro`, `tpt-mat-corrosion`.

---

## 1. Summary

This RFC adds two crates to the Phase-6 Degradation cluster:

- **`tpt-mat-fatigue-micro`** — fatigue-indicator-parameter (FIP)
  analysis on a CP-FEM RVE, driving crack-initiation prediction.
- **`tpt-mat-corrosion`** — Butler–Volmer / Tafel electrochemistry,
  mixed-potential corrosion rate, polarization curves.

These complement the classical `tpt-mat-fatigue`, `tpt-mat-damage` and
`tpt-mat-creep` crates already in the repo (Phase 6 / Phase 10).

## 2. `tpt-mat-fatigue-micro`

### 2.1 FIP criteria

Four classical critical-plane FIPs are implemented:

| Criterion | Formula | Inputs |
|---|---|---|
| Findley (1957) | `F = max_θ (Δτ/2 + k σ_n_max)` | `k` (stress multiaxiality sensitivity) |
| Fatemi–Socie (1988) | `F = (Δγ/2) (1 + k ⟨σ_n_max⟩/σ_y)` | `k`, `σ_y` |
| Smith–Watson–Topper (1970) | `F = σ_max · (Δε/2)` | none |
| Tanaka–Mura (1981) | `F = max_slip γ_acc^α / γ_c` | `γ_c` (critical shear) |

### 2.2 Pipeline

1. Run a cycle-by-cycle CP-FEM simulation on the RVE
   ([`tpt-mat-crystal-plasticity::CpFemSolver`]).
2. At each load step extract
   [`tpt-mat-crystal-plasticity::CpFemResult`] (stresses + accumulated
   shear per slip system per grain).
3. Compute the FIP field with
   [`tpt-mat-fatigue-micro::fatigue_indicator_parameter`].
4. Drive
   [`tpt-mat-fatigue-micro::predict_crack_initiation`] to identify
   the critical grain and estimate cycles-to-initiation with a
   Coffin–Manson / Tanaka–Mura law.

### 2.3 Dependencies

- `tpt-math-linalg-fixed` — 3×3 / 6×6 linear algebra for the
  critical-plane search.
- `tpt-mat-crystal-plasticity` — `CpFemResult` interface.
- `tpt-mat-rve` — `Rve` data model.

## 3. `tpt-mat-corrosion`

### 3.1 Butler–Volmer kinetics

A single half-reaction is described by
[`ElectrodeKinetics`](tpt-mat-corrosion::ElectrodeKinetics):

```text
i(E) = i_0 · ( exp((E − E_eq) / b_a) − exp(−(E − E_eq) / b_c) )
```

with `b = RT / (α n F)` the Tafel slope.  The
[`butler_volmer_current_density`](tpt-mat-corrosion::butler_volmer_current_density)
helper evaluates the net current density on the full
charge-transfer regime; [`tafel_anodic`](tpt-mat-corrosion::tafel_anodic) /
[`tafel_cathodic`](tpt-mat-corrosion::tafel_cathodic) are the
high-overpotential asymptotes.

### 3.2 Mixed-potential theory (Wagner–Traud 1938)

A [`CorrosionModel`](tpt-mat-corrosion::CorrosionModel) couples an
anode + cathode + electrolyte.  [`corrosion_rate`](tpt-mat-corrosion::corrosion_rate)
solves `i_a(E) = |i_c(E)|` by bisection (Wagner–Traud intersection)
and reports the penetration rate (mm/yr) and mass-loss rate
(g/m²·day) from the Faraday law.

### 3.3 Polarization curves

[`polarization_curve`](tpt-mat-corrosion::polarization_curve) samples
the anodic / cathodic / net branch over a potential range — the
output is consumed directly by `tpt-medical` and `tpt-transport` for
implant and pipe-corrosion lifetime predictions.

## 4. References

- Findley, W. N. (1957).  *Proc. ASTM* 57, 880–886.
- Fatemi, A., & Socie, D. F. (1988).  *Fatigue Fract. Eng. Mater.
  Struct.* 11(3), 149–165.
- Smith, K. N., Watson, P., & Topper, T. H. (1970).  *J. Mater.*
  5(4), 767–778.
- Tanaka, K. & Mura, T. (1981).  *J. Appl. Mech.* 48(1), 97–103.
- Wagner, C. & Traud, W. (1938).  *Z. Elektrochem.* 44(7), 391–402.
- Bard, A. J., & Faulkner, L. R. (2001).  *Electrochemical Methods.*
  2nd ed., Wiley.