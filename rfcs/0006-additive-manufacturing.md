# RFC 0006 — Additive-Manufacturing Microstructure & Residual-Stress

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-additive`.

---

## 1. Summary

This RFC introduces [`tpt-mat-additive`], a closed-form analytical
predictor of thermal history, microstructure and residual stress for
the three dominant additive-manufacturing processes:

- Laser powder-bed fusion (LPBF)
- Directed-energy deposition (DED)
- Electron-beam melting (EBM)

## 2. Thermal history

The moving-point-source Rosenthal (1941) solution drives the
temperature field; the Eagar–Tsai (1983) Gaussian-beam correction is
included as [`GaussianBeam`] for finite-width sources.  The
[`thermal_history`] driver samples the time–temperature history at a
fixed perpendicular distance from the scan track; the cooling rate at
the peak is used downstream by [`predict_microstructure`].

## 3. Microstructure (Hunt CET map)

The Hunt (1984) columnar-to-equiaxed transition parameterises the
solidification microstructure by the thermal gradient `G` and the
growth rate `R`:

```text
columnar :  G · R² ≥ G_CET
equiaxed :  G · R² < G_CET
```

The grain-size prediction [`grain_size_from_cooling_rate`] uses the
empirical Hall–Petch-like calibration `d ∝ ε̇^(-0.2)`, with
`ε̇ = cooling_rate / ΔT_solidus`.

## 4. Residual stress

The 1-D eigenstrain model [`residual_stress`] captures the thermal-
contraction field in a build attached to a substrate: linear CTE
`α`, solidus temperature `T_solidus`, room temperature `T_room`,
elastic constants `E, ν`, and a plastic-relaxation fraction `f`:

```text
σ = E · α · (T_solidus − T_room) · (1 − f) / (1 − ν)
```

## 5. Example

[`examples/additive-manufacturing-microstructure`] drives the LPBF
Ti-6Al-4V scan: linear / volumetric energy density, Rosenthal peak
at `d = 100 µm`, Hunt CET map, and a 1-D residual-stress field along
the build direction.

## 6. References

- Rosenthal, D. (1941).  *Weld. J.* 20, 220s–234s.
- Eagar, T. W., & Tsai, N.-S. (1983).  *Weld. J.* 62(12), 346s–355s.
- Hunt, J. D. (1984).  *Mater. Sci. Eng.* 65(1), 75–83.
- Kou, S. (2003).  *Welding Metallurgy.*  2nd ed., Wiley.