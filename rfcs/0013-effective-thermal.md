# RFC 0013 — Effective Thermal & Transport-Property Homogenization

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-thermal`.

---

## 1. Summary

The `tpt-mat-thermal` crate bundles the effective-property
homogenization schemes that share the same symmetric-positive-
definite transport-tensor math as the elastic stiffness:
thermal conductivity `k`, electrical conductivity `σ`, ionic
diffusivity `D`, and linear coefficient of thermal expansion
(CTE) `α`.

## 2. Models

- **Effective conductivity** — Voigt / Reuss / VRH arithmetic and
  harmonic averages, Hashin–Shtrikman two-phase scalar bounds,
  Maxwell–Garnett.
- **Effective CTE** — Turner (bulk-modulus weighted),
  Kerner (matrix + spherical inclusion), Rosen–Hashin upper /
  lower bounds.
- **Effective specific heat** — mass-weighted rule of mixtures.
- **Effective diffusivity** — Bruggeman tortuosity-corrected
  relation for porous media.
- **Interface thermal resistance** — Kapitza laminate correction.

## 3. Generalization from `tpt-mat-homogenization`

The 6×6 / scalar bound machinery in `tpt-mat-homogenization` is
specialised to a generic symmetric-positive transport tensor
(k / σ / D) so the same code path serves mechanical and thermal
properties.

## 4. Verification

- Voigt ≥ Reuss.
- Hashin–Shtrikman lower ≤ upper.
- Maxwell–Garnett bracketed by Voigt / Reuss.
- CTE bounds collapse at zero phase contrast.
- Kapitza correction returns zero at zero interface resistance.

## 5. References

- Hashin, Z., & Shtrikman, S. (1962).  *J. Appl. Phys.* 33,
  3125–3131.
- Turner, P. S. (1946).  *J. Res. NBS* 37, 239–247.
- Rosen, B. W., & Hashin, Z. (1970).  *Int. J. Eng. Sci.* 8,
  157–173.
- Bruggeman, D. A. G. (1935).  *Ann. Phys.* 24, 636-664.
- Kapitza, P. L. (1941).  *J. Phys. USSR* 4, 181.