# RFC 0009 — Soft-Matter Constitutive Models: Polymer, Hydrogel, Hydrogen Storage

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-polymer`, `tpt-mat-hydrogel`,
`tpt-mat-hydrogen-storage`.

---

## 1. Summary

The spec §5 Domain 7 / 8 workplan did not schedule the
soft-matter crates, but the spec §4 workspace layout names
`tpt-mat-polymer`, `tpt-mat-hydrogel` and `tpt-mat-hydrogen-storage`.
This RFC documents the implementation in three crates:

- **`tpt-mat-polymer`** — Arruda–Boyce 8-chain rubber elasticity,
  Worm-Like-Chain (Marko–Siggia) force–extension, Freely-Jointed-Chain
  inverse-Langevin.
- **`tpt-mat-hydrogel`** — Flory–Rehner equilibrium swelling,
  poroelastic Fickian slab uptake.
- **`tpt-mat-hydrogen-storage`** — van 't Hoff PCT isotherms for
  metal, chemical, and porous-material hydrides.

## 2. `tpt-mat-polymer`

The Arruda–Boyce (1993) 8-chain model reduces to neo-Hookean at
small stretch and diverges at the finite-extensibility limit via the
inverse-Langevin function.  The Worm-Like-Chain Marko–Siggia
approximation handles force–extension behaviour for semi-flexible
chains.

The verification test `arruda_boyce_neo_hookean_small_stretch` is
satisfied: `arruda_boyce_stress(λ) → 3 n k T (λ − 1/λ)` for `λ → 1`.

## 3. `tpt-mat-hydrogel`

Flory–Rehner (1943) balances the mixing osmotic pressure against
the elastic osmotic pressure of the cross-linked network.  The
[`equilibrium_swelling_ratio`](tpt-mat-hydrogel::equilibrium_swelling_ratio)
helper root-finds over the polymer-volume fraction.  Poroelastic
Fickian slab uptake reuses the `tpt-science` grid.

## 4. `tpt-mat-hydrogen-storage`

[`pct_isotherm`](tpt-mat-hydrogen-storage::pct_isotherm) builds the
pressure–composition–temperature curve for a `HydrogenStorageMaterial`
with van 't Hoff equilibrium and the standard plateau-pressure
plateau.  Output is consumed by `tpt-energy` and `tpt-transport`.

## 5. References

- Arruda, E. M., & Boyce, M. C. (1993).  *J. Mech. Phys. Solids*
  41(2), 389–412.
- Marko, J. F., & Siggia, E. D. (1995).  *Macromolecules* 28,
  8759–8770.
- Flory, P. J., & Rehner, J. (1943).  *J. Chem. Phys.* 11, 521–526.
- van 't Hoff, M. J. H. (1884).  *Etudes de dynamique chimique.*  F.
  Muller.